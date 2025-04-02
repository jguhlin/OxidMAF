use flate2::bufread::GzDecoder;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Error as IoError, Lines, Read, Seek, SeekFrom};

use super::ContainsResult;

/// Represents one index entry for a block.
#[derive(Debug)]
pub struct TaiEntry {
    /// The uncompressed block start coordinate (second column)
    block_start: u64,
    /// The file offset (third column). For bgzipped files, this is the compressed offset.
    offset: u64,
}

/// The TAI index, mapping a contig name to a vector of index entries.
#[derive(Debug)]
pub struct TaiIndex {
    entries: HashMap<String, Vec<TaiEntry>>,
}

impl TaiIndex {
    /// Load a TAI index from the given file path.
    ///
    /// The file is expected to be tab-delimited with three columns.
    /// The first column holds the contig name or "*" to indicate the same contig as the previous line.
    /// The second column is the block start (uncompressed coordinate).
    /// The third column is the file offset (for bgzipped files, the compressed offset).
    pub fn from_file(path: &str) -> Result<Self, std::io::Error> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut entries: HashMap<String, Vec<TaiEntry>> = HashMap::new();
        let mut current_contig = String::new();

        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue;
            }

            // Update the current contig if the first field is not "*"
            if fields[0] != "*" {
                // Split after the first dot to get the contig name.
                let parts: Vec<&str> = fields[0].splitn(2, '.').collect();
                if parts.len() > 1 {
                    current_contig = parts[1].to_string();
                } else {
                    current_contig = fields[0].to_string();
                }
            }

            // Parse the block start and offset
            let block_start = fields[1].parse::<u64>().unwrap_or(0);
            let offset = fields[2].parse::<u64>().unwrap_or(0);
            let entry = TaiEntry {
                block_start,
                offset,
            };

            entries
                .entry(current_contig.clone())
                .or_insert_with(Vec::new)
                .push(entry);
        }
        Ok(TaiIndex { entries })
    }

    /// Given a contig and a target uncompressed coordinate, find the block entry with the highest block_start that is <= target.
    /// Returns a tuple: (block_start, file_offset).
    pub fn get_seek_info(&self, contig: &str, pos: u64) -> Option<(u64, u64)> {
        self.entries.get(contig).and_then(|entries| {
            // Since entries are inserted in order, we can iterate until block_start exceeds the target.
            let mut candidate = None;
            for entry in entries {
                if entry.block_start <= pos {
                    candidate = Some(entry);
                } else {
                    break;
                }
            }
            candidate.map(|e| (e.block_start, e.offset))
        })
    }
}

#[derive(Debug)]
pub struct TafHeader {
    pub tags: HashMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct TafColumn {
    pub raw_bases: String,
    pub alleles: Vec<char>,
    pub coordinates: Vec<CoordinateOp>,
    pub tags: HashMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct Coordinate {
    pub species: String,
    pub chrom: String,
    pub offset: u64,
    pub strand: Strand,
    pub sequence_length: Option<u64>,
}

#[derive(Debug, Clone)]
pub enum CoordinateOp {
    Insertion {
        row: usize,
        coord: Option<Coordinate>,
    },
    Deletion {
        row: usize,
    },
    Substitution {
        row: usize,
        coord: Option<Coordinate>,
    },
    Gap {
        row: usize,
        gap_length: usize,
    },
    GapString {
        row: usize,
        gap_string: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

// ---
// Now, the parser is generic over any type that implements BufRead.

pub struct TafParser {
    pub header: TafHeader,
    pub run_length_encode: bool,
    inner: BufReader<File>,
    inner_buf: BufReader<Cursor<Vec<u8>>>,
    line: String,
}

impl TafParser {
    /// Create a new parser from any reader that implements BufRead.
    pub fn from_file(path: &str) -> Result<Self, IoError> {
        let reader = File::open(path).unwrap();
        let mut reader = BufReader::new(reader);

        let mut inner_buf = Vec::new();
        // Decode the first Gz Block
        let mut gz_decoder = GzDecoder::new(reader);
        gz_decoder
            .read_to_end(&mut inner_buf)
            .map_err(|e| IoError::new(std::io::ErrorKind::Other, e))?;

        // Create a new reader from the decompressed data
        let mut reader = gz_decoder.into_inner();

        let mut inner_buf = BufReader::new(std::io::Cursor::new(inner_buf));

        // The header must be the first line.
        // Get the first line
        let mut header_line = String::new();
        inner_buf.read_line(&mut header_line)?;
        if !header_line.starts_with("#taf") {
            return Err(IoError::new(
                std::io::ErrorKind::InvalidData,
                "Missing #taf header",
            ));
        }
        let header = parse_taf_header(&header_line);
        let run_length_encode = match header.tags.get("run_length_encode_bases") {
            Some(val) if val == "1" => true,
            _ => false,
        };
        Ok(TafParser {
            header,
            run_length_encode,
            inner: reader,
            line: String::new(),
            inner_buf,
        })
    }

    pub fn seek_to(&mut self, pos: u64, block_offset: u64) {
        // Seek to the specified position in the file.
        self.inner.seek(SeekFrom::Start(pos)).unwrap();

        // Decode the next block
        let mut gz_decoder = GzDecoder::new(&mut self.inner);
        let mut inner_buf = Vec::new();
        gz_decoder
            .read_to_end(&mut inner_buf)
            .expect("Failed to read gzipped data");
        // Create a new reader from the decompressed data
        self.inner_buf = BufReader::new(Cursor::new(inner_buf));
        self.inner_buf
            .seek(SeekFrom::Start(block_offset))
            .expect("Failed to seek to block offset");
    }

    fn read_next_block(&mut self) -> Result<(), IoError> {
        // Read the next block of data from the inner buffer.
        let mut gz_decoder = GzDecoder::new(&mut self.inner);
        let mut inner_buf = Vec::new();
        gz_decoder
            .read_to_end(&mut inner_buf)
            .map_err(|e| IoError::new(std::io::ErrorKind::Other, e))?;
        // Create a new reader from the decompressed data
        self.inner_buf = BufReader::new(Cursor::new(inner_buf));
        Ok(())
    }
}

impl Iterator for TafParser {
    type Item = Result<TafColumn, String>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read lines until we find a valid TAF column or EOF.
        loop {
            self.line.clear();
            if self.inner_buf.read_line(&mut self.line).is_err() {
                // Try to read the next block
                match self.read_next_block() {
                    Ok(_) => continue,
                    Err(e) => return Some(Err(e.to_string())),
                }
            }

            let line = self.line.trim();
            if line.is_empty() || (line.starts_with('#') && !line.starts_with("#taf")) {
                continue;
            }
            return Some(parse_taf_column(line, self.run_length_encode));
        }
    }
}

// ---
// The rest of your parser (header parsing, decode_bases, coordinate parsing, etc.) remains unchanged.

fn parse_taf_header(line: &str) -> TafHeader {
    let without_hash = line.trim_start_matches('#');
    let mut parts = without_hash.split_whitespace();
    let _ = parts.next(); // Skip "taf"
    let mut tags = HashMap::new();
    for token in parts {
        if let Some((key, value)) = token.split_once(':') {
            tags.insert(key.to_string(), value.to_string());
        }
    }
    TafHeader { tags }
}

fn parse_taf_column(line: &str, run_length: bool) -> Result<TafColumn, String> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.is_empty() {
        return Err("Empty column line".to_string());
    }

    let mut i = 0;
    let raw_bases: String;
    if run_length {
        // Accumulate tokens for run-length encoded bases until we hit ";" or a tag token.
        let mut bases_tokens = Vec::new();
        while i < tokens.len() && tokens[i] != ";" && !tokens[i].starts_with('@') {
            bases_tokens.push(tokens[i]);
            i += 1;
        }
        raw_bases = bases_tokens.join("");
    } else {
        // In non-run-length mode, the first token is the entire bases string.
        raw_bases = tokens[0].to_string();
        i = 1;
    }
    let alleles = decode_bases(&raw_bases, run_length)?;

    // Now, parse optional coordinates and tag tokens starting at the current token index.
    let mut coord_tokens = Vec::new();
    let mut tag_tokens = Vec::new();
    while i < tokens.len() {
        let token = tokens[i];
        if token == ";" {
            i += 1;
            while i < tokens.len() && !tokens[i].starts_with('@') {
                coord_tokens.push(tokens[i]);
                i += 1;
            }
        } else if token.starts_with('@') {
            // Remove the leading '@'
            let tag_str = &token[1..];
            tag_tokens.push(tag_str);
            i += 1;
            while i < tokens.len() {
                tag_tokens.push(tokens[i]);
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    let coordinates = parse_coordinate_ops(&coord_tokens)?;
    let tags = parse_tags(&tag_tokens);
    Ok(TafColumn {
        raw_bases,
        alleles,
        coordinates,
        tags,
    })
}

fn decode_bases(s: &str, run_length: bool) -> Result<Vec<char>, String> {
    if !run_length {
        Ok(s.chars().collect())
    } else {
        let mut alleles = Vec::new();
        let mut chars = s.chars().peekable();
        while let Some(c) = chars.next() {
            let mut num_str = String::new();
            while let Some(&digit) = chars.peek() {
                if digit.is_ascii_digit() {
                    num_str.push(digit);
                    chars.next();
                } else {
                    break;
                }
            }
            if num_str.is_empty() {
                return Err(format!("Expected count after base '{}'", c));
            }
            let count: usize = num_str
                .parse()
                .map_err(|e: std::num::ParseIntError| e.to_string())?;
            for _ in 0..count {
                alleles.push(c);
            }
        }
        Ok(alleles)
    }
}

fn parse_coordinate_ops(tokens: &[&str]) -> Result<Vec<CoordinateOp>, String> {
    let mut ops = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let op_token = tokens[i];
        i += 1;
        match op_token {
            "i" | "s" => {
                if i + 3 >= tokens.len() {
                    return Err(format!("Not enough tokens for op '{}'", op_token));
                }
                let row = tokens[i]
                    .parse::<usize>()
                    .map_err(|e| format!("Parsing row: {}", e))?;
                i += 1;
                // Instead of taking the token as the full sequence name,
                // split it into species and chromosome.
                let full_seq_name = tokens[i].to_string();
                i += 1;
                let mut parts = full_seq_name.splitn(2, '.');
                let species = parts.next().unwrap().to_string();
                let chrom = parts.next().unwrap_or("").to_string();
                let offset = tokens[i]
                    .parse::<u64>()
                    .map_err(|e| format!("Parsing offset: {}", e))?;
                i += 1;
                let strand = match tokens[i] {
                    "+" => Strand::Plus,
                    "-" => Strand::Minus,
                    _ => return Err("Invalid strand; expected '+' or '-'".to_string()),
                };
                i += 1;
                let mut sequence_length = None;
                if i < tokens.len() && !["i", "d", "s", "g", "G"].contains(&tokens[i]) {
                    if let Ok(len) = tokens[i].parse::<u64>() {
                        sequence_length = Some(len);
                        i += 1;
                    }
                }
                let coord = Coordinate {
                    species,
                    chrom,
                    offset,
                    strand,
                    sequence_length,
                };
                if op_token == "i" {
                    ops.push(CoordinateOp::Insertion {
                        row,
                        coord: Some(coord),
                    });
                } else {
                    ops.push(CoordinateOp::Substitution {
                        row,
                        coord: Some(coord),
                    });
                }
            }
            "d" => {
                if i >= tokens.len() {
                    return Err("Not enough tokens for deletion op".to_string());
                }
                let row = tokens[i]
                    .parse::<usize>()
                    .map_err(|e| format!("Parsing row in deletion: {}", e))?;
                i += 1;
                ops.push(CoordinateOp::Deletion { row });
            }
            "g" => {
                if i + 1 >= tokens.len() {
                    return Err("Not enough tokens for gap op".to_string());
                }
                let row = tokens[i]
                    .parse::<usize>()
                    .map_err(|e| format!("Parsing row in gap: {}", e))?;
                i += 1;
                let gap_length = tokens[i]
                    .parse::<usize>()
                    .map_err(|e| format!("Parsing gap_length: {}", e))?;
                i += 1;
                ops.push(CoordinateOp::Gap { row, gap_length });
            }
            "G" => {
                if i + 1 >= tokens.len() {
                    return Err("Not enough tokens for gap string op".to_string());
                }
                let row = tokens[i]
                    .parse::<usize>()
                    .map_err(|e| format!("Parsing row in gap string: {}", e))?;
                i += 1;
                let gap_string = tokens[i].to_string();
                i += 1;
                ops.push(CoordinateOp::GapString { row, gap_string });
            }
            unknown => {
                return Err(format!("Unknown coordinate op: {}", unknown));
            }
        }
    }
    Ok(ops)
}

fn parse_tags(tokens: &[&str]) -> HashMap<String, String> {
    let mut tags = HashMap::new();
    for token in tokens {
        if let Some((key, value)) = token.split_once(':') {
            tags.insert(key.to_string(), value.to_string());
        }
    }
    tags
}

/// An alignment iterator that wraps the TafParser and maintains two mappings:
/// one for species names and one for current coordinates (updated by coordinate ops).
pub struct TafAlignmentIterator<'taf> {
    parser: &'taf mut TafParser,
    /// Mapping from row index to species name (if known).
    species_map: Vec<Option<String>>,
    /// Mapping from row index to the current coordinate for that row.
    current_coords: Vec<Option<Coordinate>>,
    col_index: usize,
}

/// A structure holding a column plus the current species mapping, coordinates, and its index.
#[derive(Debug, Clone)]
pub struct TafAlignmentColumn {
    /// The parsed column.
    pub column: TafColumn,
    /// The current mapping from row indices to species.
    pub species_map: Vec<Option<String>>,
    /// The current mapping from row indices to coordinates.
    pub coords: Vec<Option<Coordinate>>,
    /// Column index (0-based) in the alignment.
    pub col_index: usize,
}

impl<'taf> TafAlignmentIterator<'taf> {
    /// Create a new alignment iterator from the given parser.
    pub fn new(parser: &'taf mut TafParser) -> Self {
        TafAlignmentIterator {
            parser,
            species_map: Vec::new(),
            current_coords: Vec::new(),
            col_index: 0,
        }
    }
}

impl<'taf> Iterator for TafAlignmentIterator<'taf> {
    type Item = Result<TafAlignmentColumn, String>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.parser.next() {
            Some(Ok(mut col)) => {
                // Ensure our mappings have as many entries as there are rows.
                if self.species_map.len() < col.alleles.len() {
                    self.species_map.resize(col.alleles.len(), None);
                    self.current_coords.resize(col.alleles.len(), None);
                }
                // Track which rows are updated in this column.
                let mut updated_rows: HashSet<usize> = HashSet::new();

                // Update mappings based on coordinate operations.
                for op in &col.coordinates {
                    match op {
                        CoordinateOp::Insertion { row, coord }
                        | CoordinateOp::Substitution { row, coord } => {
                            if let Some(c) = coord {
                                if *row < self.species_map.len() {
                                    self.species_map[*row] = Some(c.species.clone());
                                    self.current_coords[*row] = Some(c.clone());
                                    updated_rows.insert(*row);
                                }
                            }
                        }
                        CoordinateOp::Gap { row, gap_length } => {
                            if *row < self.current_coords.len() {
                                if let Some(ref mut c) = self.current_coords[*row] {
                                    c.offset += *gap_length as u64;
                                }
                                updated_rows.insert(*row); // treat gap as an update
                            }
                        }
                        CoordinateOp::GapString { row, .. } => {
                            // You can update similarly if needed.
                            updated_rows.insert(*row);
                        }
                        _ => {}
                    }
                }
                // For any row that was not updated in this column, advance its coordinate by 1 (if present).
                for i in 0..self.current_coords.len() {
                    if !updated_rows.contains(&i) {
                        if let Some(coord) = &mut self.current_coords[i] {
                            coord.offset += 1;
                        }
                    }
                }
                let result = TafAlignmentColumn {
                    column: col,
                    species_map: self.species_map.clone(),
                    coords: self.current_coords.clone(),
                    col_index: self.col_index,
                };
                self.col_index += 1;
                Some(Ok(result))
            }
            Some(Err(e)) => Some(Err(e)),
            None => None,
        }
    }
}

impl TafAlignmentColumn {
    /// Returns true if the reference coordinate (at ref_index) matches the given position.
    /// Here we assume:
    ///   - TAF coordinates are 0-based.
    ///   - VCF positions are 1-based.
    pub fn ref_matches_pos(&self, ref_index: usize, pos: u64) -> bool {
        if let Some(Some(coord)) = self.coords.get(ref_index) {
            // Compare the coordinate offset to (pos - 1)
            coord.offset == pos - 1
        } else {
            false
        }
    }

    /// Given a species name, return the allele (and its row) present in this column.
    pub fn allele_for_species(&self, species: &str) -> Option<(usize, char)> {
        for (i, s_opt) in self.species_map.iter().enumerate() {
            if let Some(s) = s_opt {
                if s == species && i < self.column.alleles.len() {
                    return Some((i, self.column.alleles[i]));
                }
            }
        }
        None
    }

    /// Return a coordinate for the given chromosome if present.
    pub fn ref_coord(&self, chrom: &str) -> Option<&Coordinate> {
        self.coords.iter().find_map(|c| {
            if let Some(coord) = c {
                if coord.chrom == chrom {
                    Some(coord)
                } else {
                    None
                }
            } else {
                None
            }
        })
    }

    /// Checks if the reference coordinate for `chrom` contains `pos`.
    /// If the coordinate has a known length, we use that; otherwise we assume length 1.
    pub fn contains_pos(&self, chrom: &str, pos: u64) -> ContainsResult {
        if let Some(coord) = self.ref_coord(chrom) {
            let len = coord.sequence_length.unwrap_or(1);
            let start = coord.offset;
            let end = start + len;
            if pos < start {
                ContainsResult::Before
            } else if pos >= end {
                ContainsResult::After
            } else {
                ContainsResult::True
            }
        } else {
            ContainsResult::WrongChrom
        }
    }
}

