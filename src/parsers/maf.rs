use std::fmt::{Debug, Display};
use std::io::BufRead;

// Create an iterator from a bufreader
pub fn maf_parser(file: std::fs::File) -> MafParser {
    let reader = std::io::BufReader::new(file);
    let lines = reader.lines();

    MafParser {
        lines,
        current_block: Vec::new(),
        block_start_line: None,
    }
}

pub struct MafParser {
    lines: std::io::Lines<std::io::BufReader<std::fs::File>>,
    current_block: Vec<MafLine>,
    block_start_line: Option<MafLine>,
}

impl Iterator for MafParser {
    type Item = Vec<MafLine>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(line) = self.lines.next() {
            let line = line.unwrap();

            // Match on the first character
            let x = parse_maf_line(&line);

            if let MafLine::BlankLine = x {
                if self.current_block.len() > 0 {
                    let blank = Vec::with_capacity(self.current_block.len());
                    let block = std::mem::replace(&mut self.current_block, blank);
                    return Some(block);
                } else {
                    continue;
                }
            } else {
                self.current_block.push(x);
            }
        }
        return None;
    }
}

// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
#[derive(Clone)]
pub enum MafLine {
    Comment(String),
    AlignmentBlockLine(String),
    SequenceLine(String, String, u64, u64, Strand, u64, String),
    BlankLine,
}

impl Display for MafLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MafLine::Comment(x) => {
                write!(f, "#{}", x)
            }
            MafLine::AlignmentBlockLine(x) => {
                write!(f, "a{}", x)
            }
            MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) => {
                write!(
                    f,
                    "s {} {} {} {} {} {} {}",
                    species,
                    seqid,
                    start,
                    length,
                    match strand {
                        Strand::Plus => "+",
                        Strand::Minus => "-",
                    },
                    src_size,
                    text
                )
            }
            MafLine::BlankLine => {
                write!(f, "")
            }
        }
    }
}

impl Debug for MafLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MafLine::Comment(x) => {
                write!(f, "Comment: {}", x)
            }
            MafLine::AlignmentBlockLine(x) => {
                write!(f, "AlignmentBlockLine {}", x)
            }
            MafLine::SequenceLine(species, seqid, start, length, strand, src_size, _text) => {
                write!(
                    f,
                    "SequenceLine {} {} {} {} {} {}",
                    species,
                    seqid,
                    start,
                    length,
                    match strand {
                        Strand::Plus => "+",
                        Strand::Minus => "-",
                    },
                    src_size,
                )
            }
            MafLine::BlankLine => {
                write!(f, "BlankLine")
            }
        }
    }
}

impl MafLine {

    pub fn is_seqline(&self) -> bool {
        match self {
            MafLine::SequenceLine(_, _, _, _, _, _, _) => true,
            _ => false,
        }
    }

    pub fn fasta_out(&self) -> String {
        match self {
            MafLine::Comment(_x) => {
                panic!("Cannot convert comment to fasta")
            }
            MafLine::AlignmentBlockLine(_x) => {
                panic!("Cannot convert alignment block to fasta")
            }
            MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) => {
                let mut out = String::new();
                out.push_str(&format!(
                    ">{}:{}-{} {} {}",
                    seqid,
                    start,
                    start + length,
                    match strand {
                        Strand::Plus => "+",
                        Strand::Minus => "-",
                    },
                    src_size,
                ));
                out.push_str("\n");
                out.push_str(text);
                out.push_str("\n");
                return out;
            }
            MafLine::BlankLine => {
                panic!("Cannot convert blank line to fasta")
            }
        }
    }
}

pub enum Strand {
    Plus,
    Minus,
}


fn parse_maf_line(line: &str) -> MafLine {
    // Match on the first character
    match line.chars().next() {
        Some('#') => {
            // Remove first character
            return MafLine::Comment(line[1..].to_string());
        }
        Some('a') => {
            return MafLine::AlignmentBlockLine(line[1..].to_string());
        }
        Some('s') => {
            // Split by tabs and convert to:
            // String, u64, u64, Strand, u64, String
            let mut split = line.split_whitespace();
            let _ = split.next(); // Remove first character
            let seqid = split.next().unwrap().to_string();

            // Seqid is in the format species.chromosome
            let species = seqid.split('.').next().unwrap().to_string();
            let seqid = seqid.split('.').nth(1).unwrap().to_string();

            let start = split.next().unwrap().parse::<u64>().unwrap();
            let length = split.next().unwrap().parse::<u64>().unwrap();
            let strand = match split.next().unwrap() {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                _ => panic!("Invalid strand"),
            };
            let src_size = split.next().unwrap().parse::<u64>().unwrap();
            let text = split.next().unwrap().to_string();

            return MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text);
        }
        None => {
            return MafLine::BlankLine;
        }
        _ => {
            return MafLine::BlankLine;
        }
    }
}