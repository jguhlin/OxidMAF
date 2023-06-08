// use bevy_app::prelude::*;
// use bevy_app::App;
// use bevy_ecs::prelude::*;
// use bevy_tasks::TaskPool;
use clap::{Parser, Subcommand};

use std::fmt::Display;
use std::io::BufRead;
use std::io::{Seek, SeekFrom, Write};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(about = "Count Reference Gaps")]
    CountRefGaps { input: String },
    #[command(about = "Count Duplicate Reference Entries")]
    CountDupeRefs { input: String },
    #[command(about = "Remove Alignment Blocks with Duplicate Reference Entries")]
    RemoveDupeRefBlocks { input: String },
    #[command(about = "Split MAF File")]
    Split { input: String, output_path: String },
    #[command(about = "Process GERP Scores from .maf, .rates, and .elems files and export to .tsv")]
    ProcessGERP {
        maf: String,
        rates: String,
        elems: String,
        output: String,
    },
}

fn main() {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Split { input, output_path } => {
            println!("split was used")
        }
        Commands::CountRefGaps { input } => {
            count_ref_gaps(input);
        }
        Commands::CountDupeRefs { input } => {
            count_dupe_refs(input);
        }
        Commands::RemoveDupeRefBlocks { input } => {
            remove_dupe_ref_blocks(input);
        }
        Commands::ProcessGERP {
            maf,
            rates,
            elems,
            output,
        } => {
            process_gerp(maf, rates, elems, output);
        }
    }
}

// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
enum MafReadState {
    Comment,
    AlignmentBlockLine,
    SequenceLine,
    BlankLine,
}

pub enum Strand {
    Plus,
    Minus,
}

// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
pub enum MafLine {
    Comment(String),
    AlignmentBlockLine(String),
    SequenceLine(String, u64, u64, Strand, u64, String),
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
            MafLine::SequenceLine(seqid, start, length, strand, src_size, text) => {
                write!(
                    f,
                    "s {} {} {} {} {} {}",
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
            let start = split.next().unwrap().parse::<u64>().unwrap();
            let length = split.next().unwrap().parse::<u64>().unwrap();
            let strand = match split.next().unwrap() {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                _ => panic!("Invalid strand"),
            };
            let src_size = split.next().unwrap().parse::<u64>().unwrap();
            let text = split.next().unwrap().to_string();

            return MafLine::SequenceLine(seqid, start, length, strand, src_size, text);
        }
        None => {
            return MafLine::BlankLine;
        }
        _ => {
            return MafLine::BlankLine;
        }
    }
}

fn count_ref_gaps(input: &String) {
    // Open file
    let file = std::fs::File::open(input).unwrap();
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut count = 0;
    let mut state = MafReadState::BlankLine;
    let mut seqcount = 0;

    while let Some(line) = lines.next() {
        let line = line.unwrap();

        // Match on the first character
        match line.chars().next() {
            Some('#') => {
                state = MafReadState::Comment;
                continue;
            }
            Some('a') => {
                state = MafReadState::AlignmentBlockLine;
                seqcount = 0;
                continue;
            }
            Some('s') => {
                state = MafReadState::SequenceLine;
                seqcount += 1;
            }
            None => {
                state = MafReadState::BlankLine;
                continue;
            }
            _ => {
                state = MafReadState::BlankLine;
                continue;
            }
        }

        if seqcount == 1 {
            let mut gapcount = 0;
            for c in line.chars() {
                if c == '-' {
                    gapcount += 1;
                }
            }
            count += gapcount;
        }
    }

    println!("Count: {}", count);
}

fn count_dupe_refs(input: &String) {
    // Open file
    let file = std::fs::File::open(input).unwrap();
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut count = 0;
    let mut state = MafReadState::BlankLine;
    let mut seqcount = 0;
    let mut reference = String::new();
    let mut seqlengths = 0;

    while let Some(line) = lines.next() {
        let line = line.unwrap();

        // Match on the first character
        match line.chars().next() {
            Some('#') => {
                state = MafReadState::Comment;
                continue;
            }
            Some('a') => {
                state = MafReadState::AlignmentBlockLine;
                seqcount = 0;
                continue;
            }
            Some('s') => {
                state = MafReadState::SequenceLine;
                seqcount += 1;
            }
            None => {
                state = MafReadState::BlankLine;
                continue;
            }
            _ => {
                state = MafReadState::BlankLine;
                continue;
            }
        }

        if seqcount == 1 {
            // Get reference ID
            // Split: s Kakapo.S18                                14671 248
            let mut split = line.split_whitespace();
            reference = split.nth(1).unwrap().to_string();
        } else if seqcount > 1 {
            // Get sequence ID
            // Split: s Kakapo.S18                                14671 248
            let mut split = line.split_whitespace();
            let seqid = split.nth(1).unwrap().to_string();

            if seqid == reference {
                // println!("{:#?}", split);
                count += 1;
                seqlengths += split.nth(4).unwrap().len();
            }
        }
    }

    println!("Count: {}", count);
    println!("Length: {}", seqlengths);
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

// remove_dupe_ref_blocks
fn remove_dupe_ref_blocks(input: &String) {
    // Open file
    let file = std::fs::File::open(input).unwrap();
    let mut parser = maf_parser(file);
    let mut reference = String::new();
    let mut refcount = 0;
    let mut seqcount = 0;
    let mut removed_count = 0;

    'outer: while let Some(block) = parser.next() {
        seqcount = 0;
        refcount = 0;

        for line in block.iter() {
            match line {
                MafLine::SequenceLine(seqid, _, _, _, _, _) => {
                    seqcount += 1;

                    let speciesid = seqid.split('.').next().unwrap();

                    if seqcount == 1 {
                        reference = speciesid.to_string();
                    } else if seqcount > 1 && *speciesid == reference {
                        removed_count += 1;
                        continue 'outer;
                    }
                    // println!("{}, {}, {}", seqid, reference, seqcount);
                }
                _ => {}
            }
        }

        for line in block {
            println!("{}", line);
        }
        println!();
    }

    // Print to STDERR
    // eprintln!("Removed {} blocks", removed_count);
}

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

#[derive(Clone, Debug, PartialEq)]
struct GERPElem {
    region: String,
    start: u64,
    end: u64,
    score: f32,
    pval: f64,
    expected: f64,
    observed: f64,
    length: f64,
}

// Handy one liner function to create GERPElem
impl GERPElem {
    fn new(
        region: String,
        start: u64,
        end: u64,
        score: f32,
        pval: f64,
        expected: f64,
        observed: f64,
        length: f64,
    ) -> GERPElem {
        GERPElem {
            region,
            start,
            end,
            score,
            pval,
            expected,
            observed,
            length,
        }
    }

    // Is given position within this element?
    fn contains(&self, pos: u64) -> bool {
        if pos >= self.start && pos <= self.end {
            return true;
        } else {
            return false;
        }
    }
}

// process_gerp(maf, rates, elems, output);
fn process_gerp(
    maf: &String,
    rates: &String,
    elems: &String,
    output: &String,
) {
    let maf_fh = std::fs::File::open(maf).expect("Unable to open maf file");
    let rates_fh = std::fs::File::open(rates).expect("Unable to open rates file");
    let elems_fh = std::fs::File::open(elems).expect("Unable to open elems file");
    let output_fh = std::fs::File::create(output).expect("Unable to create file");

    let mut rates_fh = std::io::BufReader::new(rates_fh);
    let mut elems_fh = std::io::BufReader::new(elems_fh);
    let mut output_fh = std::io::BufWriter::new(output_fh);

    // Read first line of MAF file
    let mut maf_parser = maf_parser(maf_fh);
    let mut maf_block = maf_parser.next();

    // If none, then we're done
    if maf_block.is_none() {
        return;
    }

    let mut maf_block = maf_block.unwrap();

    // Get length of chromosome
    let mut chrom_length = 0;
    for line in maf_block.iter() {
        match line {
            // MafLine::SequenceLine(seqid, start, length, strand, src_size, text) => {
            MafLine::SequenceLine(_, _, _, _, length, _) => {
                chrom_length = *length;
                break;
            }
            _ => {}
        }
    }

    // Count lines in rates file
    let mut rates_count = 0;
    // Borrow lines iterator
    let mut lines = rates_fh.lines();
    while let Some(_) = lines.next() {
        rates_count += 1;
    }

    // GERP Sanity Check
    println!("Rates Length: {} Chrom Length: {}", rates_count, chrom_length);
    assert!(rates_count <= chrom_length, "Rates file length does not match chromosome length - This program expect GERP scores to be split by chromosome");

    // Rewind file
    let rates_fh = std::fs::File::open(rates).expect("Unable to open rates file");
    let mut rates_fh = std::io::BufReader::new(rates_fh);

    // Import elems file
    let mut elems: Vec<GERPElem> = elems_fh
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let split: Vec<&str> = line.split_whitespace().collect();

            let region = split[0].to_string();
            let start = split[1].parse::<u64>().unwrap();
            let end = split[2].parse::<u64>().unwrap();
            let score = split[3].parse::<f32>().unwrap();
            let pval = split[4].parse::<f64>().unwrap();
            let expected = split[5].parse::<f64>().unwrap();
            let observed = split[6].parse::<f64>().expect(format!("Unable to parse observed value: {}", split[6]).as_str());
            let length = split[7].parse::<f64>().expect(format!("Unable to parse length value: {}", split[7]).as_str());

            GERPElem::new(region, start, end, score, pval, expected, observed, length)
        })
        .collect();

    // Sort elems by start position
    elems.sort_by(|a, b| a.start.cmp(&b.start));

    let mut pos = 0;

    // Iterate through rates file
    let mut lines = rates_fh.lines();
    let mut in_elem = false;
    while let Some(line) = lines.next() {
        pos += 1;
        in_elem = false;

        let line = line.unwrap();
        let split: Vec<&str> = line.split_whitespace().collect();

        let neutral_rate = split[0].parse::<f64>().expect(format!("Unable to parse neutral rate: {}", split[0]).as_str());
        let score = split[1].parse::<f32>().unwrap();

        // If both are zero skip
        if neutral_rate == 0.0 && score == 0.0 {
            continue;
        }

        // Find element that contains this position
        let mut elem: Option<&GERPElem> = None;
        let mut remove = Vec::new();
        for i in 0..elems.len() {
            let e = &elems[i];
            if e.contains(pos) {
                elem = Some(e);
                in_elem = true;
                break;
            }

            // If pos is past end of element, remove it from vector
            if pos > e.end {
                remove.push(i);
            }

            // If pos is past end of chromosome, break
            if pos > chrom_length {
                break;
            }

            // If pos is before start of element, break
            if pos < e.start {
                break;
            }
        }

        // Write to output
        match elem {
            Some(e) => {
                writeln!(
                    output_fh,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    pos, neutral_rate, score, in_elem, e.start, e.end, e.score, e.pval, e.expected, e.observed, e.length
                )
                .unwrap();
            }
            None => {
                // Write to output, NA for element if not contained within a conserved element
                writeln!(
                    output_fh,
                    "{}\t{}\t{}\t{}\t\t\t\t\t\t\t",
                    pos, neutral_rate, score, in_elem
                ).unwrap();
            }
        }

        // Remove elements that are before the current position
        for i in remove.iter() {
            elems.remove(*i);            
        }
        
    }
}
