use bevy_app::prelude::*;
use bevy_app::App;
use bevy_ecs::prelude::*;
use bevy_tasks::TaskPool;
use clap::{Parser, Subcommand};

use std::fmt::Display;
use std::io::BufRead;

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

pub fn analyze() {
    App::new().add_plugin(bevy_app::ScheduleRunnerPlugin).run();
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
