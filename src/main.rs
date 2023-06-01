use bevy_app::prelude::*;
use bevy_app::App;
use bevy_ecs::prelude::*;
use bevy_tasks::TaskPool;
use clap::{Parser, Subcommand};

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
    }
}

// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
enum MafReadState {
    Comment,
    AlignmentBlockLine,
    SequenceLine,
    BlankLine,
}

enum Strand {
    Plus,
    Minus,
}

// https://genome.ucsc.edu/FAQ/FAQformat.html#format5
enum MafLine {
    Comment(String),
    AlignmentBlockLine(String),
    SequenceLine(String, u64, u64, Strand, u64, String),
    BlankLine,
}

fn parse_maf_line<'a>(line: &'a str) -> MafLine {
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
            
        }
        _ => {
            
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
