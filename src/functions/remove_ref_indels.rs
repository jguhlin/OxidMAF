//! Extract SNPs from a MAF file

use std::collections::HashSet;

pub use crate::*;

// todo optional pass in reference genome
pub fn remove_ref_indels(maf: &String, output_prefix: &String) {

    let maf_fh = std::fs::File::open(maf).expect("Unable to open maf file");
    let maf_parser = maf_parser(maf_fh);

    let mut reference = None;

    let mut alignment_block;

    let mut windows = Vec::new(); // And write the windows to a file in samtools format

    for block in maf_parser {
        alignment_block = AlignmentBlock::default();
        for line in block.iter() {
            match line {
                // Do nothing for these two...
                MafLine::Comment(_) => (),
                MafLine::BlankLine => (),

                // This should trigger the start of a new alignment block
                // i.e., process the previous alignment block, then reset the state
                MafLine::AlignmentBlockLine(_) => (),

                MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) => {

                    if reference.is_none() {
                        reference = Some(species.clone());
                        alignment_block.reference = species.clone();
                    }

                    // If the alignment block is empty, this is the first line of the block
                    if alignment_block.lines.is_empty() {
                        assert!(species == reference.as_ref().unwrap(), "First line of alignment block is not the reference genome");
                        alignment_block.seqid = seqid.clone();
                        alignment_block.start = *start;
                    }

                    alignment_block.add_line(line.clone());
                }
            }
        }

        if alignment_block.is_empty() {
            continue;
        }

        alignment_block.remove_ref_indels();

        // For each block, write the alignment block to a file
        let ref_line = &alignment_block.lines[0];
        let window;

        if let MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) = ref_line {
            if *length < 2000 {
                println!("Skipping block: {} - Short block", seqid);
                continue;
            }
            window = format!("{}:{}-{}", seqid, start+1, start+length);
            assert!(*strand == Strand::Plus, "Strand is not positive");            
        } else {
            unreachable!();
        }

        alignment_block.calc_conservation();
        if alignment_block.stats.as_ref().unwrap().conservation < 0.95 {
            println!("Skipping block: {} - Low Conservation", window);
            continue;
        }

        let output_file = format!("{}-{}.fasta", output_prefix, &window);
        windows.push(window.clone());

        // Output file is prefix + window + .fasta
        let mut output_fh = std::fs::File::create(output_file).expect("Unable to create output file");

        // For each line, print species, then sequence
        let mut used_species = HashSet::new();
        for line in alignment_block.lines.iter() {
            if let MafLine::SequenceLine(species, _, _, _, _, _, text) = line {
                if used_species.contains(species) {
                    continue;
                }
                used_species.insert(species);
                writeln!(output_fh, ">{}", species).expect("Unable to write to output file");
                writeln!(output_fh, "{}", text).expect("Unable to write to output file");
            }
        }
    }

    // Write the windows to a file
    let window_file = format!("{}.windows", output_prefix);
    let mut window_fh = std::fs::File::create(window_file).expect("Unable to create window file");
    for window in windows.iter() {
        writeln!(window_fh, "{}", window).expect("Unable to write to window file");
    }
}

// For accumulating the alignment block before processing
#[derive(Default)]
pub struct AlignmentBlock {
    reference: String,
    lines: Vec<MafLine>,
    seqid: String,
    start: u64,
    stats: Option<AlignmentBlockStats>,
}

impl AlignmentBlock {
    pub fn add_line(&mut self, line: MafLine) {
        self.lines.push(line);
    }

    pub fn remove_ref_indels(&mut self) {

        let mut columns_to_remove: Vec<usize> = Vec::new();

        let reference = &self.lines[0];
        if let MafLine::SequenceLine(_, _, _, _, _, _, reference_text) = reference {
            // We only look for indels in the reference genome
            for (i, c) in reference_text.chars().enumerate() {
                if c == '-' {
                    columns_to_remove.push(i);
                }
            }           
        } else {
            unreachable!();
        }

        // Remove the columns from the alignment block
        for line in self.lines.iter_mut() {
            if let MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) = line {
                let new_text: String = text.chars().enumerate().filter(|(i, _)| !columns_to_remove.contains(i)).map(|(_, c)| c).collect();
                *text = new_text;
            }
        }

    }

    pub fn len(&self) -> usize {
        self.lines.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    pub fn calc_conservation(&mut self) {
        let mut conservation = 0.0;
        let mut total = 0;

        let mut ref_length = None;

        for line in self.lines.iter() {
            if let MafLine::SequenceLine(_, _, _, length, _, _, text) = line {
                if ref_length.is_none() {
                    ref_length = Some(*length);
                }

                total += ref_length.unwrap();
                
                for c in text.chars() {
                    if c != '-' {
                        conservation += 1.0;
                    }
                }
            }
        }

        self.stats = Some(AlignmentBlockStats {
            conservation: conservation / total as f32,
        });
    }


}

pub struct AlignmentBlockStats {
    pub conservation: f32,
}