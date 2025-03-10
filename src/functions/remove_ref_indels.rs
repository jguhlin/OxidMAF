//! Extract SNPs from a MAF file

pub use crate::*;

// todo optional pass in reference genome
pub fn remove_ref_indels(maf: &String, output_prefix: &String) {

    let maf_fh = std::fs::File::open(maf).expect("Unable to open maf file");
    let maf_parser = maf_parser(maf_fh);

    let mut reference = None;

    let mut alignment_block = AlignmentBlock::default();

    for block in maf_parser {
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



    }
}

// For accumulating the alignment block before processing
#[derive(Default)]
pub struct AlignmentBlock {
    reference: String,
    lines: Vec<MafLine>,
    seqid: String,
    start: u64,
}

impl AlignmentBlock {
    pub fn add_line(&mut self, line: MafLine) {
        self.lines.push(line);
    }

    pub fn remove_ref_indels(&self) {

        let mut columns_to_remove: Vec<usize> = Vec::new();

        let reference = &self.lines[0];
        if let MafLine::SequenceLine(_, _, _, _, _, _, reference_text) = reference {
            for (i, line) in self.lines.iter().enumerate() {
                if let MafLine::SequenceLine(_, _, _, _, _, _, text) = line {
                    for (j, (ref_base, base)) in reference_text.chars().zip(text.chars()).enumerate() {
                        if ref_base == '-' || base == '-' {
                            columns_to_remove.push(j);
                        }
                    }
                }
            }
        } else {
            unreachable!();
        }

        // Remove the columns from the alignment block
        for line in self.lines.iter() {
            if let MafLine::SequenceLine(species, seqid, start, length, strand, src_size, text) = line {
                let new_text: String = text.chars().enumerate().filter(|(i, _)| !columns_to_remove.contains(i)).map(|(_, c)| c).collect();
                println!("{}", new_text);
            }
        }

    }


}