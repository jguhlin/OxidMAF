//! Extract SNPs from a MAF file

pub use crate::*;

// todo optional pass in reference genome
pub fn extract_snps(maf: &String, output_prefix: &String, coordinates: bool) {

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
    lines: Vec<MafLine>,
    seqid: String,
    start: u64,
}

impl AlignmentBlock {
    pub fn add_line(&mut self, line: MafLine) {
        self.lines.push(line);
    }

    pub fn extract_snps(&self) {
        // Go through each column in the alignment block
        // Find where they are not all a match, and extract the SNPs

        // let mut snps = Vec::new(); // Vec<Vec<char>> same order as the block





    }


}