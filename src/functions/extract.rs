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
                        assert!(
                            species == reference.as_ref().unwrap(),
                            "First line of alignment block is not the reference genome"
                        );
                        alignment_block.seqid = seqid.clone();
                        alignment_block.start = *start;
                    }

                    alignment_block.add_line(line.clone());
                }
            }
        }
    }
}
