use crate::parsers::*;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom, Write};

// This function now expects that a corresponding TAF index (.tai) file exists.
// It also assumes that the TAF file is uncompressed (or that seeking is supported).
pub fn annotate_ancestral_allele(taf: &String, vcf: &String, ancestors: &String, output: &String) {
    // Assert ancestors are a list of numbers, separated by commas.
    assert!(
        ancestors.chars().all(|x| x.is_numeric() || x == ','),
        "Ancestors must be a list of numbers separated by commas"
    );

    let mut out_aa = std::fs::File::create(format!("{}.tsv", output)).unwrap();
    let mut out_aa = std::io::BufWriter::new(out_aa);
    // Write header
    out_aa.write_all(b"#CHROM\tPOS\tAA\n").unwrap();

    let mut out_remove = std::fs::File::create(format!("{}_remove.tsv", output)).unwrap();
    let mut out_remove = std::io::BufWriter::new(out_remove);

    let ancestors: Vec<usize> = ancestors.split(',').map(|x| x.parse().unwrap()).collect();

    let mut vcf_reader = VcfParser::from_file(vcf); // Assume you have a VCF parser.

    // Load the TAF index from the .tai file (assumed to be at taf+".tai")
    let taf_index_path = format!("{}.tai", taf);
    let taf_index = TafIndex::load(&taf_index_path).expect("Failed to load TAF index");

    // Keep track of the current chromosome and alignment iterator.
    let mut current_chrom = String::new();
    // We use an Option so we can reinitialize when the chromosome changes.
    let mut align_iter: Option<TafAlignmentIterator> = None;
    // Cache the current alignment column.
    let mut column: Option<TafAlignmentColumn> = None;

    let mut taffy = TafParser::from_file(taf.as_str()).unwrap();
    let taf_index = TaiIndex::from_file(&format!("{}.tai", taf)).unwrap();

    // Statistics counters.
    let mut ancestral_matches = 0;
    let mut ancestral_not_found = 0;
    let mut ancestral_missing = 0;
    let mut ancestral_by_majority = 0;

    for record in vcf_reader.records() {
        // When the record's chromosome changes, seek to that contig.
        if current_chrom != record.chrom {
            current_chrom = record.chrom.clone();
            // Use the TAF index to get an iterator for this contig.
            let seek_info = taf_index.get_seek_info(&current_chrom, 1).unwrap();
            drop(align_iter);
            taffy.seek_to(seek_info.0, seek_info.1);
            align_iter = Some(TafAlignmentIterator::new(&mut taffy));
            column = align_iter
                .as_mut()
                .and_then(|it| it.next().transpose().ok().flatten());
        }

        // If we don't have a current column, skip the record.
        if column.is_none() {
            continue;
        }
        let mut col = column.clone().unwrap();

        // Advance alignment columns until the reference coordinate (row 0) matches the VCF position.
        while !col.ref_matches_pos(0, record.pos) {
            if let Some(next_result) = align_iter.as_mut().unwrap().next() {
                match next_result {
                    Ok(next_col) => {
                        // If the new column isn’t for the current chromosome,
                        // then we've reached the end of this contig's block.
                        if next_col.ref_coord(&record.chrom).is_none() {
                            col = next_col;
                            break;
                        }
                        col = next_col;
                    }
                    Err(e) => {
                        eprintln!("Error reading alignment column: {}", e);
                        break;
                    }
                }
            } else {
                // End of alignment columns for this contig.
                col = col;
                break;
            }
        }

        // If the alignment column still does not match the SNP position, skip this record.
        if !col.ref_matches_pos(0, record.pos) {
            continue;
        }

        // Try to extract an ancestral allele from one of the specified ancestors.
        let mut ancestral_allele = None;
        for anc in ancestors.iter() {
            let species = format!("Anc{}", anc);
            if let Some((_, allele)) = col.allele_for_species(&species) {
                let allele = allele.to_ascii_uppercase();
                if allele != '-' && allele != 'N' && allele != '*' {
                    if allele == record.ref_.chars().nth(0).unwrap()
                        || allele == record.alt.chars().nth(0).unwrap()
                    {
                        ancestral_allele = Some(allele);
                        break;
                    }
                }
            }
        }

        if ancestral_allele.is_none() {
            ancestral_missing += 1;
            continue;
        }

        let ancestral_allele = ancestral_allele.unwrap();
        if ancestral_allele == record.ref_.chars().nth(0).unwrap() {
            ancestral_matches += 1;
        } else {
            ancestral_not_found += 1;
        }

        // Write the annotation to output.
        out_aa
            .write_all(
                format!("{}\t{}\t{}\n", record.chrom, record.pos, ancestral_allele).as_bytes(),
            )
            .unwrap();

        // Cache the column for the next record.
        column = Some(col.clone());
    }

    println!("Done");
    println!("Ancestral matches: {}", ancestral_matches);
    println!("Ancestral not found: {}", ancestral_not_found);
    println!("Ancestral missing: {}", ancestral_missing);
    println!("Ancestral by majority: {}", ancestral_by_majority);
    let total = ancestral_matches + ancestral_not_found;
    let fraction = ancestral_matches as f64 / total as f64;
    println!("Fraction ancestral matches: {}", fraction);
}
