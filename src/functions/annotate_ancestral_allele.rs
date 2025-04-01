use std::io::Write;

use crate::parsers::*;
pub fn annotate_ancestral_allele(taf: &String, vcf: &String, ancestors: &String, output: &String) {
    // Assert ancestors are a list of numbers, separated by commas.
    assert!(
        ancestors.chars().all(|x| x.is_numeric() || x == ','),
        "Ancestors must be a list of numbers separated by commas"
    );

    let mut out_aa = std::fs::File::create(format!("{}.tsv", output)).unwrap();
    let mut out_aa = std::io::BufWriter::new(out_aa);

    // Add header
    // #CHROM  POS     AA
    out_aa.write_all(b"#CHROM\tPOS\tAA\n").unwrap();

    let mut out_remove = std::fs::File::create(format!("{}_remove.tsv", output)).unwrap();
    let mut out_remove = std::io::BufWriter::new(out_remove);

    let ancestors: Vec<usize> = ancestors.split(',').map(|x| x.parse().unwrap()).collect();

    let mut vcf_reader = VcfParser::from_file(vcf); // Assume you have a VCF parser.
    let taf_file = std::fs::File::open(taf).unwrap();
    let taf_file = flate2::read::MultiGzDecoder::new(taf_file);
    let buf_reader = std::io::BufReader::new(taf_file);

    let mut parser = TafParser::new(buf_reader).unwrap();
    let mut align_iter = TafAlignmentIterator::new(&mut parser);

    let mut ancestral_matches = 0;
    let mut ancestral_not_found = 0;
    let mut ancestral_missing = 0;
    let mut ancestral_by_majority = 0;

    // Get the first alignment column.
    let mut column = align_iter.next().unwrap().unwrap();

    for record in vcf_reader.records() {
        // First, ensure the current column is for the correct chromosome.
        if column.ref_coord(&record.chrom).is_none() {
            println!(
                "Record: {}:{} ({}), Chrom: {}",
                record.chrom, record.pos, record.id, record.chrom
            );
            println!(
                "Alignment column reference: {:?}",
                column.ref_coord(&record.chrom)
            );
            println!("Ancestral matches: {}", ancestral_matches);
            println!("Ancestral not found: {}", ancestral_not_found);
            println!("Ancestral missing: {}", ancestral_missing);
            println!("Ancestral by majority: {}", ancestral_by_majority);

            // Print fraction
            let total = ancestral_matches + ancestral_not_found;
            let fraction = ancestral_matches as f64 / total as f64;
            println!("Fraction ancestral matches: {}", fraction);

            // TODO Temp while only working on one chr
            break;
        }

        // Now, loop until the reference coordinate in row 0 matches the SNP position.
        // (Assume row 0 is the reference.)
        while !column.ref_matches_pos(0, record.pos) {
            // Advance to the next column.
            column = align_iter.next().unwrap().unwrap();
            // Also check that this new column is on the correct chromosome.
            if column.ref_coord(&record.chrom).is_none() {
                println!(
                    "Record: {}:{} ({}), Chrom: {}",
                    record.chrom, record.pos, record.id, record.chrom
                );
                println!(
                    "Alignment column reference: {:?}",
                    column.ref_coord(&record.chrom)
                );
                println!("Ancestral matches: {}", ancestral_matches);
                println!("Ancestral not found: {}", ancestral_not_found);
                println!("Ancestral missing: {}", ancestral_missing);
                println!("Ancestral by majority: {}", ancestral_by_majority);

                // Print fraction
                let total = ancestral_matches + ancestral_not_found;
                let fraction = ancestral_matches as f64 / total as f64;
                println!("Fraction ancestral matches: {}", fraction);

                // TODO Temp while only working on one chr
                break;
            }
        }

        // At this point, column.ref_matches_pos(0, record.pos) is true.
        // Try to extract the SNP for an ancestor.
        let mut ancestral_allele = None;
        for anc in ancestors.iter() {
            let species = format!("Anc{}", anc);
            if let Some((_, allele)) = column.allele_for_species(&species) {
                let allele = allele.to_ascii_uppercase();
                if allele != '-' && allele != 'N' && allele != '*' {
                    if allele == record.ref_.chars().nth(0).unwrap()
                        || allele == record.alt.chars().nth(0).unwrap()
                    {
                        ancestral_allele = Some(allele);
                        break;
                    }
                }
                // println!("Ancestral allele for {} is not a valid base: {}", species, allele);
            }
            // println!("Ancestral allele not found for {}", species);
        }

        // Majority rule: not good for tsinfer though...
        /*
        if ancestral_allele.is_none() {
            // Try to get it via majority rule
            // (But exclude the reference species.)
            let mut allele_counts = std::collections::HashMap::new();
            for (species, allele) in column.column.alleles.iter().enumerate() {
                if species == 0 {
                    continue;
                }
                let allele = allele.to_ascii_uppercase();
                if allele != '-' && allele != 'N' && allele != '*' {
                    *allele_counts.entry(allele).or_insert(0) += 1;
                }
            }

            let mut max_count = 0;
            let mut max_allele = None;
            for (allele, count) in allele_counts.iter() {
                if *count > max_count {
                    max_count = *count;
                    max_allele = Some(allele);
                }
            }

            if let Some(allele) = max_allele {
                if *allele == record.ref_.chars().nth(0).unwrap()
                    || *allele == record.alt.chars().nth(0).unwrap()
                {
                    ancestral_allele = Some(allele.clone());
                    ancestral_by_majority += 1;
                }
            }
        }  */

        if ancestral_allele.is_none() {
            // println!("No ancestral allele found for SNP at {}:{} ({}). Skipping...", record.chrom, record.pos, record.id);
            // panic!("Ancestral allele not found");
            // continue;
            ancestral_missing += 1;
            continue;
        }

        let ancestral_allele = ancestral_allele.unwrap();
        if ancestral_allele == record.ref_.chars().nth(0).unwrap() {
            ancestral_matches += 1;
        } else {
            ancestral_not_found += 1;
        }

        // Write to output file
        // You need a tab-delimited file with CHROM, POS, and the AA field as the 3rd column.
        out_aa
            .write_all(
                format!(
                    "{}\t{}\t{}\n",
                    record.chrom, record.pos, ancestral_allele
                )
                .as_bytes(),
            )
            .unwrap();

    }
    println!("Done");

    println!("Ancestral matches: {}", ancestral_matches);
    println!("Ancestral not found: {}", ancestral_not_found);
    println!("Ancestral missing: {}", ancestral_missing);
    println!("Ancestral by majority: {}", ancestral_by_majority);

    // Print fraction
    let total = ancestral_matches + ancestral_not_found;
    let fraction = ancestral_matches as f64 / total as f64;
    println!("Fraction ancestral matches: {}", fraction);
}
