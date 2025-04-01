use flate2::read::MultiGzDecoder;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
##fileformat=VCFv4.2
##filedate=20250318
##source="beagle.29Oct24.c8e.jar"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record  (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + 2*P(AA)]">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  2311-06 2311-128        2311-59 2311-60 2311-66 2311-73 2311-78 A1      A10     A11     A12     A13     A14     A15     A16     A17     A19     A2      A3      A4      A5      A6      A7      A8   A9       C1      C100C1  C100C2  C11     C12     C13     C14     C15     C16     C17     C18     C2      C20     C21     C23     C25     C26     C27     C28     C3      C30     C31     C32     C33     C34     C35     C36     C38     C39     C4      C40     C41     C42  C43      C44     C45     C48     C49     C5      C51     C54     C56     C57     C58     C59     C6      C60     C61     C62     C63     C64     C65     C66     C67     C7      C72     C73     C74     C75     C76     C8      C9      CE9     E10     E11     E12     E13  E14      E15     E18     E19     E2      E20     E21     E22     E23     E24     E25     E27     E29     E3      E31     E32     E34     E35     E36     E37     E4      E42     E43     E44     E46     E5      E6      EN1     ENA1    HC24    HC29    HC30    HC32    HC33 HC34     HC37    HC39    HC40    HO10    HO11    HO12    HO13    HO15    HO16    HO17    HO19    HO2     HO3     HO4     HO9     N23-1   N23-10  N23-12  N23-13  N23-14  N23-15  N23-16  N23-17  N23-18  N23-19  N23-2   N23-20  N23-21  N23-22  N23-23  N23-24  N23-25  N23-26N23-27  N23-3   N23-31  N23-4   N23-5   N23-6   N23-7   N23-8   N23-9   P1      P13     P14     P15     P16     P17     P18     P19     P2      P20     P21     P22     P23     P24     P25     P26     P27     P28     P29     P3      P31     P32     P33     P34     P35  P36      P40     P41     P42     P43     P5      P6      P7      P9      PP10    PP11    PP12    PP14    PP15    PP16    PP18    PP2     PP20    PP21    PP22    PP23    PP24    PP25    PP26    PP27    PP28    PP29    PP3     PP30    PP31    PP32    PP33    PP34    PP35 PP36     PP37    PP38    PP4     PP40    PP41    PP42    PP43    PP44    PP46    PP7     PP9     YEP12   YEP13   YEP14   YEP15   YEP16   YEP17   YEP18   YEP19   YEP20   YEP21   YEP22
ptg000001l      10248   ptg000001l_10248_A_C    A       C       .       PASS    .       GT      1|1     1|0     0|0     0|1     0|1     1|1     1|1     1|1     1|1     1|0     0|1     1|0     0|0     1|1     1|1     1|0     0|1     1|1     1|1     1|1     1|1     0|1  1|1      1|1     0|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|0     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1  1|1      1|1     1|1     1|1     0|1     1|1     0|1     1|1     1|0     0|1     1|1     1|0     0|1     1|1     1|1     0|0     1|1     0|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1  1|0      0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|0     0|1     0|0     1|1     0|0     1|1     1|0     0|0     0|1     1|0     0|1     0|1     1|1     0|1     1|1     1|0     0|1     0|1     1|1     0|0     0|1     1|1     0|0     0|1     1|1     1|1  1|1      1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|0     0|1     1|0     0|1     1|1     1|0     0|1     1|0     0|1     1|0     0|0     0|0     0|1  1|1      1|1     0|1     1|1     0|0     1|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|0     0|0     0|1     1|1     1|1     1|1     1|1     1|0     0|1     0|1     1|1     1|0     0|1     1|0     0|1     0|0     1|1     1|1     1|1     1|0     0|1  0|1      1|1     1|0     0|1     1|0     0|0     0|1     1|1     1|0     0|1     0|1     1|1     0|1     1|1     0|0     0|1     1|1     1|1     1|1     1|1     1|0     0|0     1|1     1|1     1|1     1|1     1|1     0|1     1|1     1|1     1|1     1|1     1|1     1|0  0|1      0|1     1|1     1|1     1|0     0|1     0|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     0|1     1|1     0|1     1|1     1|0     0|1     0|0     0|1     1|1     1|1     1|0
ptg000001l      13673   ptg000001l_13673_T_C    T       C       .       PASS    .       GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0  0|0      0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0  0|0      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0  0|0      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|1     1|0     0|0     0|0  0|0      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0  0|0      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0  0|0      0|0     0|1     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|1  0|0      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     1|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0
*/

pub struct VcfParser {
    version: String,
    info: HashMap<String, String>,
    format: HashMap<String, String>,
    samples: Vec<String>,
    reader: Box<dyn BufRead>,
}

impl VcfParser {
    pub fn from_file(file: &str) -> Self {
        // Is the file compressed? Does it end with .gz?
        let mut reader: Box<dyn BufRead> = if file.ends_with(".gz") {
            let file = File::open(file).expect("Error opening VCF file");
            let decoder = MultiGzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            let file = File::open(file).expect("Error opening VCF file");
            Box::new(BufReader::new(file))
        };

        let (version, info, format, samples) = Self::parse_header(&mut reader);

        Self {
            version,
            info,
            format,
            samples,
            reader,
        }
    }

    fn parse_header(
        reader: &mut Box<dyn BufRead>,
    ) -> (
        String,
        HashMap<String, String>,
        HashMap<String, String>,
        Vec<String>,
    ) {
        let mut version = String::new();
        let mut info = HashMap::new();
        let mut format = HashMap::new();
        let samples;

        let mut line = String::new();
        loop {
            line.clear();
            reader
                .read_line(&mut line)
                .expect("Error reading VCF header");

            if line.starts_with("##fileformat=") {
                version = line
                    .split('=')
                    .nth(1)
                    .expect("Error parsing VCF version")
                    .trim()
                    .to_string();
            } else if line.starts_with("##INFO") {
                let id = line
                    .split(',')
                    .nth(0)
                    .expect("Error parsing INFO ID")
                    .split('=')
                    .nth(1)
                    .expect("Error parsing INFO ID")
                    .to_string();
                let description = line
                    .split(',')
                    .nth(3)
                    .expect("Error parsing INFO description")
                    .split('=')
                    .nth(1)
                    .expect("Error parsing INFO description")
                    .to_string();
                info.insert(id, description);
            } else if line.starts_with("##FORMAT") {
                let id = line
                    .split(',')
                    .nth(0)
                    .expect("Error parsing FORMAT ID")
                    .split('=')
                    .nth(1)
                    .expect("Error parsing FORMAT ID")
                    .to_string();
                let description = line
                    .split(',')
                    .nth(3)
                    .expect("Error parsing FORMAT description")
                    .split('=')
                    .nth(1)
                    .expect("Error parsing FORMAT description")
                    .to_string();
                format.insert(id, description);
            } else if line.starts_with("#CHROM") {
                samples = line
                    .split('\t')
                    .skip(9)
                    .map(|x| x.trim().to_string())
                    .collect();
                break;
            }
        }

        (version, info, format, samples)
    }

    // Return a records iterator (need to borrow?)
    pub fn records(&mut self) -> VcfRecords {
        VcfRecords { parser: self }
    }
}

pub struct VcfRecords<'a> {
    parser: &'a mut VcfParser,
}

impl<'a> Iterator for VcfRecords<'a> {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();

        loop {
            line.clear();
            self.parser
                .reader
                .read_line(&mut line)
                .expect("Error reading VCF record");

            if line.is_empty() {
                return None;
            }

            if line.starts_with("#") {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();

            let chrom = fields[0].to_string();
            let pos = fields[1].parse().expect("Error parsing VCF position");
            let id = fields[2].to_string();
            let ref_ = fields[3].to_string();
            let alt = fields[4].to_string();
            let qual = fields[5].to_string();
            let filter = fields[6].to_string();

            let mut info = HashMap::new();
            for field in fields[7].split(';') {
                let parts: Vec<&str> = field.split('=').collect();
                if parts.len() == 2 {
                    info.insert(parts[0].to_string(), parts[1].to_string());
                } else {
                    info.insert(parts[0].to_string(), "".to_string());
                }
            }

            let mut format = HashMap::new();
            for field in fields[8].split(':') {
                format.insert(field.to_string(), "".to_string());
            }

            let mut samples = Vec::new();
            for sample in fields.iter().skip(9) {
                let mut sample_fields = HashMap::new();
                for (key, value) in format.iter() {
                    sample_fields.insert(key.to_string(), value.to_string());
                }
                for (key, value) in sample.split(':').zip(sample_fields.iter_mut()) {
                    *value.1 = key.to_string();
                }
                samples.push(sample_fields);
            }

            return Some(Record {
                chrom,
                pos,
                id,
                ref_,
                alt,
                qual,
                filter,
                info,
                format,
                samples,
            });
        }
    }
}

pub struct Record {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: HashMap<String, String>,
    pub format: HashMap<String, String>,
    pub samples: Vec<HashMap<String, String>>,
}
