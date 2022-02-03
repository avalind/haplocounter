use clap::Parser;
use rust_htslib::{bam, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader, Read as BcfRead, record as BcfRecord};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// input vcf file with sites to pileup
    #[clap(long)]
    vcf_path: String,

    /// input bam file to perform pileups on
    #[clap(long)]
    bam_path: String,
    
    /// the sample in the vcf file that corresponds
    /// to the mother
    #[clap(long)]
    name_mother: String,

    /// the sample in the vcf file that corresponds to
    /// the proband
    #[clap(long)]
    name_proband: String,
}

enum EvidenceType {
    Maternal,
    Transmitted,
}

// a struct that represents evidence at
// a single position for either
struct PositionalEvidence {
    cid: i32,
    pos: i64,
    evtype: EvidenceType,
}


fn main() {
    let args = Args::parse();   
    let mut vcf = Reader::from_path(&args.vcf_path).expect("Unable to open vcf");
    let mut aln = bam::IndexedReader::from_path(&args.bam_path).expect("Unable to open bam");
    let header_view = vcf.header().clone();     

    let mother_id = header_view.sample_id(&args.name_mother.as_bytes()).unwrap();
    println!("{}", mother_id);

    for (idx, rr) in vcf.records().enumerate() {
        let mut rec = rr.expect("Unable to read record!");
        let ctg = match rec.rid() {
            Some(r) => r,
            None => panic!("malformed vcf file, panicking!"),
        };

        aln.fetch((header_view.rid2name(ctg).unwrap(), rec.pos(), rec.pos()+1));
        for p in aln.pileup() {
            let pileup = p.unwrap();
            if i64::from(pileup.pos()) == rec.pos() {
                // we are now in a position with a variant. next step is to check whether
                // this variant comes from the transmitted X or the one the mother did not
                // transmit. we do this by looking at the genotypes for the maternal sample
                // if they all are 0, we know that the variant is on the transmitted chromosome
                // (because of how the input vcf file is generated)
                let gt = rec.genotypes().expect("Unable to read genotypes!");               
                if gt.get(mother_id).iter().all(|x| *x == BcfRecord::GenotypeAllele::Unphased(0)) {
                    // we know now we are on a variant present _ONLY_ on the transmitted X,
                    // and we thus proceed to count the number nucleotides at this position.
                    for alignment in pileup.alignments() {
                        let alnrec = alignment.record();
                        let barcode = match alnrec.aux(b"CB") {
                            Ok(val) => { 
                                if let Aux::String(s) = val {
                                    s
                                } else {
                                    // this really shouldn't happen.
                                    panic!("Found CB tag but it isn't a string, bailing!");
                                }
                            },
                            _ => "",
                        };
                        
                        let readpos = match alignment.qpos() {
                            Some(v) => v,
                            None => continue,
                        };
    
                        if readpos >= alnrec.seq_len() {
                            panic!("readpos >= read length, this should REALLY NOT happen!");
                        };

                        let seqdata = alnrec.seq();
                        let nt = seqdata[readpos] as char;
                        let alt = rec.alleles()[1];
                        if alt.len() > 1 {
                            continue;
                        } else {
                            let altnuc = alt[0] as char;
                            println!("[{:?}] - [{:?}]", nt, altnuc);
                        }
                    }
                } else {
                    println!("Variant counted as coming from the kept maternal homologue!");
                }
                
                println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
            }
        }
    }
}
