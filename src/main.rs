use clap::Parser;
use rust_htslib::{bam, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader, Read as BcfRead, record as BcfRecord};

mod evidence;

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
    name_mother: Option<String>,

    /// the sample in the vcf file that corresponds to
    /// the proband
    #[clap(long)]
    name_proband: Option<String>,

    /// whether we should run in paternal mode
    #[clap(long)]
    paternal: bool,

    /// if we run in paternal mode, we need the sample name of the 
    /// father in the vcf file
    #[clap(long)]
    name_father: Option<String>,
}

fn tally_evidence(evt: &mut evidence::EvidenceTable, var: &BcfRecord::Record, pileups: &bam::pileup::Pileup, vartype: evidence::EvidenceType)
{
    for aln in pileups.alignments() {
        // Unpack the cellular barcode this read/alignment came from.
        let alnrec = aln.record();
        let barcode = match alnrec.aux(b"CB") {
            Ok(v) => {
                if let Aux::String(s) = v {
                    s
                } else {
                    panic!("CB tag found but wrong type!");
                }
            },
            _ => continue,
        };

        let readpos = match aln.qpos() {
            Some(v) => v,
            None => continue, // None here means that the read is deleted/ref missing and if so we skip to the next alignment.
        };
        let nt = aln.record().seq()[readpos] as char;
        // find the first alt allele and check that its only a SNV.
        let alt = {
            let alt_first = var.alleles()[1];
            if alt_first.len() > 1 {
                continue;
            } else {
                alt_first[0] as char
            }
        };

        if alt == nt {
            evt.update_site_data(&barcode, (var.rid().unwrap(), var.pos()), 1, vartype);
        } else {
            evt.update_site_data(&barcode, (var.rid().unwrap(), var.pos()), 1, evidence::EvidenceType::Reference);
        }
    }
}


fn main() {
    let args = Args::parse();   
    let mut vcf = Reader::from_path(&args.vcf_path).expect("Unable to open vcf");
    let mut aln = bam::IndexedReader::from_path(&args.bam_path).expect("Unable to open bam");
    let header_view = vcf.header().clone();
    let mut tab = evidence::EvidenceTable::new();    

    if args.paternal {
        let father_sample_id = match(args.name_father) {
            Some(name) => name,
            None => panic!("--name_father needs to be supplied when running in paternal mode"),
        };
        let father_id = header_view.sample_id(&father_sample_id.as_bytes()).unwrap();

    } else {
        let name_mother = match(args.name_mother) {
            Some(name) => name,
            None => panic!("--name_mother needs to be supplied when running in non-paternal mode"),
        };

        let mother_id = header_view.sample_id(&name_mother.as_bytes()).unwrap();

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
                        tally_evidence(&mut tab, &rec, &pileup, evidence::EvidenceType::Transmitted);
                    } else {
                        tally_evidence(&mut tab, &rec, &pileup, evidence::EvidenceType::Maternal);
                    }
                
                }
            }
        }
        tab.classify();
    }
}
