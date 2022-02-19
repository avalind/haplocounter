use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub enum EvidenceType {
    Maternal,
    Transmitted,
    Reference,
    Paternal,
}

pub struct EvidenceTable {
    data: HashMap<String, HashMap<(u32, i64), Vec<(EvidenceType, i32)>>>,
}

#[derive(Debug)]
pub struct BarcodeTally {
    pub barcode: String,
    pub transmitted_count: i32,
    pub maternal_count: i32,
    pub reference_count: i32,
    pub paternal_count: i32,
}

impl BarcodeTally {
    pub fn new(barcode: &str) -> BarcodeTally {
        BarcodeTally{
            barcode: barcode.to_string(),
            transmitted_count: 0,
            maternal_count: 0,
            reference_count: 0,
            paternal_count: 0,
        }
    }

    #[inline(always)]
    pub fn update(&mut self, offset: (i32, i32, i32, i32)) {
        self.maternal_count += offset.0;
        self.transmitted_count += offset.1;
        self.reference_count += offset.2;
        self.paternal_count += offset.3;
    }

    pub fn haplotype_class(self) -> EvidenceType {
        if self.maternal_count <= self.transmitted_count {
            EvidenceType::Transmitted
        } else {
            EvidenceType::Maternal
        }
    }
}

impl EvidenceTable {
    pub fn new() -> EvidenceTable {
        let t = EvidenceTable{ data: HashMap::new() };
        t
    }

    pub fn update_site_data(&mut self, barcode: &str, location: (u32, i64), count: i32, evt: EvidenceType) {
        if !self.data.contains_key(barcode) {
            self.data.insert(barcode.to_string(), HashMap::new());
        };

        // subtab is now a mutable reference into the hashmap, so we update through it.
        let subtab = match self.data.get_mut(barcode) {
            Some(v) => v,
            None => { panic!("barcode not found in table!\n"); },
        };

        if !subtab.contains_key(&location) {
            let mut v: Vec<(EvidenceType, i32)> = Vec::new();
            v.push((evt, count));
            subtab.insert(location, v); 
        } else {
            let vecref = match subtab.get_mut(&location) {
                Some(v) => v,
                None => { panic!("location not found in subtable\n"); },
            };
            vecref.push((evt, count));
        }       
    }

    pub fn prettyprint(self) {
        for (k, v) in self.data.iter() {
            println!("{} = {:?}", k, v);
        }
    }

    pub fn tally_barcode(data: &HashMap<(u32, i64), Vec<(EvidenceType, i32)>>, barcode: &str) -> BarcodeTally {
        let mut tally = BarcodeTally::new(barcode);
        for (k, v) in data.iter() {
            let adds = {
                let mut mats: i32 = 0;
                let mut trans: i32 = 0;
                let mut refs: i32 = 0;
                let mut pat: i32 = 0;
                v.iter().for_each(|v| {
                    match *v {
                        (EvidenceType::Maternal, count) => mats += count,
                        (EvidenceType::Transmitted, count) => trans += count,
                        (EvidenceType::Reference, count) => refs += count,
                        (EvidenceType::Paternal, count) => pat += count,
                    }
                });
                (mats, trans, refs, pat)
            };
            tally.update(adds);
        }
        tally  
    }

    pub fn classify(self) {
        let mut tallies: Vec<BarcodeTally> = Vec::new();
        for (k, v) in self.data.iter() {
            let barcode_summary = EvidenceTable::tally_barcode(&v, k);
            println!("{}\t{:?}", k, barcode_summary.haplotype_class());
        }
    }
}   

/*fn main() {
    let mut tab = EvidenceTable::new();
    tab.update_site_data(&"ACGTGA", (20, 12334534), 1, EvidenceType::Maternal);   
    tab.update_site_data(&"ACGTGA", (20, 12334534), 1, EvidenceType::Reference);
    tab.update_site_data(&"ACGTTT", (20, 12334112), 1, EvidenceType::Transmitted);
    tab.update_site_data(&"ACGTGA", (20, 12312321), 1, EvidenceType::Maternal);
    tab.prettyprint()
}*/
