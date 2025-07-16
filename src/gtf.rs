use crate::{isoform::FusionCandidate, utils::trim_chr_prefix_to_upper};
use anyhow::Result;
use core::panic;
use noodles_gtf::{self as GTF, io::Reader as gtfReader, record::Strand};
use rust_lapper::{Interval, Lapper};
use std::{collections::HashMap, io::BufRead};
#[derive(Debug, Clone)]
pub struct Transcript {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub splice_junc: Vec<(u64, u64)>,
    pub exons: Vec<(u64, u64)>,
    pub gene_id: String,
    pub trans_id: String,
    pub records: Vec<GTF::Record>,
    pub strand: Strand,
}

impl Transcript {
    pub fn process(&mut self) -> &mut Self {
        self.exons.sort_by_key(|k| k.0);

        if self.exons.len() >= 2 {
            self.splice_junc = self.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        } else {
            self.splice_junc = vec![(self.exons[0].0, self.exons[0].1)];
        }
        self
    }

    pub fn from_gtf_record(record: &GTF::Record) -> Self {
        let mut trans = Transcript {
            chrom: String::new(),
            splice_junc: Vec::new(),
            start: 0,
            end: 0,
            exons: Vec::new(),
            gene_id: String::new(),
            trans_id: String::new(),
            records: Vec::new(),
            strand: record.strand().expect("GTF must have strand information"),
        };

        trans.chrom = trim_chr_prefix_to_upper(record.reference_sequence_name());
        trans.gene_id = record
            .attributes()
            .iter()
            .find(|x| x.key() == "gene_id")
            .expect("GTF must have gene_id")
            .value()
            .to_string();
        trans.trans_id = record
            .attributes()
            .iter()
            .find(|x| x.key() == "transcript_id")
            .expect("GTF must have transcript_id")
            .value()
            .to_string();
        trans.start = record.start().get() as u64 - 1;
        trans.end = record.end().get() as u64;
        trans.records.push(record.clone());
        trans
    }

    pub fn get_quieries(&self) -> Vec<(String, u64)> {
        let mut queries = Vec::new();
        for junc in &self.splice_junc {
            queries.push((self.chrom.clone(), junc.0));
            queries.push((self.chrom.clone(), junc.1));
        }
        queries
    }

    pub fn get_attributes(&self) -> String {
        let s = self
            .records
            .iter()
            .map(|x| {
                x.attributes()
                    .iter()
                    .map(|y| format!("{}:{}", y.key(), y.value()))
                    .collect::<Vec<String>>()
                    .join(";")
            })
            .collect::<Vec<String>>()
            .join(";");
        s
    }

    pub fn get_exon_count(&self) -> usize {
        self.exons.len()
    }

    pub fn get_transcript_length(&self) -> u64 {
        self.end - self.start
    }
}

pub struct TranscriptChunker<R: BufRead> {
    pub gtfreader: gtfReader<R>,
    pub cur_chrom: String,
    pub cur_pos: u64,
    pub is_end: bool,
    pub trans_count: u64,
    pub hold_transcript: Option<Transcript>,
}

impl<R: BufRead> TranscriptChunker<R> {
    pub fn new(gtf_reader: gtfReader<R>) -> TranscriptChunker<R> {
        TranscriptChunker {
            gtfreader: gtf_reader,
            cur_chrom: String::new(),
            cur_pos: 0,
            is_end: false,
            trans_count: 0,
            hold_transcript: None,
        }
    }

    pub fn get_next_transcript(&mut self) -> Option<Transcript> {
        let mut records = self.gtfreader.records();

        loop {
            let record = records.next();
            // println!("Processing record: {:?}", &record);
            match record {
                Some(ref rec) => match rec {
                    Ok(rec) => match rec.ty() {
                        "transcript" => match self.hold_transcript {
                            Some(_) => {
                                let mut trans = self.hold_transcript.clone();

                                let new_trans = Transcript::from_gtf_record(&rec);
                                self.hold_transcript = Some(new_trans);
                                trans.as_mut().unwrap().process();
                                return trans;
                            }
                            None => {
                                let transcript = Transcript::from_gtf_record(&rec);
                                self.hold_transcript = Some(transcript);
                            }
                        },
                        "exon" => match self.hold_transcript {
                            Some(_) => {
                                self.hold_transcript
                                    .as_mut()
                                    .unwrap()
                                    .exons
                                    .push((rec.start().get() as u64 - 1, rec.end().get() as u64));
                            }
                            None => {}
                        },
                        _ => {}
                    },
                    Err(e) => {
                        panic!("error in reading gtf/gff: {:?}, record: {:?}\nPlease check your gtf/gff format", e, &record);
                    }
                },
                None => {
                    if self.is_end == false {
                        self.is_end = true;
                        let mut trans = self.hold_transcript.clone();
                        trans.as_mut().unwrap().process();
                        return trans;
                    } else {
                        return None;
                    }
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for TranscriptChunker<R> {
    type Item = Transcript;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_transcript()
    }
}

#[derive(Debug, Clone)]
pub struct GeneInterval {
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub splice_sites: Vec<u64>,
}

impl PartialEq for GeneInterval {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom
            && self.start == other.start
            && self.end == other.end
            && self.gene_id == other.gene_id
    }
}

impl Eq for GeneInterval {}

impl GeneInterval {
    pub fn match2(&self, splice_sites: &Vec<u64>, flank: u64) -> usize {
        let mut match_count = 0;
        for site in splice_sites.iter() {
            for gsite in self.splice_sites.iter() {
                if (*site >= gsite.saturating_sub(flank)) && (*site <= gsite + flank) {
                    match_count += 1;
                    break;
                }
            }
        }
        match_count
    }
}

type IV = Interval<u64, GeneInterval>;

pub struct GeneIntervalTree {
    pub tree: HashMap<String, Lapper<u64, GeneInterval>>,
    pub count: usize,
    pub chroms: Vec<String>,
    // pub list: Vec<GeneInterval>,
}

impl GeneIntervalTree {
    pub fn new<R: BufRead>(gtf_reader: &mut gtfReader<R>) -> Result<Self> {
        let mut chrom_tree = HashMap::new();

        let records = gtf_reader.records();
        let mut splice_sites = Vec::new();
        // let mut splice_set = std::collections::HashSet::new();
        // let mut splice_sites_exon = Vec::new();
        let mut chrom = String::new();
        let mut gene_id = String::new();
        let mut gene_name = String::new();
        let mut gene_start = 0;
        let mut gene_end = 0;
        let mut interval_count = 0;
        let mut chrom_set = std::collections::HashSet::new();
        // let mut interval_list = Vec::new();
        for record in records {
            let record = record?;
            chrom = trim_chr_prefix_to_upper(record.reference_sequence_name());
            chrom_set.insert(chrom.clone());
            chrom_tree
                .entry(chrom.clone())
                .or_insert_with(|| Vec::<GeneInterval>::new());
            match record.ty() {
                "gene" => {
                    if splice_sites.is_empty() {
                        gene_id = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_id")
                            .expect("GTF must have gene_id")
                            .value()
                            .to_string();
                        gene_name = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_name")
                            .expect("GTF must have gene_name")
                            .value()
                            .to_string();
                        gene_start = record.start().get() as u64 - 1;
                        gene_end = record.end().get() as u64;
                    } else {
                        // make new GeneInterval and reset splice_sites

                        splice_sites.sort();
                        splice_sites.dedup();

                        let gene_interval = GeneInterval {
                            gene_id: gene_id.clone(),
                            gene_name: gene_name.clone(),
                            chrom: chrom.clone(),
                            start: gene_start,
                            end: gene_end,
                            splice_sites: splice_sites.clone(),
                        };
                        let gene_intervals = chrom_tree.get_mut(&chrom).unwrap();
                        interval_count += 1;
                        gene_intervals.push(gene_interval);

                        // reload & reset
                        splice_sites.clear();

                        gene_id = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_id")
                            .expect("GTF must have gene_id")
                            .value()
                            .to_string();
                        gene_name = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_name")
                            .expect("GTF must have gene_name")
                            .value()
                            .to_string();
                        gene_start = record.start().get() as u64 - 1;
                        gene_end = record.end().get() as u64;
                    }
                }

                "exon" => {
                    splice_sites.push(record.start().get() as u64 - 1);
                    splice_sites.push(record.end().get() as u64);
                }

                _ => {}
            }
        }

        // make last GeneInterval
        if !splice_sites.is_empty() {
            splice_sites.dedup();
            let gene_interval = GeneInterval {
                gene_id: gene_id.clone(),
                gene_name: gene_name.clone(),
                chrom: chrom.clone(),
                start: gene_start,
                end: gene_end,
                splice_sites: splice_sites.clone(),
            };
            let gene_intervals = chrom_tree.get_mut(&chrom).unwrap();
            interval_count += 1;
            gene_intervals.push(gene_interval);
        }
        // convert to Lapper
        let mut intervals = HashMap::new();
        for (chrom, gene_intervals) in chrom_tree {
            let lapper = Lapper::new(
                gene_intervals
                    .into_iter()
                    .map(|x| Interval {
                        start: x.start,
                        stop: x.end,
                        val: x,
                    })
                    .collect::<Vec<IV>>(),
            );
            intervals.insert(chrom.clone(), lapper);
        }
        Ok(GeneIntervalTree {
            tree: intervals,
            count: interval_count,
            chroms: chrom_set.into_iter().collect(),
        })
    }

    pub fn find(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Option<Vec<Interval<u64, GeneInterval>>> {
        let chrom = trim_chr_prefix_to_upper(chrom);
        match self.tree.get(&chrom) {
            Some(lapper) => {
                let result = lapper
                    .find(start, end)
                    .map(|x| x.clone())
                    .collect::<Vec<Interval<u64, GeneInterval>>>();
                Some(result)
            }
            None => None,
        }
    }

    /// return a vec of gene part that match the splice sites in the sample read segments level
    pub fn match_splice_sites(
        &self,
        sample_read_splice: &FusionCandidate,
        flank: u64,
        min_match: usize,
    ) -> u64 {
        let mut match_count = 0;
        for seg_vec in sample_read_splice.supp_segments_by_read.iter() {
            let mut read_match_count = 0;
            for seg in seg_vec.iter() {
                if let Some(iv) = self.find(&seg.chrom, seg.start, seg.end) {
                    if iv.len() > 0 {
                        // the iv is the hit that overlap the supp segment and the gene interval from gtf.
                        // next iterately check if the splice site match
                        for gene_iv in iv.iter() {
                            // check if the gene interval match the splice site with flank and min_match with the read segment.
                            read_match_count +=
                                gene_iv.val.match2(&vec![seg.start, seg.end], flank);
                        }
                    }
                }
            }

            if read_match_count >= min_match {
                match_count += 1;
            }
        }

        match_count
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_gene_interval_tree() {
        let gtf_data_path =
            "/ssd1/stix-iso-devspace/isopedia-dev/test/gencode.v47.basic.annotation.gtf";
        let file = std::fs::File::open(gtf_data_path).expect("Failed to open GTF file");
        let reader = std::io::BufReader::new(file);
        let mut gtf_reader = gtfReader::new(reader);

        let gene_tree =
            GeneIntervalTree::new(&mut gtf_reader).expect("Failed to create GeneIntervalTree");
        println!("{}", gene_tree.count);

        // test find gene intervals
        let chrom = "1".to_string();
        let intervals = gene_tree.tree.get(&chrom);
        assert!(
            intervals.is_some(),
            "Chromosome {} not found in intervals",
            chrom
        );
        let intervals = intervals.unwrap();
        let target = intervals
            .find(13221, 14409)
            .collect::<Vec<&Interval<u64, GeneInterval>>>();
        println!("Found intervals: {:?}", target);

        intervals.iter().for_each(|iv| {
            println!(
                "{}:{}-{} {} {}",
                iv.val.chrom, iv.val.start, iv.val.end, iv.val.gene_id, iv.val.gene_name
            );
            let x = intervals
                .find(iv.val.start, iv.val.end)
                .collect::<Vec<&Interval<u64, GeneInterval>>>();
            println!("  found {} intervals： {:?}", x.len(), x);
        });
    }
}
