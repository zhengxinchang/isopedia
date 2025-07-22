use std::io::BufRead;

use anyhow::Result;
use noodles_gtf::Reader as gtfReader;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;

use crate::utils::trim_chr_prefix_to_upper;

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
    pub fn match_splice_sites(&self, splice_sites: &Vec<u64>, flank: u64) -> Vec<u64> {
        // let mut match_count = 0;
        let mut matched_sites = Vec::new();
        for site in splice_sites.iter() {
            for gsite in self.splice_sites.iter() {
                if (*site >= gsite.saturating_sub(flank)) && (*site <= gsite + flank) {
                    // match_count += 1;
                    matched_sites.push(*site);
                    break;
                }
            }
        }
        matched_sites
    }
}

type IV = Interval<u64, GeneInterval>;

/// # GeneIntervalTree
/// A structure to hold gene intervals indexed by chromosome using Lapper.
/// It used in indexing the GTF file
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

    pub fn match2(
        &self,
        chrom: &str,
        splice_sites: &Vec<u64>,
        flank: u64,
    ) -> Option<(String, Vec<u64>)> {
        let chrom = trim_chr_prefix_to_upper(chrom);
        match self.tree.get(&chrom) {
            Some(lapper) => {
                let result = lapper
                    .find(splice_sites[0], splice_sites[1])
                    .map(|x| x.clone())
                    .collect::<Vec<Interval<u64, GeneInterval>>>();

                if result.is_empty() {
                    return None;
                } else if result.len() != 1 {
                    return None;
                } else {
                    // only process the first match
                    let gene_interval = &result[0].val;
                    let matched_sites =
                        gene_interval.match_splice_sites(&splice_sites.to_vec(), flank);
                    if matched_sites.is_empty() {
                        return None;
                    } else {
                        return Some((chrom, matched_sites));
                    }
                }
            }
            None => None,
        }
    }
}
