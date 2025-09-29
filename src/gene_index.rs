use std::io::BufRead;

use anyhow::Result;
use indexmap::{IndexMap, IndexSet};
use noodles_gtf::Reader as gtfReader;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;

use crate::utils::trim_chr_prefix_to_upper;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneType {
    ProteinCoding,
    NonCoding,
    Pseudogene,
    MiRna,
    Other,
}

#[derive(Debug, Clone)]
pub struct GeneInterval {
    pub gene_id: String,
    pub gene_name: String,
    pub gene_type: GeneType,
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
    pub fn default() -> Self {
        GeneInterval {
            gene_id: String::new(),
            gene_name: String::new(),
            gene_type: GeneType::Other,
            chrom: String::new(),
            start: 0,
            end: 0,
            splice_sites: Vec::new(),
        }
    }

    pub fn new(
        gene_id: String,
        gene_name: String,
        gene_type: GeneType,
        chrom: String,
        start: u64,
        end: u64,
    ) -> Self {
        GeneInterval {
            gene_id,
            gene_name,
            gene_type,
            chrom,
            start,
            end,
            splice_sites: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.gene_id.clear();
        self.gene_name.clear();
        self.chrom.clear();
        self.start = 0;
        self.end = 0;
        self.splice_sites.clear();
    }

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
#[derive(Debug, Clone)]
pub struct GeneIntervalTree {
    pub tree: HashMap<String, Lapper<u64, GeneInterval>>,
    pub count: usize,
    pub chroms: Vec<String>,
}

impl GeneIntervalTree {
    pub fn new<R: BufRead>(gtf_reader: &mut gtfReader<R>) -> Result<Self> {
        let mut chrom_tree: IndexMap<String, Vec<GeneInterval>> = IndexMap::new();

        let records = gtf_reader.records();
        let mut splice_sites = Vec::new();

        let mut interval_count = 0;
        let mut chrom_set = IndexSet::new();
        let mut current_gene = GeneInterval::default();

        for record in records {
            let record = record?;
            let chrom = trim_chr_prefix_to_upper(record.reference_sequence_name());
            chrom_set.insert(chrom.clone());
            chrom_tree
                .entry(chrom.clone())
                .or_insert_with(|| Vec::<GeneInterval>::new());
            match record.ty() {
                "gene" => {
                    if splice_sites.is_empty() {
                        current_gene.chrom = chrom.clone();
                        current_gene.gene_id = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_id")
                            .expect("GTF must have gene_id")
                            .value()
                            .to_string();
                        current_gene.gene_name = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_name")
                            .expect("GTF must have gene_name")
                            .value()
                            .to_string();

                        current_gene.gene_type = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_type")
                            .map(|x| match x.value() {
                                "protein_coding" => GeneType::ProteinCoding,
                                "lncRNA" => GeneType::NonCoding,
                                "transcribed_unprocessed_pseudogene" => GeneType::Pseudogene,
                                _ => GeneType::Other,
                            })
                            .unwrap_or(GeneType::Other);

                        current_gene.start = record.start().get() as u64 - 1;
                        current_gene.end = record.end().get() as u64;
                    } else {
                        // make new GeneInterval and reset splice_sites

                        splice_sites.sort();
                        splice_sites.dedup();

                        current_gene.splice_sites = splice_sites.clone();

                        // dbg!(&current_gene.gene_name);

                        let gene_intervals = chrom_tree.get_mut(&current_gene.chrom).unwrap();
                        interval_count += 1;
                        gene_intervals.push(current_gene.clone());

                        // reload & reset
                        splice_sites.clear();
                        current_gene.clear();

                        current_gene.chrom = chrom.clone();
                        current_gene.gene_id = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_id")
                            .expect("GTF must have gene_id")
                            .value()
                            .to_string();
                        current_gene.gene_name = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_name")
                            .expect("GTF must have gene_name")
                            .value()
                            .to_string();
                        current_gene.gene_type = record
                            .attributes()
                            .iter()
                            .find(|x| x.key() == "gene_type")
                            .map(|x| match x.value() {
                                "protein_coding" => GeneType::ProteinCoding,
                                "lncRNA" => GeneType::NonCoding,
                                "transcribed_unprocessed_pseudogene" => GeneType::Pseudogene,
                                _ => GeneType::Other,
                            })
                            .unwrap_or(GeneType::Other);
                        current_gene.start = record.start().get() as u64 - 1;
                        current_gene.end = record.end().get() as u64;
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

            current_gene.splice_sites = splice_sites.clone();

            let gene_intervals = chrom_tree.get_mut(&current_gene.chrom).unwrap();
            interval_count += 1;
            gene_intervals.push(current_gene.clone());
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

    /// match splice sites between fusion candidates and gene intervals.
    pub fn match2(
        &self,
        chrom: &str,
        splice_sites: &Vec<u64>,
        flank: u64,
    ) -> Option<(GeneInterval, Vec<u64>, CandidateMatchStatus)> {
        let chrom = trim_chr_prefix_to_upper(chrom);
        match self.tree.get(&chrom) {
            Some(lapper) => {
                let result = lapper
                    .find(splice_sites[0], *splice_sites.last().unwrap())
                    .map(|x| x.clone())
                    .collect::<Vec<Interval<u64, GeneInterval>>>();

                if result.is_empty() {
                    // skip if no matches
                    // log::debug!("No matches found for splice sites: {:?}", splice_sites);
                    return None;
                } else {
                    let mut potential_maches = Vec::new();

                    for gene in result.iter() {
                        if gene.val.gene_type != GeneType::ProteinCoding {
                            // debug!("Skipping non-protein coding gene: {}", gene.val.gene_name);
                            continue;
                        } else {
                            // make sure splice sites match
                            let matched_sites = gene.val.match_splice_sites(splice_sites, flank);
                            if !matched_sites.is_empty() {
                                potential_maches.push((gene.val.clone(), matched_sites));
                            }
                        }
                    }

                    if potential_maches.is_empty() {
                        // debug!(
                        //     "No potential matches found for splice sites: {:?}",
                        //     splice_sites
                        // );
                        return None;
                    } else {
                        potential_maches.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

                        if potential_maches.len() == 1 {
                            // unique match
                            let (gene, splice_sites) = potential_maches[0].clone();
                            // debug!("Unique match found: {}", gene.gene_name);
                            return Some((gene, splice_sites, CandidateMatchStatus::Unique));
                        } else {
                            return Some((
                                potential_maches[0].0.clone(),
                                potential_maches[0].1.clone(),
                                CandidateMatchStatus::Ambiguous,
                            ));
                        }
                    }

                    // sort by the number of matched splice sites
                }
            }
            None => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateMatchStatus {
    Unique,
    Ambiguous,
}
