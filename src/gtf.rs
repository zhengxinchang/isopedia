use crate::utils::trim_chr_prefix_to_upper;
use core::panic;
use noodles_gtf::{self as GTF, io::Reader as gtfReader, record::Strand};
use std::io::BufRead;

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
            .unwrap()
            .value()
            .to_string();
        trans.trans_id = record
            .attributes()
            .iter()
            .find(|x| x.key() == "transcript_id")
            .unwrap()
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
            match record {
                Some(rec) => match rec {
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
                        panic!("error in reading gtf: {:?}", e);
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
