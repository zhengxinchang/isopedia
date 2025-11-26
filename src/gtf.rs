use crate::{reads::SingleRead, utils::trim_chr_prefix_to_upper};
use bio_types::strand::ReqStrand;
use core::panic;
use flate2::read::MultiGzDecoder;
use noodles_gtf::{self as GTF, io::Reader as gtfReader, record::Strand};
use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Seek},
};

// need pack into a u64
pub struct TranscriptMini {
    // pub idx: u32, // 27bit range: 0 ~ 134,217,727
    // pub is_mono_exonic: bool, // 1bit range 0 ~ 1
    // pub n_sj: u32, //10bit range 0 ~ 1023
    // pub seq_len: u32, // 26bit range 0 ~ 67,108,863
    pub bits: u64,
    // pub sj_pairs: Vec<(u64,u64)>,
}

impl TranscriptMini {
    pub fn from_transcript(tx: &Transcript, tx_idx: usize) -> Self {
        let mut bits: u64 = 0;
        bits |= (tx_idx as u64) & 0x07FF_FFFF; // 27 bits for idx
        bits <<= 1;
        bits |= if tx.is_mono_exonic { 1 } else { 0 }; // 1 bit for is_mono_exonic
        bits <<= 10;
        bits |= (tx.splice_junc.len() as u64) & 0x03FF; // 10 bits for n_sj
        bits <<= 26;
        bits |= (tx.get_transcript_seq_length() as u64) & 0x03FF_FFFF; // 26 bits for seq_len

        TranscriptMini {
            // idx: tx_idx as u32,
            // is_mono_exonic: tx.is_mono_exonic,
            // n_sj: tx.splice_junc.len() as u32,
            // seq_len: tx.get_transcript_seq_length() as u32,
            bits,
        }
    }

    pub fn get_transcript_idx(&self) -> u32 {
        ((self.bits >> 37) & 0x07FF_FFFF) as u32
    }

    pub fn get_is_mono_exonic(&self) -> bool {
        ((self.bits >> 36) & 0x01) != 0
    }

    pub fn get_n_sj(&self) -> u32 {
        ((self.bits >> 26) & 0x03FF) as u32
    }

    pub fn get_seq_len(&self) -> u32 {
        (self.bits & 0x03FF_FFFF) as u32
    }
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub splice_junc: Vec<(u64, u64)>,
    pub exons: Vec<(u64, u64)>,
    pub is_mono_exonic: bool,
    pub gene_id: String,
    pub trans_id: String,
    pub records: Vec<GTF::Record>,
    pub strand: Strand,
}

impl Transcript {
    pub fn process(&mut self) -> &mut Self {
        // guard for emtpy exon length
        if self.exons.len() == 0 {
            panic!("Transcript (id={}) must have at least one exon. Check your if input GTF is sorted...",self.trans_id);
        }

        self.exons.sort_by_key(|k| k.0);

        if self.exons.len() >= 2 {
            self.splice_junc = self.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        } else {
            self.splice_junc = vec![(self.exons[0].0, self.exons[0].1)];
        }
        self.is_mono_exonic = self.exons.len() == 1;
        self
    }

    pub fn from_gtf_record(record: &GTF::Record) -> Self {
        let mut trans = Transcript {
            chrom: String::new(),
            splice_junc: Vec::new(),
            start: 0,
            end: 0,
            exons: Vec::new(),
            is_mono_exonic: false,
            gene_id: String::new(),
            trans_id: String::new(),
            records: Vec::new(),
            strand: record.strand().unwrap_or_else(||{
                panic!("GTF must have strand information for transcript/exon records, however this record is missing.\nRecord: {record:?}")
            }),
        };

        trans.chrom = trim_chr_prefix_to_upper(record.reference_sequence_name());
        trans.gene_id = record
            .attributes()
            .iter()
            .find(|x| x.key() == "gene_id")
            .unwrap_or_else(|| panic!("GTF must have gene id"))
            .value()
            .to_string();
        trans.trans_id = record
            .attributes()
            .iter()
            .find(|x| x.key() == "transcript_id")
            .unwrap_or_else(|| panic!("GTF must have transcript_id"))
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

    pub fn get_transcript_coord_length(&self) -> u64 {
        self.end - self.start
    }

    pub fn get_transcript_seq_length(&self) -> u32 {
        let mut total_len: u64 = 0;
        for exon in &self.exons {
            total_len += exon.1 - exon.0;
        }
        total_len as u32
    }

    /// Convert Transcript to SingleRead
    /// Used in indexing GTF files instead of parsing reads from BAM files
    pub fn to_single_read(&self) -> SingleRead {
        let mut sr = SingleRead::new(
            self.chrom.clone(),
            60,
            1,
            self.get_transcript_coord_length(),
            self.start,
        );

        let strand = match self.strand {
            Strand::Forward => ReqStrand::Forward,
            Strand::Reverse => ReqStrand::Reverse,
        };

        for exon in &self.exons {
            sr.add_segment(self.chrom.clone(), exon.0, exon.1, &strand, false);
        }
        sr.update_right(self.end);
        sr.process();

        sr
    }

    pub fn get_splice_junction_pairs(&self) -> Vec<(u64, u64)> {
        self.splice_junc.clone()
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
        // let mut record_no = 0;
        loop {
            let record = records.next();
            // record_no += 1;
            self.trans_count += 1;
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
                        panic!("error in reading gtf: {:?}, record (no.{}: {:?}\nPlease check your gtf format. One possible reason is the input file is GFF format.", e, self.trans_count, &record);
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

    pub fn get_all_transcripts_vec(&mut self) -> Vec<Transcript> {
        let mut transcripts = Vec::new();
        while let Some(trans) = self.get_next_transcript() {
            transcripts.push(trans);
        }
        transcripts.sort_by(|a, b| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
        transcripts
    }

    pub fn get_all_transcripts_by_chrom(&mut self) -> Vec<(String, Vec<Transcript>)> {
        let mut transcripts_by_chrom: std::collections::HashMap<String, Vec<Transcript>> =
            std::collections::HashMap::new();
        while let Some(trans) = self.get_next_transcript() {
            transcripts_by_chrom
                .entry(trans.chrom.clone())
                .or_insert_with(Vec::new)
                .push(trans);
        }
        let mut transcripts_by_chrom_vec: Vec<(String, Vec<Transcript>)> =
            transcripts_by_chrom.into_iter().collect();
        transcripts_by_chrom_vec.sort_by(|a, b| a.0.cmp(&b.0));

        // sort transcripts in each chrom by start position
        for (_chrom, txs) in transcripts_by_chrom_vec.iter_mut() {
            txs.sort_by(|a, b| a.start.cmp(&b.start));
        }

        transcripts_by_chrom_vec
    }
}

impl<R: BufRead> Iterator for TranscriptChunker<R> {
    type Item = Transcript;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_transcript()
    }
}

pub fn open_gtf_reader(path: &str) -> std::io::Result<noodles_gtf::io::Reader<Box<dyn BufRead>>> {
    let file = File::open(path)?;
    let mut buf = [0u8; 2];
    let mut peek = BufReader::new(file);
    let n = peek.read(&mut buf)?;

    let mut file = peek.into_inner();
    file.rewind()?;

    let reader: Box<dyn BufRead> = if n == 2 && buf == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(noodles_gtf::io::Reader::new(reader))
}
