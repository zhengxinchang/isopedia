use std::fmt::Display;

use bio_types::strand::ReqStrand;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

impl Strand {
    pub fn from_string(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    #[allow(unused)]
    pub fn from_num_string(s: &str) -> Self {
        match s {
            "0" => Strand::Plus,
            "1" => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    #[allow(unused)]
    pub fn from_u8(n: u8) -> Self {
        match n {
            0 => Strand::Plus,
            1 => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            Strand::Plus => "+".to_string(),
            Strand::Minus => "-".to_string(),
        }
    }

    pub fn from_bam_strand(s: ReqStrand) -> Self {
        match s {
            ReqStrand::Forward => Strand::Plus,
            ReqStrand::Reverse => Strand::Minus,
        }
    }
}
