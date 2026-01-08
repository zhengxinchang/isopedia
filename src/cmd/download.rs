//! download remote index files

use std::path::PathBuf;

use clap::Parser;
use log::error;
use serde::Deserialize;
use serde::Serialize;

use crate::utils::greetings2;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia download")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[Download indexes] Download pre-built index files.
", long_about = None)]
#[clap(after_long_help = "
Example: 

# (step 1) fetch pre-built index manifest file
isopedia download 

# (step 2) download pre-built index files to the specified directory
isopedia download -n demo -o /path/to/output

")]

pub struct DownloadCli {
    /// List of available pre-built index files
    #[arg(short, long)]
    pub list: bool,

    /// Name of the index files to download from official repository
    #[arg(short, long)]
    pub name: Option<String>,

    /// Download a custom URL index file
    #[arg(short, long)]
    pub url: Option<String>,

    /// Download index based on a custom manifest file. This is mutually exclusive with --name
    #[arg(short, long)]
    pub manifest: Option<PathBuf>,

    /// Output directory to save the downloaded index files
    #[arg(short, long, default_value = ".")]
    pub outdir: Option<PathBuf>,
}

impl DownloadCli {
    fn validate(&self) {
        let mut is_ok = true;

        if self.name.is_none() && self.url.is_none() && self.manifest.is_none() && !self.list {
            error!("Either --name, --url, --manifest or --list must be provided.");
            is_ok = false;
        }

        if !is_ok {
            std::process::exit(1);
        }
    }
}

pub fn run_download() {
    println!("Downloading index files...");

    let cli = DownloadCli::parse();
    greetings2(&cli);
    cli.validate();
}

#[derive(Debug, Serialize, Deserialize)]
pub struct IndexItem {
    pub name: String,
    pub url: String,
    pub size: u64,

    pub description: Option<String>,

    pub md5: Option<String>,
    pub sha256: Option<String>,

    #[serde(default)]
    pub source: SourceType,

    #[serde(default)]
    pub protocol: Protocol,
}

#[derive(Debug, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum SourceType {
    #[default]
    Direct,
    S3,
    GCS, // Google Cloud Storage
    Azure,
    FTP,
    Custom(String),
}

#[derive(Debug, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum Protocol {
    HTTP,
    #[default]
    HTTPS,
    FTP,
    S3,
}
