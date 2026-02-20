//! download remote index files

use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Write},
    path::PathBuf,
    process::Command,
    time::Duration,
};

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use log::{error, info};
use md5::Context as Md5Context;
use reqwest::blocking::Client;
use serde::Deserialize;
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::constants::RMT_MANIFEST_URL1;
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
        let selected_sources =
            self.name.is_some() as u8 + self.url.is_some() as u8 + self.manifest.is_some() as u8;

        if self.list && selected_sources > 0 {
            error!("--list cannot be combined with --name, --url, or --manifest.");
            is_ok = false;
        }

        if !self.list && self.url.is_none() && self.name.is_none() {
            error!("Provide --list, or provide --name / --url.");
            is_ok = false;
        }

        if self.name.is_some() && self.url.is_some() {
            error!("--name and --url are mutually exclusive.");
            is_ok = false;
        }

        if let Some(manifest_path) = &self.manifest {
            if !manifest_path.exists() {
                error!(
                    "--manifest file does not exist: {}",
                    manifest_path.display()
                );
                is_ok = false;
            }
        }

        if let Some(outdir) = &self.outdir {
            if !outdir.exists() {
                if let Err(e) = fs::create_dir_all(outdir) {
                    error!("Failed to create outdir {}: {}", outdir.display(), e);
                    is_ok = false;
                }
            }
            if !outdir.is_dir() {
                error!("--outdir is not a directory: {}", outdir.display());
                is_ok = false;
            }
        } else {
            error!("--outdir is required.");
            is_ok = false;
        }

        if !is_ok {
            std::process::exit(1);
        }
    }
}

pub fn run_download(cli: &DownloadCli) -> Result<()> {
    greetings2(cli);
    cli.validate();

    if cli.list {
        let manifest = load_manifest(cli.manifest.as_ref())?;
        print_manifest(&manifest.index);
        return Ok(());
    }

    let item = if let Some(name) = &cli.name {
        let manifest = load_manifest(cli.manifest.as_ref())?;
        manifest
            .index
            .into_iter()
            .find(|it| &it.name == name)
            .ok_or_else(|| anyhow!("Index name not found in manifest: {}", name))?
    } else if let Some(url) = &cli.url {
        IndexItem {
            name: file_name_from_url(url),
            url: url.to_string(),
            size: 0,
            description: Some("custom-url".to_string()),
            md5: None,
            sha256: None,
            source: SourceType::Direct,
            protocol: Protocol::HTTPS,
        }
    } else {
        return Err(anyhow!("No download source selected"));
    };

    let outdir = cli.outdir.as_ref().expect("outdir validated");
    let outpath = outdir.join(&item.name);
    download_item(&item, &outpath)?;
    info!("Downloaded to {}", outpath.display());
    Ok(())
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

#[derive(Debug, Serialize, Deserialize)]
struct IndexManifest {
    index: Vec<IndexItem>,
}

fn load_manifest(manifest_path: Option<&PathBuf>) -> Result<IndexManifest> {
    if let Some(path) = manifest_path {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read manifest {}", path.display()))?;
        let parsed: IndexManifest = toml::from_str(&content)
            .with_context(|| format!("Failed to parse manifest {}", path.display()))?;
        return Ok(parsed);
    }

    let client = Client::builder()
        .timeout(Duration::from_secs(60))
        .build()
        .context("Failed to build HTTP client")?;
    let response = client
        .get(RMT_MANIFEST_URL1)
        .send()
        .with_context(|| format!("Failed to fetch remote manifest {}", RMT_MANIFEST_URL1))?;
    if !response.status().is_success() {
        return Err(anyhow!(
            "Remote manifest request failed with status {}",
            response.status()
        ));
    }
    let content = response
        .text()
        .context("Failed to decode remote manifest content")?;
    let parsed: IndexManifest =
        toml::from_str(&content).context("Failed to parse remote manifest TOML")?;
    Ok(parsed)
}

fn print_manifest(items: &[IndexItem]) {
    println!("Available index entries: {}", items.len());
    for item in items {
        let desc = item.description.clone().unwrap_or_else(|| "-".to_string());
        println!(
            "name={}\n  url={}\n  size={}\n  description={}\n",
            item.name, item.url, item.size, desc
        );
    }
}

fn file_name_from_url(url: &str) -> String {
    let no_query = url.split('?').next().unwrap_or(url);
    no_query
        .trim_end_matches('/')
        .rsplit('/')
        .next()
        .filter(|x| !x.is_empty())
        .unwrap_or("downloaded.index")
        .to_string()
}

fn download_item(item: &IndexItem, outpath: &PathBuf) -> Result<()> {
    info!("Downloading {} from {}", item.name, item.url);
    if item.url.to_ascii_lowercase().starts_with("ftp://") {
        download_ftp_with_curl(&item.url, outpath)?;
        validate_download(item, outpath)?;
        info!("Download and validation complete for {}", item.name);
        return Ok(());
    }

    download_http(&item.url, outpath)?;
    validate_download(item, outpath)?;
    info!("Download and validation complete for {}", item.name);
    Ok(())
}

fn download_http(url: &str, outpath: &PathBuf) -> Result<()> {
    let client = Client::builder()
        .timeout(Duration::from_secs(600))
        .build()
        .context("Failed to build HTTP client")?;
    let mut response = client
        .get(url)
        .send()
        .with_context(|| format!("Failed to download {}", url))?;
    if !response.status().is_success() {
        return Err(anyhow!(
            "Download request failed with status {}",
            response.status()
        ));
    }

    let file = File::create(outpath)
        .with_context(|| format!("Failed to create output file {}", outpath.display()))?;
    let mut writer = BufWriter::new(file);

    let mut buf = vec![0u8; 8 * 1024 * 1024];
    let mut total_written: u64 = 0;
    let mut next_log_mark: u64 = 100 * 1024 * 1024;

    loop {
        let n = response
            .read(&mut buf)
            .context("Failed while reading response body")?;
        if n == 0 {
            break;
        }
        let chunk = &buf[..n];
        writer
            .write_all(chunk)
            .context("Failed while writing downloaded bytes")?;
        total_written += n as u64;

        if total_written >= next_log_mark {
            info!("Downloaded {} bytes...", total_written);
            next_log_mark += 100 * 1024 * 1024;
        }
    }
    writer.flush().context("Failed to flush output file")?;
    Ok(())
}

fn download_ftp_with_curl(url: &str, outpath: &PathBuf) -> Result<()> {
    let status = Command::new("curl")
        .arg("-fL")
        .arg(url)
        .arg("-o")
        .arg(outpath)
        .status()
        .with_context(|| format!("Failed to spawn curl for {}", url))?;
    if !status.success() {
        return Err(anyhow!("curl download failed for {}", url));
    }
    Ok(())
}

fn validate_download(item: &IndexItem, outpath: &PathBuf) -> Result<()> {
    let (total_written, sha256_hex, md5_hex) = compute_file_hashes(outpath)?;
    if item.size > 0 && total_written != item.size {
        return Err(anyhow!(
            "Size mismatch for {}: expected {}, got {}",
            item.name,
            item.size,
            total_written
        ));
    }

    if let Some(expected_sha256) = &item.sha256 {
        if sha256_hex.to_lowercase() != expected_sha256.to_lowercase() {
            return Err(anyhow!(
                "SHA256 mismatch for {}: expected {}, got {}",
                item.name,
                expected_sha256,
                sha256_hex
            ));
        }
    }

    if let Some(expected_md5) = &item.md5 {
        if md5_hex.to_lowercase() != expected_md5.to_lowercase() {
            return Err(anyhow!(
                "MD5 mismatch for {}: expected {}, got {}",
                item.name,
                expected_md5,
                md5_hex
            ));
        }
    }

    Ok(())
}

fn compute_file_hashes(path: &PathBuf) -> Result<(u64, String, String)> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut sha256_hasher = Sha256::new();
    let mut md5_hasher = Md5Context::new();
    let mut buf = vec![0u8; 8 * 1024 * 1024];
    let mut total = 0u64;

    loop {
        let n = reader
            .read(&mut buf)
            .context("Failed while reading downloaded file")?;
        if n == 0 {
            break;
        }
        let chunk = &buf[..n];
        sha256_hasher.update(chunk);
        md5_hasher.consume(chunk);
        total += n as u64;
    }

    let sha256_hex = format!("{:x}", sha256_hasher.finalize());
    let md5_hex = format!("{:x}", md5_hasher.finalize());
    Ok((total, sha256_hex, md5_hex))
}
