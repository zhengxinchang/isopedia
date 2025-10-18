use anyhow::Result;
use clap::Parser;
use log::error;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::output::{GeneralTableOutput, IsoformTableOut};

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "
Merge multiple output files into one.
Examples:
isopedia-tool merge-output --input <INPUT> --output <OUTPUT> --mode <MODE>

MODE can be:
- 'merge-isoform': merge isoform outputs
- 'merge-fusionbreakpoint': merge fusion breakpoint outputs
- 'merge-fusiondiscovery': merge fusion discovery outputs
- 'merge-splice': merge splice outputs
")]
pub struct OutputArg {
    #[arg(short, long)]
    pub input: Vec<PathBuf>,
    #[arg(short, long)]
    pub output: PathBuf,
    #[arg(short, long)]
    pub mode: String,
}

impl OutputArg {
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if self.input.len() < 2 {
            error!("--input: at least two input files must be specified");
            is_ok = false;
        }

        let valid_modes = vec![
            "merge-isoform",
            "merge-fusionbreakpoint",
            "merge-fusiondiscovery",
            "merge-splice",
        ];
        if !valid_modes.contains(&self.mode.as_str()) {
            error!(
                "--mode: invalid mode '{}'. Valid modes are: {:?}",
                self.mode, valid_modes
            );
            is_ok = false;
        }

        is_ok
    }
}

pub fn merge(cli: &OutputArg) -> Result<()> {
    let base = cli.input[0].clone();
    match cli.mode.as_str() {
        "merge-isoform" => {
            let mut base_out = IsoformTableOut::load(&base)?;
            for other in &cli.input[1..] {
                let other_out = IsoformTableOut::load(other)?;
                base_out.merge(&other_out)?;
            }

            base_out.save_to_file(&cli.output)?;
        }
        _ => {}
    }

    Ok(())
}
