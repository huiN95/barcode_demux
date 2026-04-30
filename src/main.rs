// mod alignment_statistics;
mod bam_record_extention;
mod cli;
mod find_pattern;
mod get_demuxed_reads;
mod io_utils;
mod pbar;
mod reader_worker;
mod run_logger;
// mod statitcs;
mod barcode_demuxer_v1;
mod barcode_mysers;
mod core_context;
mod demux_pipeline_v1;
mod demux_pipeline_v2;
mod demux_pipeline_v3;
mod demux_pipeline_v4;
mod demux_primer;
// mod multiple_barcode_demuxer_v1;
mod writer_worker;
use crate::demux_pipeline_v1::demux_pipeline_v1;
use crate::demux_pipeline_v2::demux_pipeline_v2;
use crate::demux_pipeline_v3::demux_pipeline_v3;
use crate::demux_pipeline_v4::demux_pipeline_v4;
use clap::Parser;
use cli::Cli;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use run_logger::{init_metrics, init_tracing_log};
// use crate::io_utils::ensure_output_dir;
use tracing::info;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let log_path = &cli.log_folder;

    let run_log_prefix = make_run_log_prefix(&cli);

    let _log_writer = init_tracing_log(&cli, &run_log_prefix);
    let _metric_writer = init_metrics(log_path, &run_log_prefix)?;

    info!("Start processing");
    info!(?cli, "parsed CLI");
    info!(run_log_prefix = %run_log_prefix, "run log prefix");

    match cli.pipeline_version.as_str() {
        "1" => {
            info!("artificial barcode");
            demux_pipeline_v1(&cli).unwrap();
        }
        "3" => {
            info!("双端打分平分过滤");
            demux_pipeline_v3(&cli).unwrap();
        }
        "4" => {
            info!("删除末端序列");
            demux_pipeline_v4(&cli).unwrap();
        }
        other => panic!("Unsupported pipeline version: {}", other),
    }

    info!("End processing");

    Ok(())
}

fn make_run_log_prefix(cli: &Cli) -> String {
    let input_stem = Path::new(&cli.input_file)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("input");

    let input_stem = sanitize_filename(input_stem);

    let ts_ns = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();

    let pid = std::process::id();

    format!(
        "{}_pv{}_pid{}_{}",
        input_stem, cli.pipeline_version, pid, ts_ns
    )
}

fn sanitize_filename(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

#[cfg(test)]
mod test {
    use clap::Parser;

    use super::*;

    #[test]
    fn test_bam_format() {
        let args = vec![
            "my_program",
            "--input_file",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/third.smc_all_reads.bam",
            "-o",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/art",
            "--log_folder",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/art",
            "--output_format",
            "bam",
            "--min_q",
            "20",
            "--search_bound",
            "50",
            "-b",
            "/mnt/data_7t/adam/adapter_bcs/fake_barcode.fasta",
        ];
        let matches = Cli::parse_from(args);
        let result = demux_pipeline_v1(&matches);
        assert!(result.is_ok(), "Pipeline failed: {:?}", result.err());
    }

    #[test]
    fn test_tail_cutoff() {
        let args = vec![
            "my_program",
            "--input_file",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/third.smc_all_reads.bam",
            "-o",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/art",
            "--log_folder",
            "/mnt/data_7t/adam/barcode_demux/20260205_250302Y0004_Run0006/art",
            "--output_format",
            "bam",
            "--min_q",
            "20",
            "--search_bound",
            "50",
            "-b",
            "/mnt/data_7t/adam/adapter_bcs/fake_barcode.fasta",
        ];
        let matches = Cli::parse_from(args);
        let result = demux_pipeline_v1(&matches);
        assert!(result.is_ok(), "Pipeline failed: {:?}", result.err());
    }

    #[test]
    fn test_tail_cutoff_rna() {
        let args = vec![
            "my_program",
            "--input_file",
            "/mnt/data_7t/adam/rna_data/RNA2K/ch_201356.bam",
            "-o",
            "/mnt/data_7t/adam/rna_data/RNA2K/barcode",
            "--log_folder",
            "/mnt/data_7t/adam/rna_data/RNA2K",
            "--output_format",
            "fq",
            "--search_bound",
            "200",
            "--barcode",
            "/mnt/data_7t/adam/rna_data/polyt.fasta",
            "--max_distance",
            "3",
            "--min_pair_len",
            "20",
            "--min_pair_score",
            "0.9",
            "--pipeline_version",
            "4",
        ];
        let matches = Cli::parse_from(args);
        let result = demux_pipeline_v4(&matches);
        assert!(result.is_ok(), "Pipeline failed: {:?}", result.err());
    }
}
