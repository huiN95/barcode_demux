use crate::core_context::BarcodePair;
use crate::get_demuxed_reads::RecordType;
use crate::io_utils::{make_writers, normalize_barcode_name};

use crate::pbar::{get_spin_pb, DEFAULT_INTERVAL};
// use clap::error::ContextKind;
use crossbeam::channel::Receiver;
use tracing::error;
// use flate2::write;
use metrics::counter;
// use std::error::Error;
// use std::path::Path;

pub fn write_barcode_results(
    input_file: &str,
    output_format: &str,
    output_folder: &str,
    barcode_demux_info: Receiver<(RecordType, Option<BarcodePair>)>,

    barcode_names: Vec<String>,
) -> anyhow::Result<()> {
    print!("enter write process");
    // let path: &Path = output_folder.as_ref();

    // let stem = path.with_extension(""); // 去掉 .fasta / .fastq / .bam

    // let distance_scale: u16 = 10;
    // static r#EPS: f32 = 0.001;

    let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);
    let mut writers = make_writers(output_folder, barcode_names, input_file, output_format)?;
    match output_format {
        "fa" | "fasta" | "fastq" | "fq" | "bam" => {
            while let Ok(demuxed_record) = barcode_demux_info.recv() {
                counter!("writer_received").increment(1 as u64);

                let (demuxed_reads, barcode_pair) = demuxed_record;
                let Some(barcode_pair) = barcode_pair else {
                    counter!("barcode_not_found_channels").increment(1 as u64); // 如果两个都没有，跳过
                    let w = writers
                        .get_mut("uncertain")
                        .ok_or_else(|| anyhow::anyhow!("writers missing key: uncertain"))?;
                    w.write_demuxed_record(&demuxed_reads)?;
                    // writers
                    //     .get_mut("uncertain")
                    //     .unwrap()
                    //     .write_demuxed_record(&demuxed_reads)?;
                    continue;
                };
                let key = barcode_pair.name.as_ref();
                let w = writers
                    .get_mut(key)
                    .ok_or_else(|| anyhow::anyhow!("writers missing key: {key}"))?;
                if let Err(e) = w.write_demuxed_record(&demuxed_reads) {
                    error!(?e, "有效subreads文件写入失败");
                }
                // if let Err(e) = writers
                //     .get_mut(barcode_pair.name.as_ref())
                //     .unwrap()
                //     .write_demuxed_record(&demuxed_reads)
                // {
                //     error!(?e, "有效subreads文件写入失败");
                // }
                counter!("wirter_ok").increment(1 as u64);

                pb.inc(1);
            }
        }
        _ => anyhow::bail!("不支持的写入格式: {output_format}"),
    }

    Ok(())
}

pub fn write_residual_results(
    input_file: &str,
    output_format: &str,
    output_folder: &str,
    barcode_demux_info: Receiver<(RecordType, Option<BarcodePair>)>,

    barcode_names: Vec<String>,
) -> anyhow::Result<()> {
    print!("enter write process");
    // let path: &Path = output_folder.as_ref();

    // let stem = path.with_extension(""); // 去掉 .fasta / .fastq / .bam

    // let distance_scale: u16 = 10;
    // static r#EPS: f32 = 0.001;

    let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);
    let mut writers = make_writers(output_folder, barcode_names, input_file, output_format)?;
    match output_format {
        "fa" | "fasta" | "fastq" | "fq" | "bam" => {
            while let Ok(demuxed_record) = barcode_demux_info.recv() {
                counter!("writer_received").increment(1 as u64);

                let (demuxed_reads, barcode_pair) = demuxed_record;
                let Some(barcode_pair) = barcode_pair else {
                    counter!("writer_uncertain").increment(1 as u64); // 如果两个都没有，跳过
                    let w = writers
                        .get_mut("uncertain")
                        .ok_or_else(|| anyhow::anyhow!("writers missing key: uncertain"))?;
                    w.write_demuxed_record(&demuxed_reads)?;
                    // writers
                    //     .get_mut("uncertain")
                    //     .unwrap()
                    //     .write_demuxed_record(&demuxed_reads)?;
                    continue;
                };
                let key = barcode_pair.name.as_ref();
                let w = writers
                    .get_mut(key)
                    .ok_or_else(|| anyhow::anyhow!("writers missing key: {key}"))?;
                if let Err(e) = w.write_demuxed_record(&demuxed_reads) {
                    error!(?e, "有效subreads文件写入失败");
                }
                // if let Err(e) = writers
                //     .get_mut(barcode_pair.name.as_ref())
                //     .unwrap()
                //     .write_demuxed_record(&demuxed_reads)
                // {
                //     error!(?e, "有效subreads文件写入失败");
                // }
                counter!("wirter_ok").increment(1 as u64);

                pb.inc(1);
            }
        }
        _ => anyhow::bail!("不支持的写入格式: {output_format}"),
    }

    Ok(())
}
