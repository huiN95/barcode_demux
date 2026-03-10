use crate::bam_record_extention::ReadRecord;
use crate::core_context::BarcodePair;

use anyhow::{anyhow, Result};
use bio::io::{fasta, fastq};
use metrics::{counter, histogram};
use rust_htslib::bam::{
    self,
    record::{Aux, Record},
};
// use std::cmp::min;
// use std::{ops::Deref, path::Path};
// use tracing_subscriber::layer::SubscriberExt;
pub enum RecordType {
    Fastq(fastq::Record),
    Fasta(fasta::Record),
    BAM(Record),
}

pub fn prepare_record_to_writer(
    rec: &ReadRecord,
    barcode_paris: Option<BarcodePair>,
    min_subread_len: usize,
    output_format: &str,
    q_threshold: usize,
) -> Result<(RecordType, Option<BarcodePair>)> {
    // let primer_records = get_primer_record(rec, &barcode_paris, output_file)?;
    let subread_records = get_subreads_from_demuxed_result(
        rec,
        barcode_paris,
        min_subread_len,
        // None, // 质量阈值
        q_threshold,
        output_format,
    )?;

    Ok(subread_records)
    // Ok((primer_records, subread_records))
}

fn make_name(rec_id: &str, primer_id: &str, s: usize, e: usize) -> String {
    // 预估容量：id+primer+数字
    let mut name = String::with_capacity(rec_id.len() + 32);
    use std::fmt::Write;
    write!(&mut name, "{}/{}-{}", rec_id, s, e).unwrap();
    name
}

fn get_subreads_from_demuxed_result(
    rec: &ReadRecord,
    barcode_position: Option<BarcodePair>,
    min_subread_len: usize,
    q_threshold: usize,
    output_format: &str,
) -> anyhow::Result<(RecordType, Option<BarcodePair>)> {
    match output_format {
        // ---------- FASTA ----------
        "fa" | "fasta" => {
            let (name, start, end, bc_opt) = match barcode_position {
                Some(bc) => {
                    let (s, e) = bc.inner_position;
                    if e.saturating_sub(s) >= min_subread_len {
                        (make_name(&rec.id, bc.name.as_ref(), s, e), s, e, Some(bc))
                    } else {
                        counter!("reads_too_short").increment(1);
                        (
                            make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                            0,
                            rec.sequence.len(),
                            None,
                        )
                    }
                }
                None => (
                    make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                    0,
                    rec.sequence.len(),
                    None,
                ),
            };

            let fasta_rec = fasta::Record::with_attrs(&name, None, &rec.sequence[start..end]);
            return Ok((RecordType::Fasta(fasta_rec), bc_opt));
        }
        // ---------- FASTQ ----------
        "fq" | "fastq" => {
            let (name, start, end, bc_opt) = match barcode_position {
                Some(bc) => {
                    let (s, e) = bc.inner_position;
                    let len_ok = e.saturating_sub(s) >= min_subread_len;

                    let q_ok = if len_ok {
                        let q_slice = rec
                            .quality
                            .as_deref()
                            .ok_or_else(|| anyhow::anyhow!("missing quality"))?
                            .get(s..e)
                            .ok_or_else(|| anyhow::anyhow!("bad range: {}..{}", s, e))?;
                        let channel_q = channel_q_from_fastq_ascii(q_slice)?;
                        channel_q >= q_threshold as f32
                    } else {
                        // counter!("q_lower_than_threshold").increment(1);
                        // histogram!("q_threshold_value").record(q_threshold as f64);
                        false
                    };
                    if q_ok && len_ok {
                        (make_name(&rec.id, bc.name.as_ref(), s, e), s, e, Some(bc))
                    } else {
                        counter!("filtered_by_length_or_q").increment(1);

                        (
                            make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                            0,
                            rec.sequence.len(),
                            None,
                        )
                    }
                }
                None => (
                    make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                    0,
                    rec.sequence.len(),
                    None,
                ),
            };
            let fastq_rec = fastq::Record::with_attrs(
                &name,                     // &str
                None,                      // description/comment
                &rec.sequence[start..end], // &[u8]
                rec.quality
                    .as_deref() // Option<&[u8]>
                    .expect("missing quality")
                    .get(start..end) // Option<&[u8]> 取子片段
                    .ok_or_else(|| anyhow::anyhow!("bad range: {}..{}", start, end))?, // quality: Option<&[u8]>
                                                                                       // .as_ref(),
            );
            return Ok((RecordType::Fastq(fastq_rec), bc_opt));
        }
        // ---------- BAM ----------
        "bam" => {
            // let channel_num = rec.id.split('_').nth(1).unwrap_or("0").parse().unwrap_or(0);
            let channel_num = rec.ch.unwrap_or(0);
            let mut read_q = 7.0;
            let (name, start, end, bc_opt) = match barcode_position {
                Some(bc) => {
                    let (s, e) = bc.inner_position;
                    let len_ok = e.saturating_sub(s) >= min_subread_len;

                    let q_ok = if len_ok {
                        let q_slice = rec
                            .quality
                            .as_deref()
                            .ok_or_else(|| anyhow::anyhow!("missing quality"))?
                            .get(s..e)
                            .ok_or_else(|| anyhow::anyhow!("bad range: {}..{}", s, e))?;
                        read_q = channel_q_from_base_q(q_slice)?;
                        read_q >= q_threshold as f32
                    } else {
                        // counter!("q_lower_than_threshold").increment(1);
                        // histogram!("q_threshold_value").record(q_threshold as f64);
                        false
                    };
                    if q_ok && len_ok {
                        (make_name(&rec.id, bc.name.as_ref(), s, e), s, e, Some(bc))
                    } else {
                        counter!("filtered_by_length_or_q").increment(1);

                        (
                            make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                            0,
                            rec.sequence.len(),
                            None,
                        )
                    }
                }
                None => (
                    make_name(&rec.id, "uncertain", 0, rec.sequence.len()),
                    0,
                    rec.sequence.len(),
                    None,
                ),
            };
            // 头部的数据够长才行

            let mut bam_rec = bam::Record::new();

            bam_rec.set(
                name.as_bytes(),
                None,
                &rec.sequence[start..end],
                rec.quality
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("missing quality"))?
                    .get(start..end)
                    .ok_or_else(|| anyhow::anyhow!("bad range: {}..{}", start, end))?,
            ); // 第 2 个参数是 CIGAR 向量
            bam_rec.set_flags(0x4);

            bam_rec.push_aux(b"ch", Aux::I32(channel_num)).unwrap();
            bam_rec.push_aux(b"RG", Aux::String("0425")).unwrap();
            bam_rec
                .push_aux(b"rq", Aux::Float(rec.rq.unwrap_or(0.777)))
                .unwrap();
            bam_rec.push_aux(b"cq", Aux::Float(read_q as f32)).unwrap();
            bam_rec
                .push_aux(b"np", Aux::I32(rec.np.unwrap_or(0)))
                .unwrap();

            let start_end_position = vec![start as u32, end as u32];

            bam_rec
                .push_aux(b"be", Aux::ArrayU32((&start_end_position[..]).into()))
                .unwrap();
            return Ok((RecordType::BAM(bam_rec), bc_opt));
        }

        // ---------- 其他格式 ----------
        _ => anyhow::bail!("不支持的输出格式: {output_format}"),
    }
}

pub fn channel_q_from_fastq_ascii(q_ascii: &[u8]) -> anyhow::Result<f32> {
    anyhow::ensure!(!q_ascii.is_empty(), "empty qual");
    let mut sum_p = 0.0f64;

    for &c in q_ascii {
        anyhow::ensure!(c >= 33, "invalid FASTQ qual byte: {c}");
        let q = (c - 33) as f64;
        let p = 10f64.powf(-q / 10.0);
        sum_p += p;
    }

    let p_mean = (sum_p / q_ascii.len() as f64).max(1e-300);
    Ok((-10.0 * p_mean.log10()) as f32)
}

pub fn channel_q_from_base_q(qs: &[u8]) -> Result<f32> {
    if qs.is_empty() {
        return Err(anyhow!("empty Q array"));
    }

    // 平均错误概率
    let mut sum_p = 0.0f32;
    for &q in qs {
        // p = 10^(-q/10)
        let p = 10f32.powf(-(q as f32) / 10.0);
        sum_p += p;
    }
    let p_mean = sum_p / (qs.len() as f32);

    // 防止 log10(0)
    let p_mean = p_mean.max(1e-300);
    Ok(-10.0 * p_mean.log10())
}
