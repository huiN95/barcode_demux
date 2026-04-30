use crate::cli::Cli;

use crate::barcode_mysers::{Direction, get_myers_from_barcodes};
use crate::core_context::BarcodePair;
use crate::find_pattern::get_pattern_keys;
use crate::get_demuxed_reads::RecordType;
// use crate::multiple_barcode_demuxer_v1::multiple_barcode_demuxer_v1;
use crate::writer_worker::write_barcode_results;
use crossbeam::channel::{Receiver, Sender};
use std::error::Error;

use crate::bam_record_extention::ReadRecord;
use crate::core_context::BarcodeCandidate;
// use crate::find_pattern::merge_non_overlapping_no_copy;
use crate::get_demuxed_reads::prepare_record_to_writer;

// use metrics::{self, counter, histogram};

use std::thread;
// use crate::core_context::BarcodeCandidate;
// use crate::find_pattern::merge_non_overlapping_no_copy;
use std::sync::Arc;

use crate::barcode_mysers::{MayersPattern, barcode_alignment};

use bio::alignment::distance;
use metrics::{self, counter, histogram};
// use rayon::result;

use crate::io_utils::{ensure_output_dir, read_sequences};
use crate::reader_worker::read_sequences_to_queue;
use anyhow::Result;
use crossbeam::channel::bounded;
use std::collections::HashMap;

pub fn demux_pipeline_v4(cli: &Cli) -> Result<(), Box<dyn Error>> {
    ensure_output_dir(&cli.output_folder)?;
    let patterns = read_sequences(&cli.barcode)?;

    let output_names = get_pattern_keys(&patterns);

    type WriteMsg = (RecordType, Option<BarcodePair>);
    let queue_len = 10000;
    std::thread::scope(|scope| {
        let (seq_sender, seq_receiver) = bounded(queue_len);
        // let (primer_sender, primer_receiver) = bounded(100000);
        let (barcode_sender, barcode_receiver): (Sender<WriteMsg>, Receiver<WriteMsg>) =
            bounded(queue_len);
        // let
        // let min_pair_len = cli.min_pair_len;
        scope.spawn(move || {
            if let Err(e) = read_sequences_to_queue(&cli.input_file, seq_sender) {
                eprintln!("[Producer] 错误: {}", e);
            }
            println!("[Producer] 已完成发送。");
        });

        // b) 多个 demux 线程
        // ---------------------------
        // let worker_count = 20;
        let worker_count: usize = match cli.threads {
            0 => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
                .saturating_sub(2)
                .max(1),
            n => n as usize,
        };
        // .min(1)
        // .max(1); // .unwrap_or(1); // 获取失败时默认 1
        //                // .max(1)
        //                // .min(1);
        // .min(5); // 确保最终值不小于 10

        println!("[Producer] 启动 {} 个 demux 线程...\n", worker_count);
        let patterns = &patterns;
        let min_pari_len = cli.min_pair_len.unwrap();
        let min_pair_score = cli.min_pair_score.unwrap();
        for _ in 0..worker_count {
            // 如果 patterns 很大，需要 Arc::clone；若直接 patterns.clone() 也可以

            let seq_receiver_clone = seq_receiver.clone();
            let barcode_sender_clone = barcode_sender.clone();
            let max_distance = cli.max_distance;

            scope.spawn(
                move || {
                    // println!("Demuxing with primers...");
                    if let Err(e) = demux_reads_by_multiple_barcode(
                        patterns,
                        max_distance,
                        seq_receiver_clone,
                        barcode_sender_clone,
                        // &cli.output_folder,
                        cli.min_subread_len,
                        cli.min_q,
                        &cli.output_format,
                        cli.search_bound,
                        cli.single_end_filter,
                        min_pari_len,
                        min_pair_score,
                    ) {
                        eprintln!("[Demux] error: {e}");
                    }

                    // println!("[Demux] finished");
                }, // 同理，如果需要让 writer 知道没有更多数据，可在所有 demux 结束前
                   // 最后一个线程里 drop(primer_sender_clone)。
            );
        }
        drop(seq_receiver);
        drop(barcode_sender);
        // print!(
        //     "save barcode flag {:?} after multiprcoessing",
        //     cli.keep_primer
        // );

        write_barcode_results(
            &cli.input_file,
            &cli.output_format,
            &cli.output_folder,
            barcode_receiver,
            output_names,
        )
        .unwrap();
    });

    Ok(())
}

fn demux_reads_by_multiple_barcode(
    patterns: &[ReadRecord],
    max_distance: u8,
    receiver: Receiver<ReadRecord>,
    sender: Sender<(RecordType, Option<BarcodePair>)>,
    // ouput_folder: &str,
    min_subread_len: usize,
    min_q: usize,
    output_format: &str,
    search_bound: usize,
    single_end_filter: bool,
    min_pair_len: usize,
    min_pari_score: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let cfg = LocalAlignmentConfig {
        min_identity: min_pari_score,
        min_alignment_length: min_pair_len,
        match_score: 1.0,
        mismatch_penalty: 1.0,
        gap_penalty: 1.0,
        allow_ambiguous_match: true,
    };
    for current_record in receiver {
        // let mut current_primer_distance: u8 = max_primer_distance;
        // // TODO 在这里限制reads的长度
        //     counter!("filtered_reads_too_short").increment(1 as u64);

        // print!("max primer tolerance: {}", max_primer_distance);
        let barcode_info = multiple_barcode_demuxer_v4(
            min_subread_len,
            patterns,
            &current_record,
            search_bound,
            single_end_filter,
            cfg,
        )?;
        let demuxed_reads = prepare_record_to_writer(
            &current_record,
            barcode_info,
            min_subread_len,
            // ouput_folder,
            output_format,
            min_q,
            false,
        )?;
        // println!(demuxed_reads);
        if let Err(e) = sender.send(demuxed_reads) {
            eprintln!("Failed to send demux: {}", e);
            // Optionally return early or break here
        }
        // }
    }
    Ok(())
}

pub fn multiple_barcode_demuxer_v4(
    min_read_length: usize,
    // max_edit_distance: u8,
    patterns: &[ReadRecord],
    target: &ReadRecord,
    search_bound: usize,
    _single_end_filter: bool,
    config: LocalAlignmentConfig,
) -> Result<Option<BarcodePair>> {
    let read_len = target.sequence.len();
    if read_len <= 2 * min_read_length {
        counter!("len_fail").increment(1);
        return Ok(None);
    }
    let updated_bound = read_len.saturating_sub(search_bound);
    // let mut query_candidates = Vec::new();
    for query in patterns.iter() {
        if let Some(pair_resutl) = longest_similar_local_alignment(
            &query.sequence,
            &target.sequence[updated_bound..read_len],
            config,
        ) {
            let result = BarcodePair {
                name: query.id.clone(),
                distance: (0, pair_resutl.edit_like_errors as i32),
                inner_position: (0, updated_bound + pair_resutl.b_range.0),
                outter_position: (0, updated_bound + pair_resutl.b_range.1),
            };
            return Ok(Some(result));
        }
    }

    // counter!("len_ok_score_fail").increment(1);
    Ok(None)
}

#[derive(Debug, Clone)]
pub struct LocalAlignmentResult {
    pub a_range: (usize, usize), // 左闭右开
    pub b_range: (usize, usize), // 左闭右开

    pub a_substring: String,
    pub b_substring: String,

    pub aligned_a: String,
    pub aligned_b: String,

    pub alignment_length: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub gaps_in_a: usize,
    pub gaps_in_b: usize,
    pub edit_like_errors: usize,

    pub identity: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct LocalAlignmentConfig {
    pub min_identity: f64,
    pub min_alignment_length: usize,

    pub match_score: f64,
    pub mismatch_penalty: f64,
    pub gap_penalty: f64,

    pub allow_ambiguous_match: bool,
}

impl Default for LocalAlignmentConfig {
    fn default() -> Self {
        Self {
            min_identity: 0.90,
            min_alignment_length: 1,
            match_score: 1.0,
            mismatch_penalty: 1.0,
            gap_penalty: 1.0,
            allow_ambiguous_match: true,
        }
    }
}

const TRACE_STOP: u8 = 0;
const TRACE_DIAG: u8 = 1;
const TRACE_UP: u8 = 2; // B 中是 gap
const TRACE_LEFT: u8 = 3; // A 中是 gap

const EPS: f64 = 1e-12;

#[inline]
fn idx(i: usize, j: usize, width: usize) -> usize {
    i * width + j
}

#[inline]
fn iupac_mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,

        b'R' => 0b0101, // A/G
        b'Y' => 0b1010, // C/T
        b'S' => 0b0110, // G/C
        b'W' => 0b1001, // A/T
        b'K' => 0b1100, // G/T
        b'M' => 0b0011, // A/C

        b'B' => 0b1110, // C/G/T
        b'D' => 0b1101, // A/G/T
        b'H' => 0b1011, // A/C/T
        b'V' => 0b0111, // A/C/G

        b'N' => 0b1111, // any
        _ => 0,
    }
}

#[inline]
fn chars_match(a: u8, b: u8, allow_ambiguous_match: bool) -> bool {
    if allow_ambiguous_match {
        let ma = iupac_mask(a);
        let mb = iupac_mask(b);
        if ma != 0 && mb != 0 {
            return (ma & mb) != 0;
        }
    }
    a.eq_ignore_ascii_case(&b)
}

#[inline]
fn better_cell(
    cand_score: f64,
    cand_len: usize,
    cand_matches: usize,
    best_score: f64,
    best_len: usize,
    best_matches: usize,
) -> bool {
    if cand_score > best_score + EPS {
        return true;
    }
    if (cand_score - best_score).abs() <= EPS {
        if cand_len > best_len {
            return true;
        }
        if cand_len == best_len && cand_matches > best_matches {
            return true;
        }
    }
    false
}

#[inline]
fn better_global(
    cand_len: usize,
    cand_score: f64,
    cand_matches: usize,
    best_len: usize,
    best_score: f64,
    best_matches: usize,
) -> bool {
    if cand_len > best_len {
        return true;
    }
    if cand_len == best_len {
        if cand_score > best_score + EPS {
            return true;
        }
        if (cand_score - best_score).abs() <= EPS && cand_matches > best_matches {
            return true;
        }
    }
    false
}

/// 在 A 和 B 中寻找“最长的相似局部片段”。
///
/// - 允许 mismatch
/// - 允许 indel
/// - 支持 IUPAC 模糊碱基匹配
/// - 只返回 identity >= min_identity 的最佳结果
///
/// 注意：
/// 这是一个很实用的局部比对版本；它优先找“长而且分数高”的局部对齐，
/// 再按 identity 过滤。对于大多数 DNA/RNA 工程场景已经很好用。
pub fn longest_similar_local_alignment(
    a: &[u8],
    b: &[u8],
    cfg: LocalAlignmentConfig,
) -> Option<LocalAlignmentResult> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    if !(0.0 < cfg.min_identity && cfg.min_identity <= 1.0) {
        return None;
    }

    let n = a.len();
    let m = b.len();
    let width = m + 1;
    let size = (n + 1) * (m + 1);

    // Smith-Waterman DP
    let mut score = vec![0.0_f64; size];
    let mut trace = vec![TRACE_STOP; size];

    // 记录“当前路径”的对齐列数与 match 数，便于在线判断 identity
    let mut aln_len = vec![0_usize; size];
    let mut match_cnt = vec![0_usize; size];

    // 最佳合法终点
    let mut best_end_i = 0_usize;
    let mut best_end_j = 0_usize;
    let mut best_len = 0_usize;
    let mut best_score = 0.0_f64;
    let mut best_matches = 0_usize;

    for i in 1..=n {
        for j in 1..=m {
            let cur = idx(i, j, width);
            let diag = idx(i - 1, j - 1, width);
            let up = idx(i - 1, j, width);
            let left = idx(i, j - 1, width);

            let is_match = chars_match(a[i - 1], b[j - 1], cfg.allow_ambiguous_match);

            let diag_score = score[diag]
                + if is_match {
                    cfg.match_score
                } else {
                    -cfg.mismatch_penalty
                };
            let diag_len = aln_len[diag] + 1;
            let diag_matches = match_cnt[diag] + usize::from(is_match);

            let up_score = score[up] - cfg.gap_penalty;
            let up_len = aln_len[up] + 1;
            let up_matches = match_cnt[up];

            let left_score = score[left] - cfg.gap_penalty;
            let left_len = aln_len[left] + 1;
            let left_matches = match_cnt[left];

            let mut best_here_score = 0.0_f64;
            let mut best_here_len = 0_usize;
            let mut best_here_matches = 0_usize;
            let mut best_here_trace = TRACE_STOP;

            if diag_score > 0.0
                && better_cell(
                    diag_score,
                    diag_len,
                    diag_matches,
                    best_here_score,
                    best_here_len,
                    best_here_matches,
                )
            {
                best_here_score = diag_score;
                best_here_len = diag_len;
                best_here_matches = diag_matches;
                best_here_trace = TRACE_DIAG;
            }

            if up_score > 0.0
                && better_cell(
                    up_score,
                    up_len,
                    up_matches,
                    best_here_score,
                    best_here_len,
                    best_here_matches,
                )
            {
                best_here_score = up_score;
                best_here_len = up_len;
                best_here_matches = up_matches;
                best_here_trace = TRACE_UP;
            }

            if left_score > 0.0
                && better_cell(
                    left_score,
                    left_len,
                    left_matches,
                    best_here_score,
                    best_here_len,
                    best_here_matches,
                )
            {
                best_here_score = left_score;
                best_here_len = left_len;
                best_here_matches = left_matches;
                best_here_trace = TRACE_LEFT;
            }

            score[cur] = best_here_score;
            aln_len[cur] = best_here_len;
            match_cnt[cur] = best_here_matches;
            trace[cur] = best_here_trace;

            if best_here_trace != TRACE_STOP && best_here_len >= cfg.min_alignment_length {
                let identity = best_here_matches as f64 / best_here_len as f64;
                if identity + EPS >= cfg.min_identity
                    && better_global(
                        best_here_len,
                        best_here_score,
                        best_here_matches,
                        best_len,
                        best_score,
                        best_matches,
                    )
                {
                    best_end_i = i;
                    best_end_j = j;
                    best_len = best_here_len;
                    best_score = best_here_score;
                    best_matches = best_here_matches;
                }
            }
        }
    }

    if best_len == 0 {
        return None;
    }

    // traceback
    let mut i = best_end_i;
    let mut j = best_end_j;

    let mut aligned_a = Vec::<u8>::with_capacity(best_len);
    let mut aligned_b = Vec::<u8>::with_capacity(best_len);

    let mut matches = 0_usize;
    let mut mismatches = 0_usize;
    let mut gaps_in_a = 0_usize;
    let mut gaps_in_b = 0_usize;

    while i > 0 && j > 0 {
        let cur = idx(i, j, width);
        match trace[cur] {
            TRACE_STOP => break,
            TRACE_DIAG => {
                let ca = a[i - 1];
                let cb = b[j - 1];
                aligned_a.push(ca);
                aligned_b.push(cb);

                if chars_match(ca, cb, cfg.allow_ambiguous_match) {
                    matches += 1;
                } else {
                    mismatches += 1;
                }

                i -= 1;
                j -= 1;
            }
            TRACE_UP => {
                let ca = a[i - 1];
                aligned_a.push(ca);
                aligned_b.push(b'-');
                gaps_in_b += 1;
                i -= 1;
            }
            TRACE_LEFT => {
                let cb = b[j - 1];
                aligned_a.push(b'-');
                aligned_b.push(cb);
                gaps_in_a += 1;
                j -= 1;
            }
            _ => break,
        }
    }

    let start_i = i;
    let start_j = j;

    aligned_a.reverse();
    aligned_b.reverse();

    let aligned_a = String::from_utf8(aligned_a).ok()?;
    let aligned_b = String::from_utf8(aligned_b).ok()?;

    let alignment_length = aligned_a.len();
    if alignment_length == 0 {
        return None;
    }

    let identity = matches as f64 / alignment_length as f64;
    if identity + EPS < cfg.min_identity || alignment_length < cfg.min_alignment_length {
        return None;
    }

    Some(LocalAlignmentResult {
        a_range: (start_i, best_end_i),
        b_range: (start_j, best_end_j),

        a_substring: String::from_utf8_lossy(&a[start_i..best_end_i]).into_owned(),
        b_substring: String::from_utf8_lossy(&b[start_j..best_end_j]).into_owned(),

        aligned_a,
        aligned_b,

        alignment_length,
        matches,
        mismatches,
        gaps_in_a,
        gaps_in_b,
        edit_like_errors: mismatches + gaps_in_a + gaps_in_b,

        identity,
        score: best_score,
    })
}

pub fn longest_similar_local_alignment_str(
    a: &str,
    b: &str,
    cfg: LocalAlignmentConfig,
) -> Option<LocalAlignmentResult> {
    longest_similar_local_alignment(a.as_bytes(), b.as_bytes(), cfg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_local_alignment_with_indel() {
        let a = "TTTACGTTACGATTT";
        let b = "GGGACGTACGATGGG";

        let cfg = LocalAlignmentConfig {
            min_identity: 0.80,
            min_alignment_length: 6,
            match_score: 2.0,
            mismatch_penalty: 1.0,
            gap_penalty: 1.0,
            allow_ambiguous_match: true,
        };

        let res = longest_similar_local_alignment_str(a, b, cfg).unwrap();
        assert!(res.identity >= 0.80);
        assert!(res.alignment_length >= 6);
    }

    #[test]
    fn test_iupac_ambiguous_match() {
        let a = "ACGTRACGA";
        let b = "ACGTAACGA"; // R 可匹配 A/G

        let cfg = LocalAlignmentConfig {
            min_identity: 0.95,
            min_alignment_length: 5,
            match_score: 2.0,
            mismatch_penalty: 1.0,
            gap_penalty: 1.0,
            allow_ambiguous_match: true,
        };

        let res = longest_similar_local_alignment_str(a, b, cfg).unwrap();
        assert!(res.identity >= 0.95);
        assert!(res.matches >= 8);
    }

    #[test]
    fn test_disable_ambiguous_match() {
        let a = "ACGTRACGA";
        let b = "ACGTAACGA";

        let cfg = LocalAlignmentConfig {
            min_identity: 0.95,
            min_alignment_length: 5,
            match_score: 2.0,
            mismatch_penalty: 1.0,
            gap_penalty: 1.0,
            allow_ambiguous_match: false,
        };

        let res = longest_similar_local_alignment_str(a, b, cfg);
        if let Some(r) = res {
            assert!(r.identity < 1.0);
        }
    }
}
