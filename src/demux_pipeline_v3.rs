use crate::barcode_mysers::{Direction, get_myers_from_barcodes};
use crate::cli::Cli;
use crate::core_context::BarcodePair;
use crate::find_pattern::get_pattern_keys;
use crate::get_demuxed_reads::RecordType;
use std::collections::HashSet;
use tracing::info;
// use crate::multiple_barcode_demuxer_v1::multiple_barcode_demuxer_v1;
use crate::bam_record_extention::ReadRecord;
use crate::core_context::BarcodeCandidate;
use crate::io_utils::normalize_barcode_name;
use crate::writer_worker::write_barcode_results;
use crossbeam::channel::{Receiver, Sender};
use std::error::Error;
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

pub fn demux_pipeline_v3(cli: &Cli) -> Result<(), Box<dyn Error>> {
    ensure_output_dir(&cli.output_folder)?;
    let mut patterns = read_sequences(&cli.barcode)?;
    let output_names = get_pattern_keys(&patterns);

    add_r_counterparts_to_records(&mut patterns);

    type WriteMsg = (RecordType, Option<BarcodePair>);
    let queue_len = 10000;
    std::thread::scope(|scope| {
        let (seq_sender, seq_receiver) = bounded(queue_len);
        // let (primer_sender, primer_receiver) = bounded(100000);
        let (barcode_sender, barcode_receiver): (Sender<WriteMsg>, Receiver<WriteMsg>) =
            bounded(queue_len);
        // let
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
                        cli.keep_barcode,
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
    keep_barcode: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut barcode_f_myers = get_myers_from_barcodes(patterns, Direction::Forward);
    let mut barcode_r_myers = get_myers_from_barcodes(patterns, Direction::Reverse);
    // println!("{:?}", barcode_f_myers.keys());
    for current_record in receiver {
        // let mut current_primer_distance: u8 = max_primer_distance;
        // // TODO 在这里限制reads的长度
        //     counter!("filtered_reads_too_short").increment(1 as u64);

        // print!("max primer tolerance: {}", max_primer_distance);
        let barcode_info = multiple_barcode_demuxer_v3(
            min_subread_len,
            max_distance,
            &mut barcode_f_myers,
            &mut barcode_r_myers,
            &current_record,
            search_bound,
            single_end_filter,
        )?;
        let demuxed_reads = prepare_record_to_writer(
            &current_record,
            barcode_info,
            min_subread_len,
            // ouput_folder,
            output_format,
            min_q,
            keep_barcode,
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

pub fn multiple_barcode_demuxer_v3(
    min_read_length: usize,
    max_edit_distance: u8,
    barcode_f_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    barcode_r_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    target: &ReadRecord,
    search_bound: usize,
    single_end_filter: bool,
) -> Result<Option<BarcodePair>> {
    use std::collections::HashMap;

    let read_len = target.sequence.len();
    if read_len <= 2 * min_read_length {
        counter!("len_fail").increment(1);
        return Ok(None);
    }
    let updated_search_bound = search_bound.min(read_len / 2);
    let leading_keys: Vec<_> = barcode_f_pattern.keys().cloned().collect();
    let mut leading_candidates = Vec::new();
    for name in leading_keys.iter() {
        let myers = barcode_f_pattern.get_mut(name).unwrap();
        if let Some(candidate) = barcode_alignment(
            name.clone(),
            myers,
            &target.sequence[0..updated_search_bound],
            max_edit_distance,
        )? {
            leading_candidates.push(candidate);
        }
    }

    let trailing_keys: Vec<_> = barcode_r_pattern.keys().cloned().collect();
    let mut trailing_candidates = Vec::new();
    for name in trailing_keys.iter() {
        let myers = barcode_r_pattern.get_mut(name).unwrap();
        if let Some(candidate) = barcode_alignment(
            name.clone(),
            myers,
            &target.sequence[read_len.saturating_sub(updated_search_bound)..read_len],
            max_edit_distance,
        )? {
            trailing_candidates.push(candidate);
        }
    }

    // 只有双端都有候选时才尝试配对
    if !trailing_candidates.is_empty() && !leading_candidates.is_empty() {
        // 每个 name 在 leading 中只保留最小 distance
        let mut best_leading: HashMap<Arc<str>, BarcodeCandidate> = HashMap::new();
        for c in leading_candidates.iter() {
            match best_leading.get_mut(&c.name) {
                Some(old) => {
                    if c.distance < old.distance {
                        *old = c.clone();
                    }
                }
                None => {
                    best_leading.insert(c.name.clone(), c.clone());
                }
            }
        }

        // 每个 name 在 trailing 中只保留最小 distance
        let mut best_trailing: HashMap<Arc<str>, BarcodeCandidate> = HashMap::new();
        for c in trailing_candidates.iter() {
            match best_trailing.get_mut(&c.name) {
                Some(old) => {
                    if c.distance < old.distance {
                        *old = c.clone();
                    }
                }
                None => {
                    best_trailing.insert(c.name.clone(), c.clone());
                }
            }
        }

        // 只比较两边共有的 name
        let mut shared_pairs: Vec<(BarcodeCandidate, BarcodeCandidate, u16)> = Vec::new();
        for (name, lead) in best_leading.iter() {
            if let Some(trail) = best_trailing.get(name) {
                let total_distance = lead.distance as u16 + trail.distance as u16;
                shared_pairs.push((lead.clone(), trail.clone(), total_distance));
            }
        }

        // 没有共有 key
        if shared_pairs.is_empty() {
            let best_lead = best_leading.values().min_by_key(|c| c.distance);

            let best_trail = best_trailing.values().min_by_key(|c| c.distance);

            match (best_lead, best_trail) {
                (Some(lead), Some(trail)) => {
                    info!(
                        "{} has different barcode on each side: leading {} distance {}, trailing {} distance {}",
                        target.id,
                        lead.name.as_ref(),
                        lead.distance,
                        trail.name.as_ref(),
                        trail.distance,
                    );
                }

                _ => {
                    info!(
                        "{} has no shared barcode between leading and trailing candidates",
                        target.id,
                    );
                }
            }
            counter!("len_ok_pair_fail_no_shared_barcode").increment(1);
            return Ok(None);
        }

        shared_pairs.sort_by_key(|x| x.2);

        // 最优并列，放弃
        if shared_pairs.len() >= 2 && shared_pairs[0].2 == shared_pairs[1].2 {
            counter!("len_ok_pair_fail_score_tie").increment(1);
            return Ok(None);
        }

        let (lead, trail, _) = &shared_pairs[0];

        counter!("len_ok_pair_ok").increment(1);
        let result = BarcodePair {
            name: normalize_barcode_name(lead.name.as_ref()),
            distance: (lead.distance, trail.distance),
            inner_position: (lead.end, read_len - updated_search_bound + trail.start),
            outter_position: (lead.start, read_len - updated_search_bound + trail.end),
        };
        return Ok(Some(result));
    } else if trailing_candidates.len() >= 1 && leading_candidates.len() == 0 {
        // counter!("len_ok_pair_ok").increment(1);

        if !single_end_filter {
            trailing_candidates.sort_by_key(|n| n.distance);

            if trailing_candidates.len() >= 2
                && trailing_candidates[0].distance == trailing_candidates[1].distance
            {
                counter!("len_ok_pair_fail_trailing_only_score_tie").increment(1);
                return Ok(None);
            }
            let result = BarcodePair {
                name: trailing_candidates[0].name.clone(),
                distance: (100, trailing_candidates[0].distance),
                inner_position: (
                    0,
                    read_len - updated_search_bound + trailing_candidates[0].start,
                ),
                outter_position: (
                    0,
                    read_len - updated_search_bound + trailing_candidates[0].end,
                ),
            };
            counter!("len_ok_single_end_trailing_ok").increment(1);

            return Ok(Some(result));
        } else {
            counter!("len_ok_both_ends_fail_trailing_only").increment(1);

            return Ok(None);
        }

        // return Ok(None);
    } else if leading_candidates.len() >= 1 && trailing_candidates.len() == 0 {
        // counter!("len_ok_pair_ok").increment(1);
        if !single_end_filter {
            leading_candidates.sort_by_key(|n| n.distance);

            if leading_candidates.len() >= 2
                && leading_candidates[0].distance == leading_candidates[1].distance
            {
                counter!("len_ok_pair_fail_leading_only_score_tie").increment(1);
                return Ok(None);
            }
            let result = BarcodePair {
                name: leading_candidates[0].name.clone(),
                distance: (leading_candidates[0].distance, 100),
                inner_position: (leading_candidates[0].end, read_len),
                outter_position: (leading_candidates[0].start, read_len),
            };
            counter!("len_ok_single_end_leading_ok").increment(1);

            return Ok(Some(result));
        } else {
            counter!("len_ok_both_ends_fail_leading_only").increment(1);

            return Ok(None);
        }
    } else {
        counter!("len_ok_score_fail_no_barcode_found").increment(1);
        return Ok(None);
    }

    // 只要不是双端成功配对，一律返回 None
    // counter!("len_ok_score_fail").increment(1);
    // Ok(None)
}

fn add_r_counterparts_to_records(patterns: &mut Vec<ReadRecord>) {
    /*
     * 规则：
     *
     * xxx_F -> 补 xxx_R_R
     * xxx_R -> 补 xxx_R_F
     *
     * 注意：
     * 这里只复制原 record，并修改 id。
     * sequence / quality / tags 等内容保持不变。
     */
    let mut existing_names: HashSet<String> = patterns
        .iter()
        .map(|record| record.id.to_string())
        .collect();

    let mut extra_records: Vec<ReadRecord> = Vec::new();

    for record in patterns.iter() {
        let name = record.id.as_ref();

        if let Some(prefix) = name.strip_suffix("_F") {
            let new_name = format!("{}_R_R", prefix);

            if existing_names.insert(new_name.clone()) {
                let mut new_record = record.clone();
                new_record.id = Arc::<str>::from(new_name);
                extra_records.push(new_record);
            }
        } else if let Some(prefix) = name.strip_suffix("_R") {
            let new_name = format!("{}_R_F", prefix);

            if existing_names.insert(new_name.clone()) {
                let mut new_record = record.clone();
                new_record.id = Arc::<str>::from(new_name);
                extra_records.push(new_record);
            }
        }
    }

    patterns.extend(extra_records);
}
