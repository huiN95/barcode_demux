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
// use crate::core_context::{Annotated, BarcodeCandidate, PrimerMeta};
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

pub fn demux_pipeline_v2(cli: &Cli) -> Result<(), Box<dyn Error>> {
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
        scope.spawn(move || {
            if let Err(e) = read_sequences_to_queue(&cli.input_file, seq_sender) {
                eprintln!("[Producer] 错误: {}", e);
            }
            println!("[Producer] 已完成发送。");
        });

        // b) 多个 demux 线程
        // ---------------------------
        // let worker_count = 20;
        // let worker_count = thread::available_parallelism()
        //     .map(|n| n.get().saturating_sub(cli.reservesed_threads.into())) // 减去 4，不小于 0
        //     .unwrap_or(1);
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

pub fn demux_reads_by_multiple_barcode(
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
) -> Result<(), Box<dyn std::error::Error>> {
    let mut barcode_f_myers = get_myers_from_barcodes(patterns, Direction::Forward);
    let mut barcode_r_myers = get_myers_from_barcodes(patterns, Direction::Reverse);
    // println!("{:?}", barcode_f_myers.keys());
    for current_record in receiver {
        // let mut current_primer_distance: u8 = max_primer_distance;
        // // TODO 在这里限制reads的长度
        //     counter!("filtered_reads_too_short").increment(1 as u64);

        // print!("max primer tolerance: {}", max_primer_distance);
        let barcode_info = multiple_barcode_demuxer_v1(
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

pub fn multiple_barcode_demuxer_v1(
    min_read_length: usize,
    max_edit_distance: u8,

    barcode_f_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    barcode_r_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    target: &ReadRecord,
    search_bound: usize,
    single_end_filter: bool,
) -> Result<Option<BarcodePair>> {
    // let mut result: Vec<BarcodePair> = Vec::<BarcodePair>::new();
    // 如果序列长度小于等于2倍的intial_primer_check_len, 直接返回空结果
    let read_len = target.sequence.len();
    if read_len <= 2 * min_read_length {
        counter!("len_fail").increment(1 as u64);
        return Ok(None);
    }

    // let mut all_primer_pos = IndexMap::<Arc<str>, Vec<BarcodeCandidate>>::new();
    // let mut leading_primer_pos = Vec::<BarcodeCandidate>::new();

    // loop all the patterns, to see which one is the best match
    let leading_keys: Vec<_> = barcode_f_pattern.keys().cloned().collect();
    // leading_keys.sort();

    // println!("primer demux: {:?} patterns found, ", keys,);
    // for (name, myers) in inside_patterns_myers.iter_mut() {
    let mut leading_candidates = Vec::new();
    for name in leading_keys.iter() {
        // println!("primer name {name}");
        let myers = barcode_f_pattern.get_mut(name).unwrap();

        if let Some(candidate) = barcode_alignment(
            name.clone(),
            myers,
            &target.sequence[0..search_bound],
            max_edit_distance,
        )? {
            leading_candidates.push(candidate);
        }
    }
    // let leading_primer_pos = merge_non_overlapping_no_copy(&mut candidates);

    let trailing_keys: Vec<_> = barcode_r_pattern.keys().cloned().collect();
    // trailing_keys.sort();
    // let mut trailing_primer_pos = Vec::<BarcodeCandidate>::new();
    let mut trailing_candidates = Vec::new();

    for name in trailing_keys.iter() {
        let myers: &mut MayersPattern = barcode_r_pattern.get_mut(name).unwrap();

        if let Some(candidate) = barcode_alignment(
            name.clone(),
            myers,
            &target.sequence[target.sequence.len() - search_bound..target.sequence.len()],
            max_edit_distance,
        )? {
            trailing_candidates.push(candidate);
        }
        // 理论上都不应该有重叠
    }

    if trailing_candidates.len() >= 1 && leading_candidates.len() >= 1 {
        trailing_candidates.sort_by_key(|n| n.distance);
        leading_candidates.sort_by_key(|n| n.distance);
        if trailing_candidates[0].name == leading_candidates[0].name {
            counter!("len_ok_pair_ok").increment(1);
            let result = BarcodePair {
                name: trailing_candidates[0].name.clone(),
                distance: (
                    leading_candidates[0].distance,
                    trailing_candidates[0].distance,
                ),
                inner_position: (
                    leading_candidates[0].end,
                    read_len - search_bound + trailing_candidates[0].start,
                ),
                outter_position: (
                    leading_candidates[0].start,
                    read_len - search_bound + trailing_candidates[0].end,
                ),
            };
            return Ok(Some(result));
        } else {
            counter!("len_ok_pair_fail").increment(1);
            return Ok(None);
        }
    } else if trailing_candidates.len() >= 1 && leading_candidates.len() == 0 && !single_end_filter
    {
        trailing_candidates.sort_by_key(|n| n.distance);
        counter!("len_ok_pair_ok").increment(1);

        counter!("len_ok_pair_ok_trailing_only").increment(1);
        let result = BarcodePair {
            name: trailing_candidates[0].name.clone(),
            distance: (0, trailing_candidates[0].distance),
            inner_position: (0, read_len - search_bound + trailing_candidates[0].start),
            outter_position: (0, read_len - search_bound + trailing_candidates[0].end),
        };
        return Ok(Some(result));

        // return Ok(None);
    } else if leading_candidates.len() >= 1 && trailing_candidates.len() == 0 && !single_end_filter
    {
        leading_candidates.sort_by_key(|n| n.distance);
        counter!("len_ok_pair_ok").increment(1);
        counter!("len_ok_pair_ok_leading_only").increment(1);
        let result = BarcodePair {
            name: leading_candidates[0].name.clone(),
            distance: (leading_candidates[0].distance, 0),
            inner_position: (leading_candidates[0].end, read_len),
            outter_position: (leading_candidates[0].start, read_len),
        };
        return Ok(Some(result));
    } else {
        counter!("len_ok_score_fail").increment(1);
        return Ok(None);
    }

    // if trailing_candidates.len() == 1 && leading_candidates.len() == 1 {
    //     if trailing_candidates[0].name == leading_candidates[0].name {
    //         counter!("reads with barcodes").increment(1);
    //         let result = BarcodePair {
    //             name: trailing_candidates[0].name.clone(),
    //             distance: (
    //                 leading_candidates[0].distance,
    //                 trailing_candidates[0].distance,
    //             ),
    //             inner_position: (
    //                 leading_candidates[0].end,
    //                 read_len - search_bound + trailing_candidates[0].start,
    //             ),
    //             outter_position: (
    //                 leading_candidates[0].start,
    //                 read_len - search_bound + trailing_candidates[0].end,
    //             ),
    //         };
    //         return Ok(Some(result));
    //     } else {
    //         counter!("barcode_not_paired").increment(1);
    //         return Ok(None);
    //     }
    // } else {
    //     counter!("failed_to_demux_barcode").increment(1);
    //     return Ok(None);
    // }
}
