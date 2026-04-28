use crate::cli::Cli;

use crate::core_context::BarcodePair;
use crate::demux_primer::demux_reads_by_barcode;
use crate::find_pattern::get_pattern_keys;
use crate::get_demuxed_reads::RecordType;
use crate::writer_worker::write_barcode_results;
use crossbeam::channel::{Receiver, Sender};
use std::error::Error;

use std::thread;
// use tracing::info;

// use tracing_appender;
// use tracing_subscriber::{fmt, EnvFilter};
use crate::io_utils::{ensure_output_dir, read_sequences};
use crate::reader_worker::read_sequences_to_queue;
use crossbeam::channel::bounded;

pub fn demux_pipeline_v1(cli: &Cli) -> Result<(), Box<dyn Error>> {
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
        // let worker_count = cli.threads;
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
                    if let Err(e) = demux_reads_by_barcode(
                        patterns,
                        max_distance,
                        seq_receiver_clone,
                        barcode_sender_clone,
                        // &cli.output_folder,
                        cli.min_subread_len,
                        cli.min_q,
                        &cli.output_format,
                        cli.search_bound,
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
