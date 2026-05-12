use crate::barcode_demuxer_v1::barcode_demuxer_v1;
use crate::barcode_mysers::{Direction, get_myers_from_barcodes};

use crate::bam_record_extention::ReadRecord;
// use crate::core_context::{Annotated, BarcodeCandidate, PrimerMeta};
// use crate::find_pattern::merge_non_overlapping_no_copy;
use crate::get_demuxed_reads::{RecordType, prepare_record_to_writer};
use crossbeam::channel::{Receiver, Sender};

// use metrics::{self, counter, histogram};

use crate::core_context::BarcodePair;

// use tracing::{debug, info, level_filters};

pub fn demux_reads_by_barcode(
    patterns: &[ReadRecord],
    max_distance: u8,
    receiver: Receiver<ReadRecord>,
    sender: Sender<(RecordType, Option<BarcodePair>)>,
    // ouput_folder: &str,
    min_subread_len: usize,
    min_q: usize,
    output_format: &str,
    search_bound: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut barcode_f_myers = get_myers_from_barcodes(patterns, Direction::Forward);
    let mut barcode_r_myers = get_myers_from_barcodes(patterns, Direction::Reverse);
    // println!("{:?}", barcode_f_myers.keys());
    for current_record in receiver {
        // let mut current_primer_distance: u8 = max_primer_distance;
        // // TODO 在这里限制reads的长度
        //     counter!("filtered_reads_too_short").increment(1 as u64);

        // print!("max primer tolerance: {}", max_primer_distance);
        let barcode_info = barcode_demuxer_v1(
            min_subread_len,
            max_distance,
            &mut barcode_f_myers,
            &mut barcode_r_myers,
            &current_record,
            search_bound,
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
