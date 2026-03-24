use crate::bam_record_extention::ReadRecord;
// use crate::core_context::BarcodeCandidate;
// use crate::find_pattern::merge_non_overlapping_no_copy;
use std::sync::Arc;

use crate::barcode_mysers::{barcode_alignment, MayersPattern};

use bio::alignment::distance;
use metrics::{self, counter, histogram};
// use rayon::result;

use crate::core_context::BarcodePair;
use anyhow::Result;
use std::collections::HashMap;
// use tracing::{debug, info, level_filters};

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

        counter!("len_ok_pair_pass_trailing").increment(1);
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
        counter!("len_ok_pair_pass_leading").increment(1);
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
