use metrics::{Unit, describe_counter};

/// Register descriptions for all demux metrics.
///
/// These descriptions will be emitted as `# HELP` lines in the
/// Prometheus-format metrics output.
///
/// This function should be called once after the metrics recorder
/// is initialized and before processing reads.
pub fn describe_metrics() {
    describe_counter!(
        "len_fail",
        Unit::Count,
        "Number of reads that failed the minimum length requirement. These reads usually do not enter the barcode matching step."
    );

    describe_counter!(
        "len_ok_score_fail_no_barcode_found",
        Unit::Count,
        "Number of reads that passed the length filter but no valid barcode was found within the leading or trailing search regions."
    );

    describe_counter!(
        "len_ok_single_end_leading_ok",
        Unit::Count,
        "Number of reads that passed the length filter and only matched a barcode on the leading end, while no barcode was found on the trailing end."
    );

    describe_counter!(
        "len_ok_single_end_trailing_ok",
        Unit::Count,
        "Number of reads that passed the length filter and only matched a barcode on the trailing end, while no barcode was found on the leading end."
    );

    describe_counter!(
        "len_ok_pair_ok",
        Unit::Count,
        "Number of reads that passed the length filter, matched barcodes on both leading and trailing ends, and formed a valid barcode pair."
    );

    describe_counter!(
        "len_ok_pair_ok_q_or_length_failed",
        Unit::Count,
        "Number of reads with a valid paired barcode match, but the extracted subread failed the quality or length filter."
    );

    describe_counter!(
        "len_ok_pair_fail_no_shared_barcode",
        Unit::Count,
        "Number of reads that passed the length filter and had barcode candidates on both ends, but the leading and trailing barcodes could not form a valid shared barcode pair."
    );

    describe_counter!(
        "len_ok_pair_fail_score_tie",
        Unit::Count,
        "Number of reads that passed the length filter but had a tied best score during paired barcode matching, making the barcode pair ambiguous."
    );

    describe_counter!(
        "len_ok_pair_fail_leading_only_score_tie",
        Unit::Count,
        "Number of reads that passed the length filter but had a tied best score on the leading-end barcode match, making the leading barcode ambiguous."
    );

    describe_counter!(
        "len_ok_pair_fail_trailing_only_score_tie",
        Unit::Count,
        "Number of reads that passed the length filter but had a tied best score on the trailing-end barcode match, making the trailing barcode ambiguous."
    );

    describe_counter!(
        "writer_received",
        Unit::Count,
        "Number of reads received by the writer module. This should usually equal the sum of all writer_barcode_count values."
    );

    describe_counter!(
        "writer_barcode_count",
        Unit::Count,
        "Number of reads written by the writer module, grouped by barcode name. The barcode=\"uncertain\" label indicates reads that could not be uniquely assigned or did not pass demux classification."
    );
}
