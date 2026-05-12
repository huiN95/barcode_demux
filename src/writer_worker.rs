use crate::core_context::BarcodePair;
use crate::get_demuxed_reads::RecordType;
use crate::io_utils::{make_writers, normalize_barcode_name};

use crate::pbar::{DEFAULT_INTERVAL, get_spin_pb};
use crossbeam::channel::Receiver;
use metrics::counter;
use tracing::{error, info};

pub fn write_barcode_results(
    input_file: &str,
    output_format: &str,
    output_folder: &str,
    barcode_demux_info: Receiver<(RecordType, Option<BarcodePair>)>,
    _barcode_names: Vec<String>,
) -> anyhow::Result<()> {
    println!("enter write process");

    let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);

    /*
     * 这里不再预初始化 metrics。
     *
     * 原因：
     * counter!(...).increment(0) 会注册 metric label，
     * 导致没有 reads 的 barcode 也出现在 metrics log 中。
     *
     * 现在只在真正写入成功后 increment(1)。
     */

    /*
     * 这里不再为所有 barcode_names 创建 writer。
     *
     * 原因：
     * make_writers(output_folder, barcode_names, ...)
     * 会提前为所有 barcode 创建文件，导致没有 reads 的 barcode
     * 也产生 0B 文件。
     *
     * 这里先创建空 writers，后面遇到某个 barcode 需要写入时，
     * 再为这个 barcode 懒创建 writer。
     */
    let mut writers = make_writers(
        output_folder,
        Vec::<String>::new(),
        input_file,
        output_format,
    )?;

    match output_format {
        "fa" | "fasta" | "fastq" | "fq" | "bam" => {
            while let Ok(demuxed_record) = barcode_demux_info.recv() {
                counter!("writer_received").increment(1_u64);

                let (demuxed_reads, barcode_pair) = demuxed_record;

                /*
                 * writer_key:
                 *   用于查找 / 创建 writer。
                 *   这里保持原始 barcode name。
                 *
                 * metric_barcode_name:
                 *   用于 metrics。
                 *   这里使用 normalize_barcode_name，避免 _F / _R / _R_F / _R_R
                 *   之类后缀导致统计被拆散。
                 */
                let (writer_key, metric_barcode_name) = match barcode_pair.as_ref() {
                    Some(barcode_pair) => {
                        let raw_name = barcode_pair.name.to_string();
                        let normalized_name = normalize_barcode_name(&raw_name).to_string();

                        (raw_name, normalized_name)
                    }
                    None => ("uncertain".to_string(), "uncertain".to_string()),
                };

                /*
                 * 懒创建 writer。
                 *
                 * 先看 writers 里是否已经有当前 barcode 的 writer。
                 * 如果没有，就只为当前 barcode 创建 writer。
                 */
                let writer_lookup_key: String = if writers.contains_key(writer_key.as_str()) {
                    writer_key.clone()
                } else if writers.contains_key(metric_barcode_name.as_str()) {
                    metric_barcode_name.clone()
                } else {
                    let new_writers = make_writers(
                        output_folder,
                        vec![writer_key.clone()],
                        input_file,
                        output_format,
                    )?;

                    let new_lookup_key = if new_writers.contains_key(writer_key.as_str()) {
                        writer_key.clone()
                    } else if new_writers.contains_key(metric_barcode_name.as_str()) {
                        metric_barcode_name.clone()
                    } else {
                        anyhow::bail!(
                            "new writer missing key: raw_key={}, normalized_key={}",
                            writer_key,
                            metric_barcode_name
                        );
                    };

                    writers.extend(new_writers);

                    new_lookup_key
                };

                let w = writers.get_mut(writer_lookup_key.as_str()).ok_or_else(|| {
                    anyhow::anyhow!(
                        "writers missing key after lazy creation: {}",
                        writer_lookup_key
                    )
                })?;

                match w.write_demuxed_record(&demuxed_reads) {
                    Ok(_) => {
                        /*
                         * 只在真正写入成功后统计。
                         *
                         * 因此：
                         *   1. 没有 reads 的 barcode 不会出现在 metric 里；
                         *   2. 写入失败的 read 不会被统计到 writer_barcode_count；
                         *   3. metric 中的名字是 normalize 之后的名字。
                         */
                        counter!(
                            "writer_barcode_count",
                            "barcode" => metric_barcode_name
                        )
                        .increment(1_u64);
                    }
                    Err(e) => {
                        counter!("writer_failed").increment(1_u64);

                        error!(
                            ?e,
                            barcode = metric_barcode_name.as_str(),
                            writer_key = writer_key.as_str(),
                            "demuxed read 写入失败"
                        );
                    }
                }

                pb.inc(1);
            }
        }
        _ => anyhow::bail!("不支持的写入格式: {output_format}"),
    }

    pb.finish_and_clear();

    info!(
        "barcode writing finished; writers were created lazily; metrics only include non-zero barcode counts"
    );

    Ok(())
}
// use crate::core_context::BarcodePair;
// use crate::get_demuxed_reads::RecordType;
// use crate::io_utils::{make_writers, normalize_barcode_name};

// use crate::pbar::{DEFAULT_INTERVAL, get_spin_pb};
// use crossbeam::channel::Receiver;
// use metrics::counter;
// use tracing::{error, info};

// pub fn write_barcode_results(
//     input_file: &str,
//     output_format: &str,
//     output_folder: &str,
//     barcode_demux_info: Receiver<(RecordType, Option<BarcodePair>)>,
//     barcode_names: Vec<String>,
// ) -> anyhow::Result<()> {
//     println!("enter write process");

//     let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);

//     /*
//      * 这里只初始化 metrics，不创建输出文件。
//      *
//      * 这样 metric log 里仍然可以看到所有 barcode，
//      * 但是不会因为初始化 writer 而生成 0B 文件。
//      */
//     for name in barcode_names.iter() {
//         let metric_barcode_name = normalize_barcode_name(name).to_string();

//         counter!(
//             "writer_barcode_count",
//             "barcode" => metric_barcode_name
//         )
//         .increment(0_u64);
//     }

//     counter!(
//         "writer_barcode_count",
//         "barcode" => "uncertain"
//     )
//     .increment(0_u64);

//     /*
//      * 关键修改：
//      *
//      * 不再使用全部 barcode_names 创建 writer。
//      * 先创建一个空 writers map。
//      *
//      * 后面只有真正遇到某个 barcode 时，才为它创建 writer。
//      */
//     let mut writers = make_writers(
//         output_folder,
//         Vec::<String>::new(),
//         input_file,
//         output_format,
//     )?;

//     match output_format {
//         "fa" | "fasta" | "fastq" | "fq" | "bam" => {
//             while let Ok(demuxed_record) = barcode_demux_info.recv() {
//                 counter!("writer_received").increment(1_u64);

//                 let (demuxed_reads, barcode_pair) = demuxed_record;

//                 /*
//                  * writer_key:
//                  *   用于输出文件 writer，保持原始 barcode 名字。
//                  *
//                  * metric_barcode_name:
//                  *   用于 metrics，使用 normalize_barcode_name。
//                  */
//                 let (writer_key, metric_barcode_name) = match barcode_pair.as_ref() {
//                     Some(barcode_pair) => {
//                         let raw_name = barcode_pair.name.to_string();
//                         let normalized_name = normalize_barcode_name(&raw_name).to_string();

//                         (raw_name, normalized_name)
//                     }
//                     None => ("uncertain".to_string(), "uncertain".to_string()),
//                 };

//                 /*
//                  * 懒创建 writer：
//                  *
//                  * 如果 writers 里还没有这个 barcode 对应的 writer，
//                  * 就只为当前这个 barcode 创建 writer。
//                  */
//                 let writer_lookup_key: String = if writers.contains_key(writer_key.as_str()) {
//                     writer_key.clone()
//                 } else if writers.contains_key(metric_barcode_name.as_str()) {
//                     metric_barcode_name.clone()
//                 } else {
//                     /*
//                      * 这里只创建当前要写入的这个 barcode 的文件。
//                      *
//                      * 不会再为全部 barcode_names 创建文件，
//                      * 因此不会产生大量 0B 文件。
//                      */
//                     let new_writers = make_writers(
//                         output_folder,
//                         vec![writer_key.clone()],
//                         input_file,
//                         output_format,
//                     )?;

//                     let new_lookup_key = if new_writers.contains_key(writer_key.as_str()) {
//                         writer_key.clone()
//                     } else if new_writers.contains_key(metric_barcode_name.as_str()) {
//                         metric_barcode_name.clone()
//                     } else {
//                         anyhow::bail!(
//                             "new writer missing key: raw_key={}, normalized_key={}",
//                             writer_key,
//                             metric_barcode_name
//                         );
//                     };

//                     writers.extend(new_writers);

//                     new_lookup_key
//                 };

//                 let w = writers.get_mut(writer_lookup_key.as_str()).ok_or_else(|| {
//                     anyhow::anyhow!(
//                         "writers missing key after lazy creation: {}",
//                         writer_lookup_key
//                     )
//                 })?;

//                 match w.write_demuxed_record(&demuxed_reads) {
//                     Ok(_) => {
//                         counter!(
//                             "writer_barcode_count",
//                             "barcode" => metric_barcode_name
//                         )
//                         .increment(1_u64);
//                     }
//                     Err(e) => {
//                         counter!("writer_failed").increment(1_u64);

//                         error!(
//                             ?e,
//                             barcode = metric_barcode_name.as_str(),
//                             writer_key = writer_key.as_str(),
//                             "demuxed read 写入失败"
//                         );
//                     }
//                 }

//                 pb.inc(1);
//             }
//         }
//         _ => anyhow::bail!("不支持的写入格式: {output_format}"),
//     }

//     pb.finish_and_clear();

//     info!("barcode writing finished; writers were created lazily");

//     Ok(())
// }
// use crate::core_context::BarcodePair;
// use crate::get_demuxed_reads::RecordType;
// use crate::io_utils::{make_writers, normalize_barcode_name};

// use crate::pbar::{DEFAULT_INTERVAL, get_spin_pb};
// use crossbeam::channel::Receiver;
// use metrics::counter;
// use tracing::{error, info};

// pub fn write_barcode_results(
//     input_file: &str,
//     output_format: &str,
//     output_folder: &str,
//     barcode_demux_info: Receiver<(RecordType, Option<BarcodePair>)>,
//     barcode_names: Vec<String>,
// ) -> anyhow::Result<()> {
//     println!("enter write process");

//     let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);

//     /*
//      * 初始化 barcode metric。
//      *
//      * 注意：
//      * 这里只用于 metrics 展示，所以使用 normalize_barcode_name。
//      * 这样 Single-1_Double-215_F / Single-1_Double-215_R
//      * 会统一显示成 Single-1_Double-215。
//      *
//      * increment(0) 的目的是让 0 count 的 barcode 也有机会出现在 metric log 中。
//      */
//     for name in barcode_names.iter() {
//         let metric_barcode_name = normalize_barcode_name(name).to_string();

//         counter!(
//             "writer_barcode_count",
//             "barcode" => metric_barcode_name
//         )
//         .increment(0_u64);
//     }

//     counter!(
//         "writer_barcode_count",
//         "barcode" => "uncertain"
//     )
//     .increment(0_u64);

//     let mut writers = make_writers(output_folder, barcode_names, input_file, output_format)?;

//     match output_format {
//         "fa" | "fasta" | "fastq" | "fq" | "bam" => {
//             while let Ok(demuxed_record) = barcode_demux_info.recv() {
//                 counter!("writer_received").increment(1_u64);

//                 let (demuxed_reads, barcode_pair) = demuxed_record;

//                 /*
//                  * writer_key:
//                  *   用于查找 writers，原则上保持原始 barcode 名称。
//                  *
//                  * metric_barcode_name:
//                  *   用于 metrics，必须 normalize。
//                  */
//                 let (writer_key, metric_barcode_name) = match barcode_pair.as_ref() {
//                     Some(barcode_pair) => {
//                         let raw_name = barcode_pair.name.to_string();
//                         let normalized_name = normalize_barcode_name(&raw_name).to_string();

//                         (raw_name, normalized_name)
//                     }
//                     None => ("uncertain".to_string(), "uncertain".to_string()),
//                 };

//                 /*
//                  * 优先用原始 writer_key 查找 writer。
//                  *
//                  * 这里额外加一个 fallback：
//                  * 如果 make_writers 内部已经使用 normalize_barcode_name 建立 key，
//                  * 那么也可以用 normalized name 找到对应 writer。
//                  */
//                 let writer_lookup_key = if writers.contains_key(writer_key.as_str()) {
//                     writer_key.as_str()
//                 } else if writers.contains_key(metric_barcode_name.as_str()) {
//                     metric_barcode_name.as_str()
//                 } else {
//                     anyhow::bail!(
//                         "writers missing key: raw_key={}, normalized_key={}",
//                         writer_key,
//                         metric_barcode_name
//                     );
//                 };

//                 let w = writers.get_mut(writer_lookup_key).ok_or_else(|| {
//                     anyhow::anyhow!("writers missing key after lookup: {}", writer_lookup_key)
//                 })?;

//                 match w.write_demuxed_record(&demuxed_reads) {
//                     Ok(_) => {
//                         // counter!("writer_ok").increment(1_u64);

//                         /*
//                          * barcode 统计只写入 metrics，不再额外生成 TSV。
//                          *
//                          * metric 中的 barcode 名称使用 normalize 后的名字。
//                          */
//                         counter!(
//                             "writer_barcode_count",
//                             "barcode" => metric_barcode_name
//                         )
//                         .increment(1_u64);
//                     }
//                     Err(e) => {
//                         counter!("writer_failed").increment(1_u64);

//                         error!(
//                             ?e,
//                             barcode = metric_barcode_name.as_str(),
//                             writer_key = writer_key.as_str(),
//                             "demuxed read 写入失败"
//                         );
//                     }
//                 }

//                 pb.inc(1);
//             }
//         }
//         _ => anyhow::bail!("不支持的写入格式: {output_format}"),
//     }

//     pb.finish_and_clear();

//     info!("barcode writing finished; barcode counts were recorded in metrics only");

//     Ok(())
// }
