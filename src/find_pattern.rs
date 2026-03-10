use crate::bam_record_extention::ReadRecord;
// use metrics::counter;
use std::collections::BTreeSet;
// use std::sync::Arc;

pub fn get_pattern_keys(pattern_records: &Vec<ReadRecord>) -> Vec<String> {
    let mut set: BTreeSet<String> = pattern_records
        .iter()
        .map(|r| {
            // 假设字段名是 r.name: Arc<str>（按你实际字段名改）
            let s: &str = r.id.as_ref();

            // 去掉尾部 2 个字符（不足 2 个就变成空串）
            let trimmed = if s.len() >= 2 { &s[..s.len() - 2] } else { "" };

            trimmed.to_string()
        })
        .collect();
    set.insert("uncertain".into());
    set.into_iter().collect()
}
// pub fn merge_non_overlapping_no_copy(
//     candidates: &mut Vec<BarcodeCandidate>,
// ) -> Vec<BarcodeCandidate> {
//     if candidates.is_empty() {
//         return Vec::new();
//     }

//     // 1. 先按 start 排序
//     candidates.sort_by_key(|c| c.start);

//     // 2. 用下标 i 来管理“合并后区间”末尾
//     let mut i = 0;
//     for j in 1..candidates.len() {
//         if candidates[j].start > candidates[i].end {
//             // 不重叠，指针 i 往后走一格
//             i += 1;
//             candidates[i] = candidates[j].clone();
//         } else {
//             counter!("primer overlap cnts").increment(1);
//             // 有重叠，若距离更小，则替换
//             if candidates[j].distance < candidates[i].distance {
//                 candidates[i] = candidates[j].clone();
//             } else if candidates[j].distance == candidates[i].distance {
//                 if candidates[j].end - candidates[j].start > candidates[i].end - candidates[i].start
//                 {
//                     // 如果距离相同，且新区间更长，则替换
//                     candidates[i] = candidates[j];
//                 }
//             }
//         }
//     }

//     // 3. truncate 保留前 i+1 个元素，即合并结果
//     candidates.truncate(i + 1);

//     // 4. 利用 std::mem::take “搬走”这 i+1 个元素，避免克隆
//     std::mem::take(candidates)
// }

#[cfg(test)]
mod test {
    use bio::alignment::Alignment;
    use bio::pattern_matching::myers::MyersBuilder;

    #[test]
    fn test_ambig_pattern() {
        let text: &'static [u8; 15] = b"GGATGAGCGCCATAG";
        let pattern = b"GAGGC";

        let mut myers = MyersBuilder::new().ambig(b'N', b"ACGT").build_64(pattern);
        let mut reuslt = myers.find_all(text, 2);
        let mut aln = Alignment::default();
        while reuslt.next_alignment(&mut aln) {
            println!(
                "start: {}, end: {}, distance: {}",
                aln.ystart, aln.yend, aln.score
            );
        }
    }
}
