use crate::bam_record_extention::ReadRecord;
use crate::core_context::BarcodeCandidate;
use anyhow::Result;
use bio::alignment::Alignment;
use bio::pattern_matching::myers::long::Myers as MyersLong;
use bio::pattern_matching::myers::Myers as Myers64;
use bio::pattern_matching::myers::MyersBuilder;
use std::collections::HashMap;
use std::sync::Arc;
pub enum MayersPattern {
    Myers64 {
        myers: Myers64<u64>,
        pattern: String,
    },
    MyersLong {
        myers: MyersLong<u8>,
        pattern: String,
    },
}

impl MayersPattern {
    #[inline]
    pub fn pattern(&self) -> &str {
        match self {
            MayersPattern::Myers64 { pattern, .. } => pattern,
            MayersPattern::MyersLong { pattern, .. } => pattern,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub enum Direction {
    Forward,
    Reverse,
}

const AMBIGS: &[(u8, &[u8])] = &[
    (b'M', b"AC"),
    (b'R', b"AG"),
    (b'W', b"AT"),
    (b'S', b"CG"),
    (b'Y', b"CT"),
    (b'K', b"GT"),
    (b'V', b"ACG"),
    (b'H', b"ACT"),
    (b'D', b"AGT"),
    (b'B', b"CGT"),
    (b'N', b"ACGT"),
];

pub fn barcode_alignment(
    name: Arc<str>,
    myers: &mut MayersPattern,
    target: &[u8],
    max_distance: u8,
) -> Result<Option<BarcodeCandidate>> {
    let mut aln = Alignment::default();
    match myers {
        MayersPattern::Myers64 {
            myers,
            pattern: _pat,
        } => {
            let mut matches = myers.find_all_lazy(target, max_distance.into());
            let Some((best_end, _)) = matches.by_ref().min_by_key(|&(_, dist)| dist) else {
                return Ok(None); // 没匹配到，正常返回
            };
            matches.alignment_at(best_end, &mut aln);
            // 假设 BarcodeCandidate 里有这些字段（按你真实定义改字段名）
            let barcode_result = BarcodeCandidate {
                start: aln.ystart,
                end: aln.yend,
                distance: aln.score,
                name: name, // 如果你结构体里需要距离（推荐保留）
            };

            Ok(Some(barcode_result))
        }

        MayersPattern::MyersLong { myers, pattern: _ } => {
            let mut matches = myers.find_all_lazy(target, max_distance.into());
            let Some((best_end, _)) = matches.by_ref().min_by_key(|&(_, dist)| dist) else {
                return Ok(None); // 没匹配到，正常返回
            };
            matches.alignment_at(best_end, &mut aln);
            // 假设 BarcodeCandidate 里有这些字段（按你真实定义改字段名）
            let barcode_result = BarcodeCandidate {
                start: aln.ystart,
                end: aln.yend,
                distance: aln.score,
                name: name, // 如果你结构体里需要距离（推荐保留）
            };

            Ok(Some(barcode_result))
        }
    }
}

// pub fn get_seq_reverse_complement(seq: &[u8]) -> Vec<u8> {
//     seq.iter()
//         .rev()
//         .map(|&base| match base {
//             b'A' => b'T',
//             b'T' => b'A',
//             b'C' => b'G',
//             b'G' => b'C',
//             _ => base, // 保持其他字符不变
//         })
//         .collect()
// }

/// 统一把 builder 的结果包成 MayersType
#[inline]
fn build_myers_enum(builder: &MyersBuilder, pat: &[u8]) -> MayersPattern {
    // 规则：短到 64 的用 Myers64，更长的用 MyersLong
    let pattern = String::from_utf8_lossy(pat).into_owned();

    if pat.len() <= 64 {
        let m64 = builder.build_64(pat); // ← 按你的实际 API 名字改
        MayersPattern::Myers64 {
            myers: m64,
            pattern,
        }
    } else {
        let ml = builder.build_long(pat); // ← 按你的实际 API 名字改
        MayersPattern::MyersLong { myers: ml, pattern }
    }
}

fn compl_iupac(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => b,
    }
}

fn revcomp_iupac(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(compl_iupac(b.to_ascii_uppercase()));
    }
    out
}

// pub fn get_myers_from_barcodes(
//     patterns: &[ReadRecord],
//     direction: Direction,
// ) -> HashMap<Arc<str>, MayersPattern> {
//     let mut builder = MyersBuilder::new();
//     for &(base, equivalents) in AMBIGS {
//         builder.ambig(base, equivalents);
//     }

//     let suffix = match direction {
//         Direction::Forward => "_F",
//         Direction::Reverse => "_R",
//     };

//     patterns
//         .iter()
//         .filter_map(|r| {
//             let id = r.id.as_ref(); // &str
//             if !id.ends_with(suffix) {
//                 return None;
//             }

//             let pattern = r.sequence.as_slice();
//             let myers = build_myers_enum(&builder, pattern); // 假设它会返回 MayersPattern

//             Some((Arc::<str>::from(id.to_owned()), myers))
//             // 如果 r.id 本身就是 Arc<str>，更推荐：Some((r.id.clone(), myers))
//         })
//         .collect()
// }

pub fn get_myers_from_barcodes(
    patterns: &[ReadRecord],
    direction: Direction,
) -> HashMap<Arc<str>, MayersPattern> {
    let mut builder = MyersBuilder::new();
    for &(base, equivalents) in AMBIGS {
        builder.ambig(base, equivalents);
    }

    let suffix = match direction {
        Direction::Forward => "_F",
        Direction::Reverse => "_R",
    };

    patterns
        .iter()
        .filter_map(|r| {
            let id = r.id.as_ref(); // &str

            // 既做过滤，也拿到去掉后缀后的 key
            let key = id.strip_suffix(suffix)?;

            let pattern = r.sequence.as_slice();
            let myers = build_myers_enum(&builder, pattern);

            Some((Arc::<str>::from(key), myers))
            // 如果你想避免分配，也可以：Arc::<str>::from(key.to_owned())
        })
        .collect()
}

pub fn get_two_directions_myers_from_barcodes(
    patterns: &[ReadRecord],
    direction: Direction,
) -> HashMap<Arc<str>, MayersPattern> {
    let mut builder = MyersBuilder::new();
    for &(base, equivalents) in AMBIGS {
        builder.ambig(base, equivalents);
    }

    let suffix = match direction {
        Direction::Forward => "_F",
        Direction::Reverse => "_R",
    };

    patterns
        .iter()
        .filter_map(|r| {
            let id = r.id.as_ref(); // &str

            // 既做过滤，也拿到去掉后缀后的 key
            let key = id.strip_suffix(suffix)?;

            let pattern = r.sequence.as_slice();
            let myers = build_myers_enum(&builder, pattern);

            Some((Arc::<str>::from(key), myers))
            // 如果你想避免分配，也可以：Arc::<str>::from(key.to_owned())
        })
        .collect()
}
