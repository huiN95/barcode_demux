// use serde::Serialize;
use std::sync::Arc;

#[derive(Debug, Clone)]

pub struct BarcodeCandidate {
    pub start: usize,
    pub end: usize,
    pub distance: i32,
    pub name: Arc<str>,
}

#[derive(Debug, Clone)]

pub struct BarcodePair {
    pub name: Arc<str>,
    pub distance: (i32, i32),
    pub outter_position: (usize, usize),
    pub inner_position: (usize, usize),
}
// #[derive(Debug, Clone)]
// pub struct PrimerPosition {
//     pub primer_id: Arc<str>,
//     pub distance: Vec<i32>,
//     pub positions: Vec<(usize, usize)>,
//     pub single_end: bool,
// }

// #[derive(Debug, Clone, Serialize)]

// pub struct BarcodeMeta {
//     pub primer_id: String,
//     pub primer_distance: f32,
//     pub channel_idx: i32,
//     // pub single_end: bool,
// }

// #[derive(Clone, Debug)]
// pub struct Annotated<T> {
//     pub inner: T,
//     pub primer: Option<BarcodeMeta>,
// }

// pub type DemuxedResult = Vec<Annotated<PrimerPosition>>;
