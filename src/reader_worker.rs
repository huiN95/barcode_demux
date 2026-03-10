use crate::bam_record_extention::AuxExt;
use crate::bam_record_extention::ReadRecord;
use bio::io::fastq::Reader as FastqReader;
use crossbeam::channel::Sender;
// use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Read, Reader as BAMReader}; // <─ import the trait!
use std::sync::Arc;
// use rust_htslib::htslib::{ hts_readrec_func};
use seq_io::fasta::{Reader, Record};
use std::fs::File;
use std::{path::Path, path::PathBuf, thread};
use walkdir::WalkDir;
pub fn read_fasta_2_queue<P: AsRef<Path>>(
    input_path: P,
    sender: Sender<ReadRecord>,
) -> Result<(), Box<dyn std::error::Error>> {
    // let file = File::open(&input_path).unwrap();
    // let buf_size = 8 * 1024 * 1024; // 8MB
    // let mut reader = Reader::with_capacity(file, buf_size);

    let mut reader = Reader::from_path(&input_path).unwrap();

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let mut read_record = ReadRecord::default();
        read_record.id = record.id().unwrap().into();
        read_record.sequence = record.seq().to_vec();
        read_record.quality = None;

        // let read_record = ReadRecord {
        //     id: record.id().unwrap().to_string(),
        //     sequence: record.seq().to_vec(),
        //     quality: None,
        // };
        // let read_record_clone = read_record.clone();
        if let Err(e) = sender.send(read_record) {
            // print!("Failed to send: {:?}", read_record_clone);
            eprintln!("Failed to send: {}", e);
            break;
            // Optionally return early or break here
        }
    }
    Ok(())
}

pub fn read_fastq_2_queue<P: AsRef<Path>>(
    input_path: P,
    sender: Sender<ReadRecord>,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(input_path)?;
    let reader = FastqReader::new(file);
    for record in reader.records() {
        let record = record?;
        let mut read_record = ReadRecord::default();
        read_record.id = record.id().into();
        read_record.sequence = record.seq().to_vec();
        read_record.quality = Some(record.qual().to_vec());
        if let Err(e) = sender.send(read_record) {
            eprintln!("Failed to send: {}", e);
            break;
        }
    }
    Ok(())
}

fn read_bam_2_queue<P: AsRef<Path>>(
    input_path: P,
    sender: Sender<ReadRecord>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = BAMReader::from_path(&input_path).unwrap();
    reader.set_threads(4)?; // 解压 4 线程

    for record in reader.records() {
        let record = record?;
        let mut read_record = ReadRecord::default();
        read_record.id = Arc::<str>::from(String::from_utf8(record.qname().to_vec()).unwrap());
        read_record.sequence = record.seq().as_bytes().to_vec();
        let qual_ascii: Vec<u8> = record // &bam::Record
            .qual() // &[u8]（裸 PHRED）
            .iter()
            .map(|&q| if q == 255 { b'!' } else { q }) // 255 = “缺失质量”
            .collect();
        read_record.quality = Some(qual_ascii);
        // if let Some(dw) = record.array_u8(b"dw") {
        //     read_record.dw = Some(dw);
        // }
        // if let Some(cr) = record.array_u8(b"cr") {
        //     read_record.cr = Some(cr);
        // }
        if let Some(rq) = record.float(b"rq") {
            read_record.rq = Some(rq);
        }
        if let Some(np) = record.i32(b"np") {
            read_record.np = Some(np);
        }
        // if let Some(cx) = record.i32(b"cx") {
        //     read_record.cx = Some(cx);
        // }
        // if let Ok(a) = record.aux(b"ch") {
        //     eprintln!("raw ch aux = {:?}", a);
        // } else {
        //     eprintln!("no ch tag");
        // }
        if let Some(ch) = record.i32(b"ch") {
            read_record.ch = Some(ch);
        }

        if let Some(cq) = record.float(b"cq") {
            read_record.cq = Some(cq);
        }
        // if let Some(ar) = record.array_u32(b"ar") {
        //     read_record.ar = Some(ar);
        // }

        if let Err(e) = sender.send(read_record) {
            eprintln!("Failed to send: {}", e);
            break;
        }
    }
    Ok(())
}

pub fn read_bam_folder(
    bam_file_path: &Path,
    sender: Sender<ReadRecord>,
) -> Result<(), Box<dyn std::error::Error>> {
    let bam_files = collect_bam_files(bam_file_path);
    for bam_file in bam_files {
        let sender_clone = sender.clone();
        let bam_file_path_clone = bam_file.clone();
        thread::spawn(move || {
            if let Err(e) = read_bam_2_queue(bam_file_path_clone.to_str().unwrap(), sender_clone) {
                eprintln!("Error reading BAM file: {}", e);
            }
        });
    }
    Ok(())
    // 等待所有线程完成
}

pub fn collect_bam_files<P: Into<PathBuf>>(dir: P) -> Vec<PathBuf> {
    WalkDir::new(dir.into())
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "bam"))
        .map(|e| e.into_path())
        .collect()
}

pub fn read_sequences_to_queue<P: AsRef<Path>>(
    input_path: P,
    sender: Sender<ReadRecord>,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = input_path.as_ref();

    // === ① 先判断是不是目录 ===
    if path.is_dir() {
        // 这里把 `Path` 转回 `&Path` 传给 folder 函数
        return read_bam_folder(path, sender);
    }

    // === ② 再按扩展名区分文件格式 ===
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .ok_or("unsupported input format")?
        .to_ascii_lowercase();

    match ext.as_str() {
        "fa" | "fasta" => read_fasta_2_queue(path, sender),
        "fq" | "fastq" => read_fastq_2_queue(path, sender),
        "bam" => read_bam_2_queue(path, sender),
        _ => Err("unsupported input format".into()),
    }
}
