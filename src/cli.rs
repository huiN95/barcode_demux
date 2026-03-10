
use clap::Parser;
// use tracing_subscriber::{EnvFilter, fmt, layer::SubscriberExt, util::SubscriberInitExt};

// / 定义 CLI 结构体
#[derive(Parser, Debug)]
#[command(
    name = "artificial barcode demux",
    version = "0.0.1",
    about = "匹配给定模式在序列文件中的近似匹配"
)]
pub struct Cli {

    #[arg(long = "pipeline_version", default_value ="1")]
    pub pipeline_version: String,

    /// 指定 primer 序列文件 (FASTA 格式)
    #[arg(short = 'b', long = "barcode")]
    pub barcode: String,


    // 匹配pattern时允许的最大错误率
    #[arg(short = 'd', long = "max_distance",
          value_parser = clap::value_parser!(u8),
          default_value_t = 0)]
    pub max_distance: u8,


    /// 指定输入序列文件 (FASTA/FASTQ/BAM)
    #[arg(short = 'i', long = "input_file")]
    pub input_file: String,

    /// 设置日志级别（可选），例如 "info", "debug", "trace"
    #[arg(long, env = "RUST_LOG")]
    pub log: Option<String>,

    #[arg(short = 'o', long = "output_folder")]
    pub output_folder: String,
    
    #[arg(long = "log_folder")]
    pub log_folder: String,

    #[arg(long = "reservesed_threads",
          value_parser = clap::value_parser!(u8).range(0..=255),
          default_value_t = 2)]
    pub reservesed_threads:u8,
    
    #[arg(long = "output_format", default_value ="fasta" ,       
        value_parser = ["fasta", "fastq", "bam","fq", "fa"])]
    pub output_format: String,
   
    /// 最小子读长度，默认50,如果是保留primer，是包含在里面的，需要适量增大该参数
    #[arg(short = 'l', 
    long = "min_subread_len", 
    default_value_t = 50,
    value_parser = clap::value_parser!(usize),)]
    pub min_subread_len: usize,

    /// 只处理超过该阈值的reads,进行demux，否则pass。尚未实现。
    #[arg(long = "min_q",
          value_parser = clap::value_parser!(usize),
          default_value_t = 20)]
    pub min_q: usize,

        #[arg(long = "search_bound",
          value_parser = clap::value_parser!(usize),
          default_value_t = 40)]
    pub search_bound: usize,

}
