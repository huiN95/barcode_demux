use clap::Parser;

/// 定义 CLI 结构体
#[derive(Parser, Debug)]
#[command(
    name = "barcode demux",
    version = "1.0.7",
    about = "匹配给定模式在序列文件中的近似匹配",
    long_about = None
)]
pub struct Cli {
    /// 指定使用的处理流程版本号
    #[arg(
        long = "pipeline-version",
        alias = "pipeline_version",
        default_value = "1"
    )]
    pub pipeline_version: String,

    /// 指定 primer/barcode 序列文件 (FASTA 格式)
    #[arg(short = 'b', long = "barcode")]
    pub barcode: String,

    /// 匹配 pattern 时允许的最大错误率（编辑距离）
    #[arg(
        short = 'd',
        long = "max-distance",
        alias = "max_distance",
        value_parser = clap::value_parser!(u8),
        default_value_t = 0
    )]
    pub max_distance: u8,

    /// 指定输入序列文件 (支持 FASTA/FASTQ/BAM 格式)
    #[arg(short = 'i', long = "input-file", alias = "input_file")]
    pub input_file: String,

    /// 设置日志级别（可选），例如 "info", "debug", "trace"
    #[arg(long, env = "RUST_LOG", default_value = "info")]
    pub log: Option<String>,

    /// 指定输出结果的保存目录
    #[arg(short = 'o', long = "output-folder", alias = "output_folder")]
    pub output_folder: String,

    /// 指定运行日志文件的保存目录
    #[arg(long = "log-folder", alias = "log_folder")]
    pub log_folder: String,

    /// 设置运行使用的线程数 (0 表示自动根据 CPU 核心数分配)
    #[arg(
        long = "threads",
        value_parser = clap::value_parser!(u16).range(0..=1024),
        default_value_t = 0
    )]
    pub threads: u16,

    /// 指定输出序列文件的格式
    #[arg(
        long = "output-format",
        alias = "output_format",
        default_value = "fasta",
        value_parser = ["fasta", "fastq", "bam", "fq", "fa"]
    )]
    pub output_format: String,

    /// 最小子读长度，默认50。如果是保留primer，primer是包含在里面的，需要适量增大该参数
    #[arg(
        short = 'l',
        long = "min-subread-len",
        alias = "min_subread_len",
        default_value_t = 50,
        value_parser = clap::value_parser!(usize)
    )]
    pub min_subread_len: usize,

    /// 质量控制：只处理超过该阈值的 reads 进行 demux，否则 pass (尚未实现)
    #[arg(
        long = "min-q",
        alias = "min_q",
        value_parser = clap::value_parser!(usize),
        default_value_t = 20
    )]
    pub min_q: usize,

    /// 在序列两端搜索 barcode 时的边界范围 (bp)
    #[arg(
        long = "search-bound",
        alias = "search_bound",
        value_parser = clap::value_parser!(usize),
        default_value_t = 40
    )]
    pub search_bound: usize,

    /// 拆分后的序列是否保留 barcode，仅在 pipeline 3 模式下生效
    #[arg(long = "keep-barcode", alias = "keep_barcode", default_value_t = false)]
    pub keep_barcode: bool,

    /// 是否过滤掉只有单端匹配到 barcode 的序列
    #[arg(
        long = "single-end-filter",
        alias = "single_end_filter",
        default_value_t = false
    )]
    pub single_end_filter: bool,

    /// 最小 pair 长度（仅在 pipeline_version=4 时需要）
    #[arg(
        long = "min-pair-len",
        alias = "min_pair_len",
        required_if_eq("pipeline_version", "4"),
        value_parser = clap::value_parser!(usize)
    )]
    pub min_pair_len: Option<usize>,

    /// 最小 pair 匹配得分（仅在 pipeline_version=4 时需要）
    #[arg(
        long = "min-pair-score",
        alias = "min_pair_score",
        required_if_eq("pipeline_version", "4"),
        value_parser = clap::value_parser!(f64)
    )]
    pub min_pair_score: Option<f64>,
}
