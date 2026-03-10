use crate::cli::Cli;
use metrics;
use metrics_exporter_prometheus::{PrometheusBuilder, PrometheusHandle};
// use once_cell::sync::Lazy;
// use std::sync::Mutex;
// use std::thread;
// use std::time::Duration;
use std::{
    error::Error,
    fs::{create_dir_all, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};
use tracing_subscriber;
// static LOG_GUARD: Lazy<Mutex<Option<tracing_appender::non_blocking::WorkerGuard>>> =
//     Lazy::new(|| Mutex::new(None));
// use std::{fs, path::Path};
// use tracing_subscriber::fmt::time::ChronoLocal;
use tracing_subscriber::{fmt, layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

pub fn init_tracing_log(cli: &Cli) -> tracing_appender::non_blocking::WorkerGuard {
    let log_path = Path::new(&cli.log_folder);

    // 如果不存在就递归创建
    if !log_path.exists() {
        create_dir_all(log_path).unwrap();
    }
    let run_log_path = log_path.join(format!("{}_primer_demux_run.log", cli.log_folder));

    // let run_log_filepath = format!("{}_primer_demux_run.log", log_path);
    let log_file = File::create(&run_log_path).unwrap_or_else(|e| {
        panic!(
            "create run log file error: {} ({})",
            run_log_path.display(),
            e
        )
    });
    // print!("run log file: {}", &run_log_filepath);
    let (non_blocking, _guard) = tracing_appender::non_blocking(log_file);
    let filter = if let Some(s) = &cli.log {
        EnvFilter::try_new(s).unwrap()
    } else {
        // 没传参数、也没环境变量时的默认
        EnvFilter::new("info")
    };
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .with_timer(tracing_subscriber::fmt::time::ChronoLocal::rfc_3339())
                .with_ansi(false)
                .with_writer(non_blocking)
                .with_target(true)
                .with_thread_ids(true),
        )
        .try_init()
        .ok();
    // tracing_subscriber::fmt()
    //     // .with_timer(UtcTime::rfc_3339())
    //     .with_timer(ChronoLocal::rfc_3339())
    //     .with_ansi(false) //去掉颜色信息
    //     .with_writer(non_blocking)
    //     .try_init()
    //     .expect("tracing_subscriber already initialized!");
    _guard
    // *LOG_GUARD.lock().unwrap() = Some(_guard);
}

pub struct MetricsGuard {
    handle: PrometheusHandle,
    path: PathBuf,
}

impl Drop for MetricsGuard {
    fn drop(&mut self) {
        // 收集指标
        let report = self.handle.render();

        // 尝试写入；如果失败，只打印日志，不 panic
        if let Err(e) = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.path)
            .and_then(|mut f| f.write_all(report.as_bytes()))
        {
            eprintln!("metrics dump failed: {e}");
        }
    }
}

/// 在程序启动时调用，返回一个 MetricsGuard。
pub fn init_metrics<P: AsRef<Path>>(dir: P) -> Result<MetricsGuard, Box<dyn Error>> {
    // 1. 构建 recorder
    let recorder = PrometheusBuilder::new()
        .set_buckets(
            vec![
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 13.0, 21.0, 34.0, 55.0, 100.0,
                1000.0,
            ]
            .as_slice(),
        )
        .unwrap()
        .build_recorder();
    let handle = recorder.handle();
    // 2. 注册为全局。若已注册过，返回 Err，而不是 panic。
    metrics::set_global_recorder(Box::new(recorder)).unwrap();

    // 3. 生成日志路径
    // let path = dir.as_ref().join("primer_metrics.prom"); // 可根据需要改名/加时间戳
    // let mut f = std::fs::File::create(&path).unwrap();
    let in_path = dir.as_ref();
    let out_path = if in_path.is_dir() {
        // 传进来的是目录
        in_path.join("_primer_metrics.log")
    } else {
        // 传进来的是文件；改文件名
        let stem = in_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("metrics");

        let mut buf = PathBuf::from(in_path);
        buf.set_file_name(format!("{stem}_primer_metrics.log"));
        buf
    };
    Ok(MetricsGuard {
        handle,
        path: out_path,
    })
}
