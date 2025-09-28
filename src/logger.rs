use env_logger::{fmt::Formatter, Builder};
use log::{LevelFilter, Record};
use std::{env, io::Write, path::Path};

pub fn init_logger() {
    env::set_var("RUST_LOG", "info");
    // Builder::new()
    //     .filter_level(LevelFilter::Info)
    //     .format(|buf: &mut Formatter, record: &Record| {
    //         // record.file() 可能是 Option<&str>
    //         let file = record.file().unwrap_or("unknown");
    //         let filename = Path::new(file)
    //             .file_name()
    //             .and_then(|f| f.to_str().map(|s| s.trim_end_matches(".rs")))
    //             .unwrap_or(file);

    //         writeln!(
    //             buf,
    //             "[ {} Isopedia::{}] {}",
    //             filename,
    //             record.level(),
    //             record.args()
    //         )
    //     })
    //     .init();
    //     Builder::new()
    //         .filter_level(LevelFilter::Info)
    //         .format(|buf: &mut Formatter, record: &Record| {
    //             // record.file() 可能是 Option<&str>
    //             let file = record.file().unwrap_or("unknown");
    //             let filename = Path::new(file)
    //                 .file_name()
    //                 .and_then(|f| f.to_str().map(|s| s.trim_end_matches(".rs")))
    //                 .unwrap_or(file);

    //             writeln!(
    //                 buf,
    //                 "Isopedia [{}] [{}] {}",
    //                 filename,
    //                 record.level(),
    //                 record.args()
    //             )
    //         })
    //         .init();
    //     Builder::new()
    //         .filter_level(LevelFilter::Info)
    //         .format(|buf: &mut Formatter, record: &Record| {
    //             // record.file() 可能是 Option<&str>
    //             let file = record.file().unwrap_or("unknown");
    //             let filename = Path::new(file)
    //                 .file_name()
    //                 .and_then(|f| f.to_str().map(|s| s.trim_end_matches(".rs")))
    //                 .unwrap_or(file);

    //             writeln!(
    //                 buf,
    //                 "Isopedia [{}] [{}] {}",
    //                 filename,
    //                 record.level(),
    //                 record.args()
    //             )
    //         })
    //         .init();
    env_logger::init();
}

// use env_logger::{Builder, fmt::Formatter};
// use log::{Record, LevelFilter};
// use std::{io::Write, path::Path};
// use chrono::Local; // 用于时间戳

// pub fn init_logger() {
//    env_logger::init();
// }
