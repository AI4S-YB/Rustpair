pub mod cli;
pub mod dedup;
pub mod format;
pub mod header;
pub mod io;
pub mod parse;
pub mod sort;

fn main() {
    if let Err(e) = cli::run() {
        eprintln!("Error: {}", e);
        for cause in e.chain() {
            eprintln!("  Caused by: {}", cause);
        }
        std::process::exit(1);
    }
}
