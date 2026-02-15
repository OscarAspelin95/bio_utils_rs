use std::path::PathBuf;
use tempfile::TempDir;

/// Creates a temporary directory with an empty FASTQ file.
/// Returns the TempDir (must be kept alive) and the path to the file.
pub fn temp_fastq_file(filename: &str) -> (TempDir, PathBuf) {
    let tmp_dir = TempDir::new().expect("failed to create temp dir");
    let tmp_file = tmp_dir.path().join(filename);
    std::fs::File::create(&tmp_file).expect("failed to create temp file");
    (tmp_dir, tmp_file)
}

/// Creates a temporary directory with an empty FASTA file.
/// Returns the TempDir (must be kept alive) and the path to the file.
pub fn temp_fasta_file(filename: &str) -> (TempDir, PathBuf) {
    let tmp_dir = TempDir::new().expect("failed to create temp dir");
    let tmp_file = tmp_dir.path().join(filename);
    std::fs::File::create(&tmp_file).expect("failed to create temp file");
    (tmp_dir, tmp_file)
}

/// Creates a temporary directory with a file of any name.
/// Returns the TempDir (must be kept alive) and the path to the file.
pub fn temp_seq_file(filename: &str) -> (TempDir, PathBuf) {
    let tmp_dir = TempDir::new().expect("failed to create temp dir");
    let tmp_file = tmp_dir.path().join(filename);
    std::fs::File::create(&tmp_file).expect("failed to create temp file");
    (tmp_dir, tmp_file)
}
