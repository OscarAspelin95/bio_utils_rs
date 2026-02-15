use std::fs::File;
use std::path::PathBuf;
use tempfile::TempDir;

/// Creates a temporary directory with an empty file of the given name.
///
/// The returned [`TempDir`] must be kept alive for the file to remain on disk.
pub fn temp_seq_file(filename: &str) -> (TempDir, PathBuf) {
    let tmp_dir = TempDir::new().expect("failed to create temp dir");
    let tmp_file = tmp_dir.path().join(filename);
    File::create(&tmp_file).expect("failed to create temp file");

    (tmp_dir, tmp_file)
}
