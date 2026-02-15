use bio::pattern_matching::myers::MyersBuilder;
use memchr::memmem;

#[inline]
fn myers_builder(primer_seq: &[u8]) -> bio::pattern_matching::myers::Myers {
    MyersBuilder::new()
        .ambig(b'N', b"ACGT")
        .ambig(b'R', b"AG")
        .ambig(b'Y', b"CT")
        .ambig(b'S', b"GC")
        .ambig(b'W', b"AT")
        .ambig(b'K', b"GT")
        .ambig(b'M', b"AC")
        .ambig(b'B', b"CGT")
        .ambig(b'D', b"AGT")
        .ambig(b'H', b"ACT")
        .ambig(b'V', b"ACG")
        .build_64(primer_seq)
}

/// TODO - fix and improve this.
fn search_fuzzy(seq: &[u8], barcode: &[u8], max_mismatches: u8) {
    let mut myers = myers_builder(barcode);

    // do something with this.
    let s = myers.find_all_lazy(seq, max_mismatches);
}

/// TODO - fix and improve this.
fn search_exact(seq: &[u8], barcode: &[u8]) {
    /// do something with this.
    let hits: Vec<usize> = memmem::find_iter(seq, barcode).collect();
}
