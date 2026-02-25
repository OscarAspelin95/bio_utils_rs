# bio_utils_rs

[![Crates.io](https://img.shields.io/crates/v/bio_utils_rs.svg)](https://crates.io/crates/bio_utils_rs)
[![Docs.rs](https://docs.rs/bio_utils_rs/badge.svg)](https://docs.rs/bio_utils_rs)
[![CI](https://github.com/OscarAspelin95/bio_utils_rs/actions/workflows/ci.yaml/badge.svg)](https://github.com/OscarAspelin95/bio_utils_rs/actions/workflows/ci.yaml)
[![License](https://img.shields.io/crates/l/bio_utils_rs.svg)](https://github.com/OscarAspelin95/bio_utils_rs/blob/main/LICENSE)

A Rust library of bioinformatics utilities covering sequence I/O, nucleotide operations, amino acid translation, k-mer sketching, and SIMD-accelerated sequence indexing.

## Features

| Module | Description | Feature flag |
|---|---|---|
| `nucleotide` | Reverse complement, GC content, quality metrics, entropy, homopolymer detection, pattern search | _(always available)_ |
| `aminoacid` | Nucleotide-to-amino-acid translation with configurable codon tables and reading frames | _(always available)_ |
| `kmers` | FracMinHash sketching over canonical k-mers | _(always available)_ |
| `io` | FASTQ/FASTA readers and writers for plain and gzip-compressed files, with stdin/stdout support | `io` |
| `simd_sketch` | SIMD-accelerated minimizer and syncmer sketching, parallel reverse index construction | `simd` |

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
bio_utils_rs = "0.0.2"
```

To enable optional modules:

```toml
[dependencies]
bio_utils_rs = { version = "0.0.2", features = ["io", "simd"] }
```

## Usage

### Nucleotide operations

```rust
use bio_utils_rs::nucleotide::{
    reverse_complement, gc_content, nucleotide_counts,
    mean_error_and_phred, shannon_entropy, find_homopolymers,
    search_exact, search_fuzzy,
};

// Reverse complement (supports all IUPAC ambiguity codes)
let rc = reverse_complement(b"AACGTT");
assert_eq!(rc, b"AACGTT");

// GC content
let gc = gc_content(b"ATCG"); // 0.5

// Nucleotide composition
let (counts, softmasked, ambiguous) = nucleotide_counts(b"AACGttNN");
// counts = [A, C, G, T], softmasked = lowercase bases, ambiguous = N/other

// Mean Phred quality from a quality string
let (mean_error, mean_phred) = mean_error_and_phred(b"IIIII");

// Shannon entropy
let entropy = shannon_entropy(b"ACGTACGT");

// Homopolymer detection
let runs = find_homopolymers(b"AAACCCGGG", 3);

// Exact pattern search (SIMD-accelerated)
let hits = search_exact(b"ACGTACGTACGT", b"ACGT");
// hits = [0, 4, 8]

// Fuzzy pattern search (Myers bit-parallel, IUPAC-aware)
let hits = search_fuzzy(b"AACCTAA", b"CGT", 2);
```

### Amino acid translation

```rust
use bio_utils_rs::aminoacid::{translate, codon_table::CodonTable};

// Translate in reading frame 1 (zero-offset)
// Frame is re-exported from the aminoacid module
```

### K-mer sketching (FracMinHash)

```rust
use bio_utils_rs::kmers::frac_min_hash;

// Produce a FracMinHash sketch of canonical k-mers
// kmer_size = 21, downsampling factor = 10
let sketch = frac_min_hash(21, 10, b"ACGTACGTACGTACGTACGTACGT")?;
```

The sketch uses 2-bit packed `u64` k-mer encodings and retains canonical (strand-symmetric) hashes falling below `u64::MAX / ds_factor`.

### Sequence I/O (`io` feature)

```rust
use bio_utils_rs::io::{bio_fastq_reader, bio_fasta_reader, needletail_reader};
use std::path::PathBuf;

// Read a FASTQ file (plain or gzip)
let reader = bio_fastq_reader(Some(PathBuf::from("reads.fastq.gz")))?;

// Read from stdin
let reader = bio_fastq_reader(None)?;

// Auto-detect FASTA/FASTQ format with needletail
let reader = needletail_reader(Some(PathBuf::from("sequences.fa")))?;
```

Supported extensions: `.fastq`, `.fq`, `.fasta`, `.fa` — all optionally gzip-compressed (`.gz`).

### SIMD sketching and indexing (`simd` feature)

```rust
use bio_utils_rs::simd_sketch::{
    MinimizerSketch, OpenSyncmerSketch, ClosedSyncmerSketch,
    build_reverse_index, Sketcher,
};

// Minimizer sketch
let sketcher = MinimizerSketch { kmer_size: 21, window_size: 11 };
let hashes = sketcher.sketch(b"ACGTACGTACGT...");

// Open syncmer sketch
let sketcher = OpenSyncmerSketch { kmer_size: 21, smer_size: 11, offset: 0 };

// Closed syncmer sketch
let sketcher = ClosedSyncmerSketch { kmer_size: 21, smer_size: 11 };

// Build a parallel reverse index mapping hash -> sequence bitset
let seqs: Vec<&[u8]> = vec![b"ACGT...", b"TGCA..."];
let index = build_reverse_index(&seqs, &sketcher);
```

The reverse index is built in parallel with Rayon and stored in a `DashMap<u64, FixedBitSet>`, enabling efficient sequence lookup by shared k-mer hashes.

## Feature flags

| Flag | Enables | Additional dependencies |
|---|---|---|
| `io` | `io` module — FASTQ/FASTA readers and writers | `flate2`, `needletail`, `serde`, `serde_json` |
| `simd` | `simd_sketch` module — SIMD minimizer/syncmer sketching and parallel reverse index | `simd-minimizers`, `packed-seq`, `dashmap`, `fixedbitset`, `rustc-hash`, `rayon` |

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
