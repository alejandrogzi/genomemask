<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/logo.png" >
  </p>

  <span>
    <h1 align="center">
        genomemask
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.1-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.1-green">
    </a>
    <a href="https://crates.io/crates/genomemask" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/genomemask">
    </a>
    <a href="https://github.com/alejandrogzi/genomemask" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/genomemask?color=blue">
    </a>
    <a href="https://crates.io/crates/genomemask" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/genomemask">
    </a>
  </p>

  <p align="center">
    <samp>
        <span>mask specific regions of a genome in any format</span>
        <br>
        <br>
        <a href="https://docs.rs/genomemask/0.0.1/genomemask/">docs</a> .
        <a href="https://github.com/alejandrogzi/genomemask?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/genomemask?tab=readme-ov-file#Installation">install</a>
    </samp>
  </p>

</p>

## Installation
### Binary
```bash
cargo install genomemask
```

### Docker
```bash
docker pull ghcr.io/alejandrogzi/genomemask:latest
```

### Conda
```bash
conda install -c bioconda genomemask
```

## Usage

### Quick Start

```bash
# Replace all N bases with 'A' in FASTA format
genomemask ns -s genome.fa -f fasta -n A

# Replace Ns in and output 2bit format
genomemask ns -s genome.2bit -f 2bit -n A
```

### Subcommands

#### `ns` - Replace N bases

Replace ambiguous 'N' bases in FASTA or 2bit genome files:

```bash
genomemask ns \
  --sequence genome.fa \
  --output-format 2bit \
  --nucleotide A 
```

#### `seleno` - Mask selenocysteine codons

Mask TGA codons at positions specified in a BED3 file:

```bash
genomemask seleno \
  --sequence genome.fa \
  --selenocysteine sites.bed \
  --output-format fasta \
  --nucleotide A 
```

#### `mask` - Mask genomic regions

Mask bases in specified genomic regions:

```bash
# Direct masking from BED/GTF/GFF intervals
genomemask mask \
  --sequence genome.fa \
  --regions regions.bed (can be GTF, GFF, or BED) \
  --feature cds \
  --output-format 2bit \
  --nucleotide T 
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-s, --sequence` | Input genome (.fa, .fa.gz, .fna, .fasta, .2bit) or stdin (-) | - |
| `-o, --outdir` | Output directory | . |
| `-f, --output-format` | Output format: fasta, 2bit, stdout | Required |
| `-n, --nucleotide` | Replacement base (A/T/C/G) | Required* |
| `-S, --stochastic` | Use random nucleotides with seed | - |
| `--seed` | Seed for deterministic random replacement | 0 |
| `-z, --gzip` | Compress FASTA output with gzip | - |
| `-t, --threads` | Number of worker threads | CPU count |
| `-l, --level` | Log level: error, warn, info, debug, trace | info |

*Required unless `--stochastic` is enabled
