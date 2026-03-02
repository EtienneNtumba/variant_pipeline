# 🧬 Germline Variant Calling Pipeline

**FastQC → fastp → BWA-MEM2 → GATK BQSR → HaplotypeCaller → SnpEff**  
Based on [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) for germline short variant discovery.

---

## Pipeline Overview

```
Raw FASTQ (R1 + R2)
    │
    ├──► [1] FastQC (raw)           → QC reports (HTML + ZIP)
    │
    ├──► [2] fastp                  → Trimmed reads + QC JSON/HTML
    │         │
    │         └──► [3] FastQC (trimmed)
    │
    ▼
[4] BWA-MEM2 index (reference, once)
    │
    ▼
[5] BWA-MEM2 align                  → Unsorted BAM (with @RG tags)
    │
    ▼
[6] samtools sort                   → Coordinate-sorted BAM
    │
    ▼
[7] samtools index                  → BAM + .bai index
    │
    ▼
[8] GATK MarkDuplicates             → Deduplicated BAM + metrics
    │
    ├──► [9]  samtools flagstat     → Alignment stats
    ├──► [10] samtools coverage     → Per-base coverage
    │
    ▼
[11] GATK BaseRecalibrator          → Recalibration table (BQSR step 1)
    │
    ▼
[12] GATK ApplyBQSR                 → Recalibrated BAM (BQSR step 2)
    │
    ▼                               (runs per sample in parallel)
[13] GATK HaplotypeCaller (gVCF)   → per-sample .g.vcf.gz
    │
    ▼
[14] GATK CombineGVCFs              → merged cohort .g.vcf.gz
    │
    ▼
[15] GATK GenotypeGVCFs             → raw joint-called VCF
    │
    ├──► [16] SelectVariants SNP   ──► [18] VariantFiltration SNP ──┐
    │                                                                 ├──► [20] MergeVcfs
    └──► [17] SelectVariants INDEL ──► [19] VariantFiltration INDEL ┘
                                                 │
                                                 ▼
                                        [21] SnpEff annotation
                                                 │
                                                 ▼
                                        [22] Summary (TSV + PDF)
                                                 │
                                                 ▼
                                        [23] MultiQC (aggregated QC)
```

---

## Requirements

### Tools
| Tool | Version | Purpose |
|------|---------|---------|
| [Nextflow](https://nextflow.io) | ≥ 23.04 | Pipeline orchestration |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥ 0.12 | Read quality control |
| [fastp](https://github.com/OpenGene/fastp) | ≥ 0.23 | Adapter trimming |
| [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) | ≥ 2.2 | Read alignment |
| [samtools](http://www.htslib.org) | ≥ 1.17 | BAM manipulation |
| [GATK4](https://gatk.broadinstitute.org) | ≥ 4.4 | Variant calling |
| [SnpEff](https://pcingola.github.io/SnpEff/) | ≥ 5.2 | Variant annotation |
| [bcftools](http://www.htslib.org) | ≥ 1.17 | VCF manipulation |
| [MultiQC](https://multiqc.info) | ≥ 1.19 | QC aggregation |
| Python | ≥ 3.9 | Summary reports |

### Python packages
```bash
pip install pandas pyarrow matplotlib seaborn reportlab
```

### Easiest: use Docker (no manual installs)
```bash
# All tools pulled automatically via containers
nextflow run main.nf -profile docker
```

---

## Quick Start

### Option A — Synthetic test data (no downloads needed)
```bash
# Generate synthetic FASTQ + reference
python3 bin/generate_test_data.py data/

# Run pipeline (skip BQSR since no known_sites)
nextflow run main.nf \
    -profile test \
    --ref data/reference/ref.fa \
    --known_sites ""
```

### Option B — Real data (NA12878 Platinum Genome, chr22)
```bash
# Download data (~500MB, ~10 min)
bash setup_real_data.sh

# Run with Docker
nextflow run main.nf \
    -profile docker \
    --ref data/reference/chr22.fa \
    --known_sites data/reference/known_sites_chr22.vcf.gz \
    --snpeff_db GRCh38.99
```

### Option C — HPC cluster (SLURM + Singularity)
```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --ref /path/to/GRCh38.fa \
    --known_sites /path/to/dbsnp.vcf.gz \
    --snpeff_db GRCh38.99 \
    -resume
```

---

## Input

### Samplesheet (`data/samplesheet.csv`)
```csv
sample_id,fastq_R1,fastq_R2
SAMPLE_A,data/reads/SAMPLE_A_R1.fastq.gz,data/reads/SAMPLE_A_R2.fastq.gz
SAMPLE_B,data/reads/SAMPLE_B_R1.fastq.gz,data/reads/SAMPLE_B_R2.fastq.gz
SAMPLE_C,data/reads/SAMPLE_C_R1.fastq.gz,data/reads/SAMPLE_C_R2.fastq.gz
```

### Reference genome
- FASTA format (`.fa` or `.fa.gz`)
- Must match the SnpEff database version
- Indexes created automatically by the pipeline

### Known sites (for BQSR)
- VCF.gz + tabix index (`.vcf.gz.tbi`)
- Use dbSNP 138+ and/or Mills indels from GATK resource bundle
- Set `--known_sites ""` to skip BQSR (synthetic data)

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `data/samplesheet.csv` | Input sample CSV |
| `--ref` | `data/reference/ref.fa` | Reference FASTA |
| `--known_sites` | `data/reference/known_sites.vcf.gz` | dbSNP VCF for BQSR |
| `--intervals` | `` | BED/interval list (optional) |
| `--ploidy` | `2` | Sample ploidy |
| `--snpeff_db` | `GRCh38.99` | SnpEff database |
| `--outdir` | `results` | Output directory |
| `--snp_QD` | `2.0` | SNP filter: QualByDepth |
| `--snp_FS` | `60.0` | SNP filter: FisherStrand |
| `--snp_MQ` | `40.0` | SNP filter: MappingQuality |
| `--indel_QD` | `2.0` | INDEL filter: QualByDepth |
| `--indel_FS` | `200.0` | INDEL filter: FisherStrand |

---

## Output Structure

```
results/
├── qc/
│   ├── fastqc_raw/                  # Raw read QC (per sample)
│   ├── fastqc_trimmed/              # Post-trim QC (per sample)
│   ├── flagstat/                    # Alignment statistics
│   ├── coverage/                    # Per-chromosome coverage
│   └── multiqc/
│       └── multiqc_report.html      ★ Main QC report
│
├── trimmed/                         # fastp-trimmed reads (optional)
│
├── bam/
│   ├── sorted/                      # Coordinate-sorted BAMs
│   ├── dedup/                       # Deduplicated BAMs + metrics
│   └── bqsr/                        # BQSR-recalibrated BAMs
│
├── gvcf/                            # Per-sample gVCFs
│   └── SAMPLE_A/
│       └── SAMPLE_A.g.vcf.gz
│
├── vcf/
│   ├── combined/                    # CombineGVCFs output
│   ├── genotyped/                   # Raw joint-called VCF
│   └── filtered/
│       ├── snps_filtered.vcf.gz
│       ├── indels_filtered.vcf.gz
│       └── final_variants.vcf.gz    ★ Final merged VCF
│       └── final_variants_PASS.vcf.gz  ★ PASS-only VCF
│
├── annotation/
│   ├── annotated.vcf.gz             ★ SnpEff-annotated VCF
│   ├── snpeff_summary.html          # Annotation statistics
│   └── snpeff_genes.txt             # Per-gene summary
│
├── report/
│   ├── summary_variants.tsv         ★ Variant table (all PASS variants)
│   ├── summary_stats.json           # JSON statistics
│   └── summary_report.pdf           ★ PDF report with plots
│
└── pipeline_info/
    ├── timeline.html                # Task timeline
    ├── report.html                  # Resource usage report
    ├── trace.txt                    # Per-task resource usage
    └── dag.svg                      # Pipeline DAG diagram
```

---

## Understanding the Key Steps

### Why gVCF mode?
HaplotypeCaller is run with `-ERC GVCF` to produce a **genomic VCF** — this includes both variant sites AND reference confidence blocks. This allows:
- Distinguishing true homozygous-reference from missing coverage
- **Joint genotyping** across all samples simultaneously (CombineGVCFs → GenotypeGVCFs)
- Adding new samples later without rerunning everything

### Why BQSR?
Base Quality Score Recalibration corrects systematic errors in the quality scores assigned by the sequencer. It uses known variant sites (dbSNP) to learn the error model and produces more accurate quality scores, leading to better variant calls.

### Hard filtering vs VQSR
- **Hard filtering** (used here): simple threshold-based filters. Best for small cohorts (<30 samples) or when using a non-standard reference.
- **VQSR** (GATK Variant Quality Score Recalibration): machine-learning-based, recommended for large cohorts (>30 WGS or >100 WES samples). Requires large call sets.

### Ti/Tv ratio (quality check)
A key QC metric for SNP calls:
- **WGS**: expected Ti/Tv ≈ 2.0–2.1
- **WES**: expected Ti/Tv ≈ 3.0–3.3
- Values far outside these ranges indicate poor variant calls.

---

## Scaling to Thousands of Samples

### HPC parallelism
- `HAPLOTYPE_CALLER` runs independently per sample → N samples = N parallel SLURM jobs
- No code changes needed: just add samples to the samplesheet

### Scatter-gather for large genomes
For whole-genome sequencing, parallelize per chromosome:
```nextflow
// Split by chromosome
ch_intervals = Channel.from(["chr1","chr2",...,"chr22","chrX","chrY"])
HAPLOTYPE_CALLER(ch_final_bam.combine(ch_intervals))
```

### GenomicsDBImport (>30 samples)
Replace `CombineGVCFs` with `GenomicsDBImport` for cohorts >30 samples:
- Much faster and more memory-efficient
- Supports incremental addition of new samples

### Cloud (AWS/GCP)
```bash
nextflow run main.nf -profile aws \
    --ref s3://my-bucket/ref/GRCh38.fa \
    --samplesheet s3://my-bucket/samples.csv
```
Nextflow handles automatic job submission to AWS Batch or Google Life Sciences.

---

## Benchmarks (typical WGS 30× coverage)

| Step | Wall time (1 sample, 8 CPUs) | Memory |
|------|------------------------------|--------|
| FastQC | ~5 min | 2 GB |
| fastp | ~10 min | 4 GB |
| BWA-MEM2 | ~45 min | 12 GB |
| MarkDuplicates | ~20 min | 8 GB |
| BQSR | ~30 min | 8 GB |
| HaplotypeCaller | ~4 h | 16 GB |
| GenotypeGVCFs | ~1 h | 16 GB |
| SnpEff | ~10 min | 8 GB |
| **Total** | **~7 h** | **16 GB peak** |

---

## Troubleshooting

**"Unable to determine SAM format"**  
→ Ensure BAM has proper `@RG` read group tags (added by `BWA_MEM` process).

**GATK out of memory**  
→ Increase `memory` in `nextflow.config` for the failing process label.

**BQSR fails: "No overlapping features"**  
→ Ensure your known_sites VCF chromosome names match your reference (e.g., both use `chr22` or both use `22`).

**SnpEff database not found**  
→ Download manually: `snpEff download GRCh38.99` or set `--snpeff_db` to your genome build.

**Resume a failed run**  
```bash
nextflow run main.nf -resume -profile docker
```

---

*Pipeline developed for the Bioinformatician position — Dr. Jacquemont's Lab, CHU Sainte-Justine*
