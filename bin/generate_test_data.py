#!/usr/bin/env python3
"""
Generate synthetic paired-end FASTQ test data for pipeline validation.
Simulates 2 samples with realistic read characteristics.
"""

import random
import gzip
import os
import sys

random.seed(42)

# ─── Reference-like sequence (simulated chr1 mini-region, 10kb) ─────────────
BASES = "ACGT"
WEIGHTS = [0.30, 0.20, 0.20, 0.30]  # slight AT bias like human genome

def random_sequence(length):
    return "".join(random.choices(BASES, weights=WEIGHTS, k=length))

def complement(base):
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}.get(base, "N")

def reverse_complement(seq):
    return "".join(complement(b) for b in reversed(seq))

def introduce_snp(seq, rate=0.001):
    """Introduce random SNPs to simulate a variant sample."""
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < rate:
            alt = [b for b in BASES if b != seq[i]]
            seq[i] = random.choice(alt)
    return "".join(seq)

def phred_quality(length, mean_q=35, end_drop=True):
    """Simulate realistic Phred quality scores (Illumina-like)."""
    quals = []
    for i in range(length):
        if end_drop and i > length * 0.85:
            q = max(10, mean_q - int((i - length * 0.85) * 2))
        else:
            q = random.randint(mean_q - 5, mean_q + 3)
        quals.append(chr(min(q + 33, 126)))
    return "".join(quals)

def simulate_reads(ref, n_reads=5000, read_length=150, insert_size=400,
                   snp_rate=0.001, sample_name="sample"):
    """
    Simulate paired-end reads from a reference sequence.
    Returns (R1_records, R2_records)
    """
    r1_records = []
    r2_records = []
    ref_len = len(ref)

    for i in range(n_reads):
        # Random start position
        max_start = ref_len - insert_size - read_length
        if max_start <= 0:
            continue
        start = random.randint(0, max_start)

        # Extract insert + introduce variants
        insert = ref[start: start + insert_size + read_length]
        insert = introduce_snp(insert, rate=snp_rate)

        # R1: forward strand
        r1_seq = insert[:read_length]
        # R2: reverse complement of end of insert
        r2_seq = reverse_complement(insert[insert_size: insert_size + read_length])

        # Quality scores
        r1_qual = phred_quality(read_length)
        r2_qual = phred_quality(read_length)

        # Read names
        read_id = f"@{sample_name}.{i+1} 1:N:0:ATCACG"

        r1_records.append(f"@{sample_name}.{i+1} 1:N:0:ATCACG\n{r1_seq}\n+\n{r1_qual}")
        r2_records.append(f"@{sample_name}.{i+1} 2:N:0:ATCACG\n{r2_seq}\n+\n{r2_qual}")

    return r1_records, r2_records

def write_fastq_gz(records, filepath):
    with gzip.open(filepath, "wt") as f:
        f.write("\n".join(records) + "\n")
    print(f"  ✓ Written: {filepath} ({len(records)} reads)")

def write_reference(ref, filepath):
    """Write reference as FASTA."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(">chr1_synthetic\n")
        # Write 60 chars per line
        for i in range(0, len(ref), 60):
            f.write(ref[i:i+60] + "\n")
    print(f"  ✓ Reference written: {filepath} ({len(ref)} bp)")

def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(out_dir, exist_ok=True)
    reads_dir = os.path.join(out_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)

    print("🧬 Generating synthetic reference genome (50kb)...")
    ref = random_sequence(50_000)
    write_reference(ref, os.path.join(out_dir, "reference", "ref.fa"))

    # Sample metadata
    samples = [
        {"id": "SAMPLE_A", "snp_rate": 0.002, "n_reads": 8000},
        {"id": "SAMPLE_B", "snp_rate": 0.0015, "n_reads": 8000},
        {"id": "SAMPLE_C", "snp_rate": 0.003, "n_reads": 8000},
    ]

    print(f"\n📄 Generating reads for {len(samples)} samples...")
    csv_lines = ["sample_id,fastq_R1,fastq_R2"]

    for s in samples:
        sid = s["id"]
        print(f"\n  Sample: {sid}")
        r1_records, r2_records = simulate_reads(
            ref,
            n_reads=s["n_reads"],
            read_length=150,
            insert_size=350,
            snp_rate=s["snp_rate"],
            sample_name=sid
        )
        r1_path = os.path.join(reads_dir, f"{sid}_R1.fastq.gz")
        r2_path = os.path.join(reads_dir, f"{sid}_R2.fastq.gz")
        write_fastq_gz(r1_records, r1_path)
        write_fastq_gz(r2_records, r2_path)
        csv_lines.append(f"{sid},{r1_path},{r2_path}")

    # Write samplesheet
    samplesheet = os.path.join(out_dir, "samplesheet.csv")
    with open(samplesheet, "w") as f:
        f.write("\n".join(csv_lines) + "\n")
    print(f"\n  ✓ Samplesheet: {samplesheet}")
    print("\n✅ Test data generation complete!")

if __name__ == "__main__":
    main()
