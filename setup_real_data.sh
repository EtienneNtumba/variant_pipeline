#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
#  setup_real_data.sh
#  Downloads a small real-world dataset from 1000 Genomes Project
#  + the required GATK resource bundle files
#
#  Requirements: curl, samtools, bwa-mem2, gatk (in PATH or modules loaded)
#  Time: ~5-10 minutes on a good connection
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Variant Pipeline — Real Data Setup                      ║"
echo "╚══════════════════════════════════════════════════════════╝"

# ── Directories ─────────────────────────────────────────────────────────────
DATA_DIR="data"
REF_DIR="${DATA_DIR}/reference"
READS_DIR="${DATA_DIR}/reads"
mkdir -p "${REF_DIR}" "${READS_DIR}"

# ═══════════════════════════════════════════════════════════════════════════
# 1. Download reference — chr22 from GRCh38 (small, ~50MB)
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "📥 [1/5] Downloading chr22 reference (GRCh38)..."
REF="${REF_DIR}/chr22.fa"

if [ ! -f "${REF}" ]; then
    curl -L -o "${REF}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr22.fna.gz"
    echo "  Decompressing..."
    gunzip -c "${REF}.gz" | sed 's/>CM000684.2/>chr22/' > "${REF}"
    rm "${REF}.gz"
    echo "  ✓ Reference: ${REF}"
else
    echo "  ✓ Reference already exists, skipping."
fi

# ── Index reference ──────────────────────────────────────────────────────────
echo "  Indexing reference with bwa-mem2..."
bwa-mem2 index "${REF}"
samtools faidx  "${REF}"
gatk CreateSequenceDictionary -R "${REF}"
echo "  ✓ Reference indexed."

# ═══════════════════════════════════════════════════════════════════════════
# 2. Download GATK resource bundle — known sites for BQSR (chr22 subset)
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "📥 [2/5] Downloading dbSNP known sites (chr22, GRCh38)..."
KNOWN_SITES="${REF_DIR}/dbsnp_chr22.vcf.gz"

if [ ! -f "${KNOWN_SITES}" ]; then
    curl -L -o "${KNOWN_SITES}" \
        "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    # Subset to chr22 to keep it small
    echo "  Subsetting to chr22..."
    bcftools view -r chr22 "${KNOWN_SITES}" -Oz -o "${REF_DIR}/known_sites_chr22.vcf.gz"
    tabix -p vcf "${REF_DIR}/known_sites_chr22.vcf.gz"
    rm "${KNOWN_SITES}"
    KNOWN_SITES="${REF_DIR}/known_sites_chr22.vcf.gz"
    echo "  ✓ Known sites: ${KNOWN_SITES}"
else
    echo "  ✓ Known sites already exist, skipping."
fi

# ═══════════════════════════════════════════════════════════════════════════
# 3. Download real WGS reads from 1000 Genomes (chr22 region only)
#    NA12878 — HG001, well-characterized platinum genome
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "📥 [3/5] Downloading chr22 reads for NA12878 (HG001)..."

# Use a public BAM from 1000 Genomes and extract a small region
REMOTE_BAM="https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12878_S1.bam"
REGION="chr22:20000000-21000000"  # 1Mb region of chr22

SAMPLE="NA12878"
OUT_BAM="${READS_DIR}/${SAMPLE}_chr22_region.bam"

if [ ! -f "${READS_DIR}/${SAMPLE}_R1.fastq.gz" ]; then
    echo "  Extracting reads from region ${REGION}..."
    samtools view -b -h "${REMOTE_BAM}" "${REGION}" \
        -o "${OUT_BAM}" --reference "${REF}"

    echo "  Converting BAM → FASTQ (paired-end)..."
    samtools sort -n "${OUT_BAM}" -o "${OUT_BAM%.bam}.namesorted.bam"
    samtools fastq \
        -1 "${READS_DIR}/${SAMPLE}_R1.fastq.gz" \
        -2 "${READS_DIR}/${SAMPLE}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n \
        "${OUT_BAM%.bam}.namesorted.bam"
    rm "${OUT_BAM}" "${OUT_BAM%.bam}.namesorted.bam"
    echo "  ✓ Reads: ${READS_DIR}/${SAMPLE}_R1/R2.fastq.gz"
else
    echo "  ✓ Reads already exist, skipping."
fi

# Add HG002 (second sample for joint calling demo)
echo ""
echo "📥 [4/5] Downloading chr22 reads for NA12877 (father of NA12878)..."
SAMPLE2="NA12877"
REMOTE_BAM2="https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12877_S1.bam"

if [ ! -f "${READS_DIR}/${SAMPLE2}_R1.fastq.gz" ]; then
    samtools view -b -h "${REMOTE_BAM2}" "${REGION}" \
        -o "${READS_DIR}/${SAMPLE2}_chr22.bam" --reference "${REF}"
    samtools sort -n "${READS_DIR}/${SAMPLE2}_chr22.bam" \
        -o "${READS_DIR}/${SAMPLE2}_ns.bam"
    samtools fastq \
        -1 "${READS_DIR}/${SAMPLE2}_R1.fastq.gz" \
        -2 "${READS_DIR}/${SAMPLE2}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n \
        "${READS_DIR}/${SAMPLE2}_ns.bam"
    rm "${READS_DIR}/${SAMPLE2}_chr22.bam" "${READS_DIR}/${SAMPLE2}_ns.bam"
    echo "  ✓ Reads: ${READS_DIR}/${SAMPLE2}_R1/R2.fastq.gz"
else
    echo "  ✓ Reads already exist, skipping."
fi

# ═══════════════════════════════════════════════════════════════════════════
# 4. Write samplesheet
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "📝 [5/5] Writing samplesheet..."
cat > "${DATA_DIR}/samplesheet.csv" <<EOF
sample_id,fastq_R1,fastq_R2
NA12878,${READS_DIR}/NA12878_R1.fastq.gz,${READS_DIR}/NA12878_R2.fastq.gz
NA12877,${READS_DIR}/NA12877_R1.fastq.gz,${READS_DIR}/NA12877_R2.fastq.gz
EOF
echo "  ✓ Samplesheet: ${DATA_DIR}/samplesheet.csv"

# ═══════════════════════════════════════════════════════════════════════════
# 5. Summary
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Setup complete! Run the pipeline with:                  ║"
echo "║                                                          ║"
echo "║  # With Docker (recommended):                            ║"
echo "║  nextflow run main.nf -profile docker \\                  ║"
echo "║    --ref data/reference/chr22.fa \\                       ║"
echo "║    --known_sites data/reference/known_sites_chr22.vcf.gz ║"
echo "║                                                          ║"
echo "║  # On SLURM cluster:                                     ║"
echo "║  nextflow run main.nf -profile slurm,singularity \\       ║"
echo "║    --ref data/reference/chr22.fa \\                       ║"
echo "║    --known_sites data/reference/known_sites_chr22.vcf.gz ║"
echo "║                                                          ║"
echo "║  # Resume after failure:                                 ║"
echo "║  nextflow run main.nf -profile docker -resume            ║"
echo "╚══════════════════════════════════════════════════════════╝"
