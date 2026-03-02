#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
╔══════════════════════════════════════════════════════════════════════════════╗
║          GERMLINE VARIANT CALLING PIPELINE — DSL2                           ║
║          FastQC → Trim → Align → BQSR → HaplotypeCaller → Annotate         ║
║          Based on GATK4 Best Practices                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

WORKFLOW:
  1.  FASTQC_RAW          – QC on raw reads
  2.  FASTP               – Adapter trimming + QC
  3.  FASTQC_TRIMMED      – QC post-trim
  4.  BWA_INDEX           – Index reference genome (once)
  5.  BWA_MEM             – Align paired-end reads
  6.  SAMTOOLS_SORT       – Sort BAM by coordinate
  7.  SAMTOOLS_INDEX      – Index sorted BAM
  8.  MARK_DUPLICATES     – Flag PCR duplicates (GATK/Picard)
  9.  FLAGSTAT            – Alignment statistics
 10.  BASE_RECALIBRATOR   – BQSR: compute recalibration table
 11.  APPLY_BQSR          – BQSR: apply recalibration
 12.  HAPLOTYPE_CALLER    – Call variants in gVCF mode (per sample)
 13.  COMBINE_GVCFS       – Combine per-sample gVCFs
 14.  GENOTYPE_GVCFS      – Joint genotyping → raw VCF
 15.  SELECT_SNP          – Extract SNPs
 16.  SELECT_INDEL        – Extract INDELs
 17.  FILTER_SNP          – Hard-filter SNPs
 18.  FILTER_INDEL        – Hard-filter INDELs
 19.  MERGE_VCFS          – Merge filtered SNPs + INDELs
 20.  SNPEFF              – Functional annotation
 21.  SUMMARIZE           – TSV + PDF summary report
 22.  MULTIQC             – Aggregate all QC reports
*/

// ─── PARAMETERS ────────────────────────────────────────────────────────────
params {
    // Input
    samplesheet   = "data/samplesheet.csv"
    ref           = "data/reference/ref.fa"

    // Known variant sites for BQSR (required for real data; use empty for synthetic)
    known_sites   = "data/reference/known_sites.vcf.gz"  // set to "" to skip BQSR

    // GATK HaplotypeCaller settings
    intervals     = ""          // optional BED/interval list
    ploidy        = 2

    // Filtering thresholds — SNPs
    snp_QD        = 2.0
    snp_FS        = 60.0
    snp_MQ        = 40.0
    snp_MQRankSum = -12.5
    snp_ReadPosRankSum = -8.0

    // Filtering thresholds — INDELs
    indel_QD      = 2.0
    indel_FS      = 200.0
    indel_ReadPosRankSum = -20.0

    // SnpEff genome database (must match your reference; use 'GRCh38.99' for real data)
    snpeff_db     = "GRCh38.99"

    // Output
    outdir        = "results"

    // Resources
    max_cpus      = 8
    max_memory    = "32.GB"
    max_time      = "24.h"
}

// ─── PROCESSES ─────────────────────────────────────────────────────────────

// ── 1. FastQC on raw reads ──────────────────────────────────────────────────
process FASTQC_RAW {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/qc/fastqc_raw/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*.html"), path("*.zip"), emit: reports
    path "*.zip",                                         emit: zip

    script:
    """
    fastqc --threads ${task.cpus} \\
           --outdir . \\
           ${r1} ${r2}
    """
}

// ── 2. Adapter trimming with fastp ─────────────────────────────────────────
process FASTP {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: "copy",
        saveAs: { fn -> fn.endsWith(".html") || fn.endsWith(".json") ? fn : null }

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1_trimmed.fastq.gz"),
          path("${sample_id}_R2_trimmed.fastq.gz"),  emit: reads
    path "${sample_id}_fastp.html",                   emit: html
    path "${sample_id}_fastp.json",                   emit: json

    script:
    """
    fastp \\
        --in1  ${r1} \\
        --in2  ${r2} \\
        --out1 ${sample_id}_R1_trimmed.fastq.gz \\
        --out2 ${sample_id}_R2_trimmed.fastq.gz \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --correction \\
        --cut_tail \\
        --cut_tail_window_size 4 \\
        --cut_tail_mean_quality 20 \\
        --length_required 36 \\
        --qualified_quality_phred 20
    """
}

// ── 3. FastQC on trimmed reads ─────────────────────────────────────────────
process FASTQC_TRIMMED {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/qc/fastqc_trimmed/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
}

// ── 4. Index reference genome ──────────────────────────────────────────────
process BWA_INDEX {
    label "process_medium"
    publishDir "${params.outdir}/reference", mode: "copy"

    input:
    path ref

    output:
    tuple path(ref), path("${ref}.*"), emit: indexed_ref

    script:
    """
    # BWA-MEM2 index
    bwa-mem2 index ${ref}

    # Samtools FASTA index (for GATK)
    samtools faidx ${ref}

    # GATK sequence dictionary
    gatk CreateSequenceDictionary -R ${ref}
    """
}

// ── 5. Align reads with BWA-MEM2 ──────────────────────────────────────────
process BWA_MEM {
    tag "$sample_id"
    label "process_high"
    publishDir "${params.outdir}/bam/raw/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(r1), path(r2)
    tuple path(ref), path(ref_idx)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    // Read group tags — essential for GATK
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${ref} ${r1} ${r2} \\
    | samtools view -bS - > ${sample_id}.bam
    """
}

// ── 6. Sort BAM ────────────────────────────────────────────────────────────
process SAMTOOLS_SORT {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/bam/sorted/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o ${sample_id}.sorted.bam \\
        ${bam}
    """
}

// ── 7. Index sorted BAM ────────────────────────────────────────────────────
process SAMTOOLS_INDEX {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/bam/sorted/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_bai

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

// ── 8. Mark duplicates (GATK MarkDuplicates / Picard) ─────────────────────
process MARK_DUPLICATES {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/bam/dedup/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.dedup.bam"),
          path("${sample_id}.dedup.bam.bai"), emit: bam_bai
    path "${sample_id}.duplicate_metrics.txt", emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        --INPUT  ${bam} \\
        --OUTPUT ${sample_id}.dedup.bam \\
        --METRICS_FILE ${sample_id}.duplicate_metrics.txt \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
        --TMP_DIR ./tmp \\
        --java-options "-Xmx${task.memory.toGiga()}g"

    # Rename index to .bai convention
    mv ${sample_id}.dedup.bai ${sample_id}.dedup.bam.bai || true
    samtools index ${sample_id}.dedup.bam 2>/dev/null || true
    """
}

// ── 9. Alignment statistics ────────────────────────────────────────────────
process FLAGSTAT {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/qc/flagstat", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.flagstat.txt", emit: stats

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${sample_id}.flagstat.txt
    """
}

// ── 10. Coverage statistics ────────────────────────────────────────────────
process COVERAGE_STATS {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/qc/coverage", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.coverage.txt", emit: stats

    script:
    """
    samtools coverage ${bam} > ${sample_id}.coverage.txt
    """
}

// ── 11. Base Quality Score Recalibration — Step 1 ─────────────────────────
process BASE_RECALIBRATOR {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/bam/bqsr/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_idx
    path known_sites      // .vcf.gz
    path known_sites_tbi  // .vcf.gz.tbi

    output:
    tuple val(sample_id), path("${sample_id}.recal.table"), emit: table

    script:
    def intervals_arg = params.intervals ? "--intervals ${params.intervals}" : ""
    """
    gatk BaseRecalibrator \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -I  ${bam} \\
        -R  ${ref} \\
        --known-sites ${known_sites} \\
        -O  ${sample_id}.recal.table \\
        ${intervals_arg}
    """
}

// ── 12. Apply BQSR ────────────────────────────────────────────────────────
process APPLY_BQSR {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/bam/bqsr/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple val(sample_id2), path(recal_table)
    path ref
    path ref_idx

    output:
    tuple val(sample_id),
          path("${sample_id}.bqsr.bam"),
          path("${sample_id}.bqsr.bam.bai"), emit: bam_bai

    script:
    def intervals_arg = params.intervals ? "--intervals ${params.intervals}" : ""
    """
    gatk ApplyBQSR \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -I  ${bam} \\
        -R  ${ref} \\
        --bqsr-recal-file ${recal_table} \\
        -O  ${sample_id}.bqsr.bam \\
        ${intervals_arg}

    samtools index ${sample_id}.bqsr.bam
    """
}

// ── 13. HaplotypeCaller — gVCF mode (per sample) ─────────────────────────
process HAPLOTYPE_CALLER {
    tag "$sample_id"
    label "process_high"
    publishDir "${params.outdir}/gvcf/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_idx

    output:
    tuple val(sample_id),
          path("${sample_id}.g.vcf.gz"),
          path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def intervals_arg = params.intervals ? "--intervals ${params.intervals}" : ""
    """
    gatk HaplotypeCaller \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -I   ${bam} \\
        -R   ${ref} \\
        -O   ${sample_id}.g.vcf.gz \\
        -ERC GVCF \\
        --sample-ploidy ${params.ploidy} \\
        --native-pair-hmm-threads ${task.cpus} \\
        ${intervals_arg}
    """
}

// ── 14. CombineGVCFs (small cohorts) or GenomicsDBImport (large) ──────────
process COMBINE_GVCFS {
    label "process_medium"
    publishDir "${params.outdir}/vcf/combined", mode: "copy"

    input:
    path gvcfs    // list of .g.vcf.gz
    path tbis     // list of .g.vcf.gz.tbi
    path ref
    path ref_idx

    output:
    tuple path("combined.g.vcf.gz"), path("combined.g.vcf.gz.tbi"), emit: combined

    script:
    def v_args = gvcfs.collect { "-V $it" }.join(" \\\n        ")
    """
    gatk CombineGVCFs \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} \\
        ${v_args} \\
        -O combined.g.vcf.gz
    """
}

// ── 15. Joint genotyping ──────────────────────────────────────────────────
process GENOTYPE_GVCFS {
    label "process_high"
    publishDir "${params.outdir}/vcf/genotyped", mode: "copy"

    input:
    tuple path(combined_gvcf), path(combined_tbi)
    path ref
    path ref_idx

    output:
    tuple path("genotyped.vcf.gz"), path("genotyped.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk GenotypeGVCFs \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} \\
        -V ${combined_gvcf} \\
        -O genotyped.vcf.gz \\
        --sample-ploidy ${params.ploidy}
    """
}

// ── 16. Select SNPs ────────────────────────────────────────────────────────
process SELECT_SNP {
    label "process_low"
    publishDir "${params.outdir}/vcf/filtered", mode: "copy"

    input:
    tuple path(vcf), path(tbi)
    path ref
    path ref_idx

    output:
    tuple path("snps_raw.vcf.gz"), path("snps_raw.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk SelectVariants \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} -V ${vcf} \\
        --select-type-to-include SNP \\
        -O snps_raw.vcf.gz
    """
}

// ── 17. Select INDELs ─────────────────────────────────────────────────────
process SELECT_INDEL {
    label "process_low"
    publishDir "${params.outdir}/vcf/filtered", mode: "copy"

    input:
    tuple path(vcf), path(tbi)
    path ref
    path ref_idx

    output:
    tuple path("indels_raw.vcf.gz"), path("indels_raw.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk SelectVariants \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O indels_raw.vcf.gz
    """
}

// ── 18. Hard-filter SNPs ──────────────────────────────────────────────────
process FILTER_SNP {
    label "process_low"
    publishDir "${params.outdir}/vcf/filtered", mode: "copy"

    input:
    tuple path(vcf), path(tbi)
    path ref
    path ref_idx

    output:
    tuple path("snps_filtered.vcf.gz"), path("snps_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk VariantFiltration \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} -V ${vcf} \\
        --filter-expression "QD < ${params.snp_QD}"         --filter-name "QD2"          \\
        --filter-expression "FS > ${params.snp_FS}"          --filter-name "FS60"         \\
        --filter-expression "MQ < ${params.snp_MQ}"          --filter-name "MQ40"         \\
        --filter-expression "MQRankSum < ${params.snp_MQRankSum}" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < ${params.snp_ReadPosRankSum}" --filter-name "ReadPosRankSum-8" \\
        -O snps_filtered.vcf.gz
    """
}

// ── 19. Hard-filter INDELs ────────────────────────────────────────────────
process FILTER_INDEL {
    label "process_low"
    publishDir "${params.outdir}/vcf/filtered", mode: "copy"

    input:
    tuple path(vcf), path(tbi)
    path ref
    path ref_idx

    output:
    tuple path("indels_filtered.vcf.gz"), path("indels_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk VariantFiltration \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -R ${ref} -V ${vcf} \\
        --filter-expression "QD < ${params.indel_QD}"                        --filter-name "QD2"      \\
        --filter-expression "FS > ${params.indel_FS}"                         --filter-name "FS200"    \\
        --filter-expression "ReadPosRankSum < ${params.indel_ReadPosRankSum}" --filter-name "ReadPosRankSum-20" \\
        -O indels_filtered.vcf.gz
    """
}

// ── 20. Merge filtered SNPs + INDELs ─────────────────────────────────────
process MERGE_VCFS {
    label "process_low"
    publishDir "${params.outdir}/vcf/final", mode: "copy"

    input:
    tuple path(snp_vcf),   path(snp_tbi)
    tuple path(indel_vcf), path(indel_tbi)
    path ref
    path ref_idx

    output:
    tuple path("final_variants.vcf.gz"), path("final_variants.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk MergeVcfs \\
        --java-options "-Xmx${task.memory.toGiga()}g" \\
        -I ${snp_vcf} \\
        -I ${indel_vcf} \\
        -O final_variants.vcf.gz

    # PASS-only VCF (useful for downstream)
    bcftools view -f PASS final_variants.vcf.gz -Oz -o final_variants_PASS.vcf.gz
    bcftools index -t final_variants_PASS.vcf.gz
    """
}

// ── 21. Functional annotation with SnpEff ────────────────────────────────
process SNPEFF {
    label "process_medium"
    publishDir "${params.outdir}/annotation", mode: "copy"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("annotated.vcf.gz"), path("annotated.vcf.gz.tbi"), emit: vcf
    path "snpeff_summary.html",                                    emit: html
    path "snpeff_genes.txt",                                       emit: genes

    script:
    """
    snpEff ann \\
        -v \\
        -stats snpeff_summary.html \\
        -csvStats snpeff_genes.txt \\
        ${params.snpeff_db} \\
        ${vcf} \\
    | bgzip > annotated.vcf.gz

    tabix -p vcf annotated.vcf.gz
    """
}

// ── 22. Summary report (TSV + PDF) ───────────────────────────────────────
process SUMMARIZE {
    label "process_low"
    publishDir "${params.outdir}/report", mode: "copy"

    input:
    tuple path(vcf), path(tbi)

    output:
    path "summary_variants.tsv", emit: tsv
    path "summary_stats.json",   emit: json
    path "summary_report.pdf",   emit: pdf, optional: true

    script:
    """
    python3 ${projectDir}/bin/summarize_variants.py ${vcf} summary
    """
}

// ── 23. MultiQC — aggregate all QC reports ───────────────────────────────
process MULTIQC {
    label "process_low"
    publishDir "${params.outdir}/qc/multiqc", mode: "copy"

    input:
    path "*"    // all QC files collected

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . \\
        --title "Variant Calling Pipeline QC" \\
        --comment "FastQC · fastp · samtools flagstat · GATK MarkDuplicates" \\
        --force \\
        -o .
    """
}

// ─── WORKFLOW ──────────────────────────────────────────────────────────────
workflow {

    // ── Load samplesheet ──────────────────────────────────────────────────
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq_R1), file(row.fastq_R2)) }

    // ── Reference ────────────────────────────────────────────────────────
    ch_ref = file(params.ref)

    // ── 1-2. QC + Trimming ───────────────────────────────────────────────
    FASTQC_RAW(ch_samples)
    FASTP(ch_samples)
    FASTQC_TRIMMED(FASTP.out.reads)

    // ── 3. Index reference ────────────────────────────────────────────────
    BWA_INDEX(ch_ref)

    // ── 4. Alignment ──────────────────────────────────────────────────────
    BWA_MEM(FASTP.out.reads, BWA_INDEX.out.indexed_ref)

    // ── 5-6. Sort + Index ─────────────────────────────────────────────────
    SAMTOOLS_SORT(BWA_MEM.out.bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    // ── 7. Mark duplicates ────────────────────────────────────────────────
    MARK_DUPLICATES(SAMTOOLS_INDEX.out.bam_bai)

    // ── 8. QC stats ───────────────────────────────────────────────────────
    FLAGSTAT(MARK_DUPLICATES.out.bam_bai)
    COVERAGE_STATS(MARK_DUPLICATES.out.bam_bai)

    // ── 9-10. BQSR (skip if no known_sites provided) ─────────────────────
    if (params.known_sites && file(params.known_sites).exists()) {
        ch_ks     = file(params.known_sites)
        ch_ks_tbi = file("${params.known_sites}.tbi")

        BASE_RECALIBRATOR(
            MARK_DUPLICATES.out.bam_bai,
            ch_ref,
            BWA_INDEX.out.indexed_ref,
            ch_ks, ch_ks_tbi
        )
        APPLY_BQSR(
            MARK_DUPLICATES.out.bam_bai,
            BASE_RECALIBRATOR.out.table,
            ch_ref,
            BWA_INDEX.out.indexed_ref
        )
        ch_final_bam = APPLY_BQSR.out.bam_bai
    } else {
        log.warn "⚠ No known_sites VCF found — skipping BQSR (OK for synthetic data)"
        ch_final_bam = MARK_DUPLICATES.out.bam_bai
    }

    // ── 11. HaplotypeCaller (per sample → gVCF) ───────────────────────────
    HAPLOTYPE_CALLER(
        ch_final_bam,
        ch_ref,
        BWA_INDEX.out.indexed_ref
    )

    // ── 12. Combine gVCFs + joint genotyping ─────────────────────────────
    COMBINE_GVCFS(
        HAPLOTYPE_CALLER.out.gvcf.map { it[1] }.collect(),
        HAPLOTYPE_CALLER.out.gvcf.map { it[2] }.collect(),
        ch_ref,
        BWA_INDEX.out.indexed_ref
    )

    GENOTYPE_GVCFS(
        COMBINE_GVCFS.out.combined,
        ch_ref,
        BWA_INDEX.out.indexed_ref
    )

    // ── 13. Split + hard-filter ───────────────────────────────────────────
    SELECT_SNP(GENOTYPE_GVCFS.out.vcf, ch_ref, BWA_INDEX.out.indexed_ref)
    SELECT_INDEL(GENOTYPE_GVCFS.out.vcf, ch_ref, BWA_INDEX.out.indexed_ref)
    FILTER_SNP(SELECT_SNP.out.vcf, ch_ref, BWA_INDEX.out.indexed_ref)
    FILTER_INDEL(SELECT_INDEL.out.vcf, ch_ref, BWA_INDEX.out.indexed_ref)

    // ── 14. Merge + annotate ──────────────────────────────────────────────
    MERGE_VCFS(
        FILTER_SNP.out.vcf,
        FILTER_INDEL.out.vcf,
        ch_ref,
        BWA_INDEX.out.indexed_ref
    )

    SNPEFF(MERGE_VCFS.out.vcf)

    // ── 15. Summary report ────────────────────────────────────────────────
    SUMMARIZE(SNPEFF.out.vcf)

    // ── 16. MultiQC ───────────────────────────────────────────────────────
    ch_multiqc = Channel.empty()
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTP.out.json)
        .mix(MARK_DUPLICATES.out.metrics)
        .mix(FLAGSTAT.out.stats)
        .collect()

    MULTIQC(ch_multiqc)
}

// ─── COMPLETION HANDLER ────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "✅ SUCCESS" : "❌ FAILED"
    def duration = workflow.duration
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║  Pipeline ${status}
    ║  Duration : ${duration}
    ║  Results  : ${params.outdir}/
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
}
