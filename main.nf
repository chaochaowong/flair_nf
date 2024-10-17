#!/usr/bin/env nextflow

log.info """\
  FLAIR_NF PIPELINE on mamba and slurm executor 
  ================================================
  Sample sheet      : ${params.sample_sheet}
  Project directory : ${projectDir}
  Work directory    : ${workDir}
  Genomic reference : ${params.genome_reference}
  GTF               : ${params.gtf}
  Output directory  : ${params.outdir}
"""
.stripIndent()

// Include modules
include { FLAIR_ALIGN } from './modules/local/flair/align/main.nf'
include { FLAIR_CORRECT } from './modules/local/flair/correct/main.nf'
include { CONCAT_FASTQ } from './modules/local/concat_fastq/main.nf'
include { CONCAT_CORRECTED_BED } from './modules/local/concat_corrected_bed/main.nf'
include { FLAIR_COLLAPSE } from './modules/local/flair/collapse/main.nf'
include { SAMPLE_MANIFEST_TSV }  from './modules/local/sample_manifest_tsv/main.nf'
include { FLAIR_QUANTIFY } from './modules/local/flair/quantify/main.nf'

workflow {
    // Create input channel from samplesheet in TSV format
    reads_ch = Channel.fromPath( params.sample_sheet )
                      .splitCsv( header: true )
                      .map { row -> 
                        [row.sample_id, file(row.reads, checkIfExists: true)] }

    // flair align                      
    FLAIR_ALIGN(reads_ch, 
                params.genome_reference, 
                params.genome_reference_index, 
                params.gtf, 
                params.min_mapq)

    // flair correct
    FLAIR_CORRECT(FLAIR_ALIGN.out, 
                  params.genome_reference, 
                  params.gtf)

    // concatnate all fastq 
    Channel.fromPath( params.sample_sheet )
      .splitCsv( header: true )
      .map { row -> file(row.reads, checkIfExists: true) } 
      .collect()
      | CONCAT_FASTQ

    // concatenate corrected bed files
    CONCAT_CORRECTED_BED(FLAIR_CORRECT.out.collect())

    // flair collapse 
    FLAIR_COLLAPSE(params.genome_reference, 
                   params.gtf, 
                   CONCAT_FASTQ.out, 
                   CONCAT_CORRECTED_BED.out)

    // sample manifest file
    SAMPLE_MANIFEST_TSV(params.sample_sheet)

    // flair quant
    FLAIR_QUANTIFY(SAMPLE_MANIFEST_TSV.out, 
                   FLAIR_COLLAPSE.out )
    // FLAIR_QUANTIFY.out.view { it }

}