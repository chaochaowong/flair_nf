#!/usr/bin/env nextflow

log.info """\
  FLAIR pipeline on mamba and slurm executor 
  ================================================
  sample sheet: ${params.sample_sheet}
"""

process FLAIR_ALIGN {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input: 
        tuple val(sample_id), val(sample), path(fastq)
        path ref_fasta
        path ref_index
        path gtf
        val(min_mapq)
    
    output:
        tuple val(sample_id), path("${sample_id}.flair.aligned.bam"), path("${sample_id}.flair.aligned.bam.bai"), path("${sample_id}.flair.aligned.bed")

    script:
    """
    flair align -g ${ref_fasta} \
        -r ${fastq} \
        -o ${sample_id}.flair.aligned \
        --quality ${min_mapq} \
        --threads 8
    """        
}


process FLAIR_CORRECT {
    publishDir "${params.outdir}/correct", mode: 'copy'

    input: 
        tuple val(sample_id), path(sample_bam), path(sample_bai), path(sample_bed)
        path ref_fasta
        path gtf

    output: 
        tuple val(sample_id), path("${sample_id}.flair_all_corrected.bed"), path("${sample_id}.flair_all_inconsistent.bed")

    script:
    """
    flair correct -g ${ref_fasta} \
        --gtf ${gtf} \
        -q ${sample_bed} \
        -o ${sample_id}.flair \
        --threads 8
    """
}

process FLAIR_COLLAPSE_FASTQ {
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input: 
        path(sample_sheet)

    output:    
        path('combined_samples.fastq')

    script:
    """
        awk -F, 'NR>1 {print \$3}' ${sample_sheet} | xargs cat > combined_samples.fastq
    """    
}

workflow {
    // Create input channel from samplesheet in CSV format
    reads_ch = Channel.fromPath(params.sample_sheet)
                      .splitCsv(header: true)
                      .map { row -> [row.sample_id, row.sample, file(row.fq)] }
    // reads_ch.view{ it }            
    FLAIR_ALIGN(reads_ch, params.genome_reference, params.genome_reference_index, params.gtf, params.min_mapq)
    // FLAIR_ALIGN.out.view{ it }
    FLAIR_CORRECT(FLAIR_ALIGN.out, params.genome_reference, params.gtf)
    // FLAIR_CORRECT.out.view{ it }
    FLAIR_COLLAPSE_FASTQ(params.sample_sheet)
    // FLAIR_COLLAPSE_BED()

    // sample manifest file
    // FLAIR_QUANTIFY()

    
}