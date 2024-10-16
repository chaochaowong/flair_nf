/*
* flair correctr: corrects alignments to the annotated splice sites
*/ 
process FLAIR_CORRECT {
    tag "${sample_id}"
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
        --threads $task.cpus
    """
}