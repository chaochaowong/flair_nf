/*
* flair align: aligns fastq files to using minmap2 and converts the bam file
* to a bed file
* NOTE: flair align should allow hifi reads and use PacBio's wrapper function of minmap2
*/
process FLAIR_ALIGN {
    tag "FLAIR_ALIGN on ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input: 
        tuple val(sample_id), val(condition), val(batch), path(fastq)
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
        --threads $task.cpus
    """        
}
