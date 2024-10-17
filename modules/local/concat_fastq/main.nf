/*
* Concatenate all fastq files in the 'fq' colum and deposit
# to ./collapse/combined_sample.fastq
*/

process CONCAT_FASTQ {
    tag './collapse/combined_samples.fastq'
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input: 
        path(collected_files)

    output: 
        path('combined_samples.fastq')

    script:
    """
    cat ${collected_files.join(' ')} > combined_samples.fastq
    """
}

