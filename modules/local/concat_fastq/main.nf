/*
* Concatenate all fastq files in the 'fq' colum and deposit
# to ./collapse/combined_sample.fastq
*/
process CONCAT_FASTQ {
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input: 
        path(sample_sheet)

    output:    
        path('combined_samples.fastq')

    script:
    """
    awk -F, 'NR>1 {print \$4}' ${sample_sheet} | xargs cat > combined_samples.fastq
    """    
}
