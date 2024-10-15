/*
* Concatenate all corrected bed from the FLAIR_CORRECT channel
* and deposit to './collapse/combined_samples.all_corrected.bed'
*/
process CONCAT_CORRECTED_BED {
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input: 
        val(correct_files)

    output: 
        path('combined_samples.all_corrected.bed')

    script:
    """
    # Use Groovy's findAll to extract .bed files
    bed_files=\$(echo ${correct_files.join(' ')} | xargs -n 1 | grep 'flair_all_corrected.bed')

    # Concatenate the filtered bed files
    cat \$bed_files > combined_samples.all_corrected.bed
    """        

}