/*
* flair collapse for isoform identification
*
* The important output files are
*  - prefix.isoforms.gtf - your custom transcriptome which you can align to if you want
*  - prefix.isoforms.bed - the easiest way to visualize your isoforms on the UCSC genome browser or IGV, can also be useful for FLAIR-quantify
*  - prefix.combined.isoform.read.map.txt - all detected isoforms associated with the reads that support them
*/

process FLAIR_COLLAPSE {
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input:
        path(ref_fasta)
        path(gtf)
        path(combined_fastq)
        path(combined_corrected_bed)

    output:
        path('combined_samples.flair.collapse*')        

    script: 
    """
    flair collapse -g ${ref_fasta} \
      --gtf ${gtf} \
      -q ${combined_corrected_bed} \
      -r ${combined_fastq} \
      --annotation_reliant generate --generate_map --check_splice --stringent \
      --output combined_samples.flair.collapse \
      --threads $task.cpus
    """
}
