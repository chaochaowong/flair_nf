process FLAIR_ALIGN {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input: tuple val(sample_id), val(sample), path(fastq)
        path ref_fasta
        path ref_index
        path gtf
    
    output:
        tuple val(sample_id), path("${sample_id}.flair.aligned")

    """
    flair align -g ${ref_fasta} \
        -r ${input_fastq} \
        -o ${sample_id}.flair.aligned
    """        
}

workflow {
    // Create input channel from samplesheet in CSV format
    reads_ch = Channel.fromPath(params.sample_sheet)
                      .splitCsv(header: true)
                      .map { row -> [row.sample_id, row.sample, file(row.fastq)] }
}