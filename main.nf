#!/usr/bin/env nextflow

log.info """\
  FLAIR_NF PIPELINE on mamba and slurm executor 
  ================================================
  sample sheet      : ${params.sample_sheet}
  project directory : ${projectDir}
  work directory    : ${workDir}
  genomic reference : ${params.genome_reference}
  gtf               : ${params.gtf}
  output directory  : ${params.outdir}
"""
.stripIndent()

/*
* flair align: align fastq files to feed to minmap2; in the future we should
* allow hifi reads and use PacBio's wrapper function of minmap2
*/
process FLAIR_ALIGN {
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

process FLAIR_CONCAT_FASTQ {
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


process FLAIR_CONCAT_CORRECT_BED {
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
      --threads 8
    """
}

process SAMPLE_MANIFEST_TSV {
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
        path(sample_sheet)

    output:
        path('sample-manifest.tsv')

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    csv_file = '${sample_sheet}'
    df = pd.read_csv(csv_file)
    df[['sample_id', 'condition', 'batch']] = df[['sample_id', 'condition', 'batch']].replace('_', '-', regex=True)
    df.to_csv('sample-manifest.tsv', sep='\\t', index=False, header=False)
    """        
}

process FLAIR_QUANTIFY {
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
        path(sample_manifest_tsv)
        path(collapse_files)

    output:
        path('flair.quantify*')

    script:
    """
    flair quantify -r ${sample_manifest_tsv} \
      -i combined_samples.flair.collapse.isoforms.fa \
      --generate_map --isoform_bed combined_samples.flair.collapse.isoforms.bed \
      --stringent --check_splice \
      --threads 8 \
      --output flair.quantify
    """        
}

workflow {
    // Create input channel from samplesheet in TSV format
    reads_ch = Channel.fromPath(params.sample_sheet)
                      .splitCsv(header: true)
                      .map { row -> [row.sample_id, row.condition, row.batch, file(row.fq)] }
    // call flair align                      
    FLAIR_ALIGN(reads_ch, params.genome_reference, params.genome_reference_index, params.gtf, params.min_mapq)
    // flair correct
    FLAIR_CORRECT(FLAIR_ALIGN.out, params.genome_reference, params.gtf)
    // concatnate all fastq and all_corrected.bed
    FLAIR_CONCAT_FASTQ(params.sample_sheet)
    FLAIR_CONCAT_CORRECT_BED(FLAIR_CORRECT.out.collect())
    // flair collapse 
    FLAIR_COLLAPSE(params.genome_reference, params.gtf, FLAIR_CONCAT_FASTQ.out, FLAIR_CONCAT_CORRECT_BED.out)
    // sample manifest file
    SAMPLE_MANIFEST_TSV(params.sample_sheet)
    // flair quant
    FLAIR_QUANTIFY(SAMPLE_MANIFEST_TSV.out, FLAIR_COLLAPSE.out )
    FLAIR_QUANTIFY.out.view { it }
}