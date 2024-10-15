/*
* flair quantify process that creates flair.quantify.*.isoform.read.map.txt and
* flair.quanitfy.counts.tsv files in the ./quant foler
* NOTE: the parameters suggested here are from the tutorial and for hg38
*/
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
      --threads $task.cpus \
      --output flair.quantify

    echo "TMPDIR is set to \$TMPDIR"      
    """        
}
