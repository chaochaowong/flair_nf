/*
*  Convert 'sample-sheet.csv' to FLAIR mandated TSV format using panda: 
*    - no header and sep by '\t'
*    - the first three columns should not have '_': convert all '_'
*      to '-' 
*/
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