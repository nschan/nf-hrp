/*
Original description:
https://github.com/AndolfoG/HRP/
All modules are in modules/HRP
*/

nextflow.enable.dsl = 2 
params.out = './results'
params.publish_dir_mode = 'copy'
params.cds_feature = "CDS"
params.exclude_pattern = "ATMG"
/* 
Define pattern for MAST, which is used to extract gene sequences
The default extracts TAIR names or evm names.
*/
params.mast_gene_pattern = "AT[1-5C]G[0-9]+.[0-9]+|evm[0-9a-z\\.]*"

include { HRP } from './subworkflows/main'

workflow {
  Channel.fromPath(params.samplesheet) 
    .splitCsv(header:true)
    .set { ch_input }
  HRP(ch_input)
}