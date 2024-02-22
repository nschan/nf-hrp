/*
Original description:
https://github.com/AndolfoG/HRP/
All modules are in modules/HRP
*/

nextflow.enable.dsl = 2 
params.out = './results'
params.publish_dir_mode = 'copy'
params.exclude_pattern = "ATMG"

include { HRP } from './subworkflows/main'

workflow {
  Channel.fromPath(params.samplesheet) 
    .splitCsv(header:true)
    .set { ch_input }
  HRP(ch_input)
}