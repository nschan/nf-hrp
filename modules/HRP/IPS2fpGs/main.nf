include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IPS2FPG {
  tag "$meta"
  label 'process_low'
  
  conda "conda-forge::python=3.9.5"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(proteins_fasta)
  
  output:
      tuple val(meta), path("*.tsv"), emit: out_tsv
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  IPS2FPGs ${proteins_fasta} -o ${meta}_ips2fp_out.tsv
  """
}