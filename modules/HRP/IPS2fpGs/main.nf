include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IPS2FPG {
  tag "$meta"
  label 'process_low'
  
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/debian:bullseye' :
        'debian:bullseye' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(proteins_fasta)
  
  output:
      tuple val(meta), path("*out.tsv"), emit: out_tsv
      tuple val(meta), path("*full_length.tsv"), emit: full_length_tsv
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
      def conf1 = file("$projectDir/assets/conf1.tsv", checkIfExists: true)
      def conf2 = file("$projectDir/assets/conf2.tsv", checkIfExists: true)
  """
  IPS2fpGs.sh -f -o ${meta}_ips2fp_out.tsv -c ${conf1} -d ${conf2} ${proteins_fasta}
  cat ${meta}_ips2fp_out.tsv | grep full- > ${meta}_ips2fp_out_full_length.tsv
  """
}