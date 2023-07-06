include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_SUBSET {
  tag "$meta"
  label 'process_low'
  
  conda "bioconda::seqtk=1.6.4"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
          'quay.io/biocontainers/seqtk:1.4--he4a0461_1' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(fasta), path(ids1), path(ids2)
  
  output:
      tuple val(meta), path("*.fasta"), emit: fasta_subset

  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
    cat ${ids1} ${ids2} > all_ids.txt
    seqtk subseq ${fasta} all_ids.txt > ${meta}_subset.fasta
  """
}