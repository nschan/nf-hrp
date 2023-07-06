include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEME {
  tag "$meta"
  label 'process_high'
  
  conda "bioconda::memesuite=5.5.3"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          '' :
          'memesuite/memesuite:5.5.3' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.txt"), emit: meme_out
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  meme ${protein_fasta}  \\
     -protein  \\
     -o ${prefix} \\
     -mod zoops \\
     -nmotifs 19 \\
     -minw 4 \\
     -maxw 7 \\
     -objfun classic \\
     -markov_order 0
  """
}

process MAST {
  tag "$meta"
  label 'process_high'
  
  conda "bioconda::memesuite=5.5.3"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          '' :
          'memesuite/memesuite:5.5.3' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(protein_fasta), path(meme_out)
  
  output:
      tuple val(meta), path("*.txt"), emit: mast_out
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  mast -o ${meta}_mast ${meme_out} ${protein_fasta}
  """
}