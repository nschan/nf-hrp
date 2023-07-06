include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AGAT_FILTER_BY_LENGTH {
  tag "$meta"
  label 'process_low'

  conda "bioconda::agat=1.1.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/agat:1.1.0--pl5321hdfd78af_1' :
          'quay.io/biocontainers/agat:1.1.0--pl5321hdfd78af_1' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(gff_file)
  
  output:
      tuple val(meta), path("*.gff"), emit: filtered_gff
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_sp_filter_gene_by_length.pl \\
  --gff ${gff_file} \\
  --size 20000 --test "<" \\
  -o ${meta}_genblastG-output_filtered.gff
  """
}

process AGAT_EXTRACT_PROTEINS {
  tag "$meta"
  label 'process_low'

  conda "bioconda::agat=1.1.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/agat:1.1.0--pl5321hdfd78af_1' :
          'quay.io/biocontainers/agat:1.1.0--pl5321hdfd78af_1' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(genome_fasta) ,path(genome_gff)
  
  output:
      tuple val(meta), path("*.fasta"), emit: extracted_proteins
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_sp_extract_sequences.pl \\
       -g ${genome_gff} \\
       -f ${genome_fasta} \\
       -p \\
       -o ${prefix}_proteins.fasta
  """
}

