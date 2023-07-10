include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQKIT_GET_LENGTH {
  tag "$meta"
  label 'process_low'

  conda "bioconda::seqkit=2.4.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0' :
          'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(fasta_file)
  
  output:
      tuple val(meta), path("*_length.txt"), emit: length_estimates
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  seqkit fx2tab --length --name ${fasta_file} > ${meta}_length.txt
  """
}

process AGAT_EXTRACT_PROTEINS {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(genome_fasta), path(genome_gff)
  
  output:
      tuple val(meta), path("*.fasta"), emit: extracted_proteins
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${genome_gff} | grep CDS | grep -v ChrM > ${genome_gff}_subset.gff3
  agat_sp_extract_sequences.pl \\
       -g ${genome_gff}_subset.gff3 \\
       -f ${genome_fasta} \\
       -p \\
       -o proteins.fa
  
  perl -pe 'if (/^>/) { \$. > 1 and print "\n" } else { chomp }' < proteins.fa  | sed 's/.\$//' > proteins_nostop.fa 
  grep -B1 "*" proteins_nostop.fa | grep -vFf - proteins_nostop.fa | fold > ${prefix}_proteins.fasta  
  """
}

