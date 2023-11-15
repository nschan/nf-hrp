include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process INTERPROSCAN_PFAM {
  tag "$meta"
  label 'process_high'

  spack 'interproscan@5.63-95.0'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/'), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
      tuple val(meta), path("*.bed"), emit: nb_bed
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  interproscan.sh \\
     -f TSV,GFF3 \\
     -appl Pfam \\
     -cpu $task.cpus \\
     -i ${protein_fasta} \\
     -b ${prefix}_proteins
  grep NB-ARC ${prefix}_proteins.tsv | cut -f1,7,8 > ${prefix}_NB.bed
  """
}

process INTERPROSCAN {
  tag "$meta"
  label 'process_high'

  spack 'interproscan@5.63-95.0'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/'), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(candidate_nb_lrrs)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
      tuple val(meta), path("*.gff3"), emit: protein_gff
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  interproscan.sh \\
     -f TSV,GFF3 \\
     -cpu $task.cpus \\
     -i ${candidate_nb_lrrs} \\
     -b ${prefix}_NBLRR_gene_candidates     
  """
}

process INTERPROSCAN_SUPERFAMILY {
  tag "$meta"
  label 'process_high'

  spack 'interproscan@5.63-95.0'

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  interproscan.sh \\
    -f TSV \\
    -cpu $task.cpus \\
    -app SUPERFAMILY \\
    -i ${protein_fasta} \\
    -b ${prefix}_superfam
  """
}