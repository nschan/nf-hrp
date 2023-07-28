include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_R_GENE_GFF {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"gff_subset", publish_id:meta) }

  input:
      tuple val(meta), path(gff), path(r_gene_list)
  
  output:
      tuple val(meta), path("*R_gene_annotations.gff"), emit: r_gene_gff
      tuple val(meta), path("*_R_genes_merged.gff"),    emit: r_genes_merged_gff
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  grep -Fw -f ${r_gene_list} ${gff} > ${meta}_R_gene_annotations.gff
  cat ${meta}_R_gene_annotations.gff \\
  | sort -k1,1 -k4,4n \\
  | bedtools merge -s -c 2,3,7,8,9 -o distinct,distinct,distinct,distinct,collapse > ${meta}_R_genes_merged.gff
  """
}