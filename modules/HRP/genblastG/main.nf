include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENBLAST_G {
  tag "$meta"
  label 'process_high'
  container "gitlab.lrz.de:5005/beckerlab/container-playground/genblastg:ff448d5c"
  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(nb_lrr_fasta), path(genome_fasta)
  
  output:
      tuple val(meta), path("*.pro"), emit: genblast_pro
      tuple val(meta), path("*.txt"), emit: length_estimates
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
      
  """
  ln -s /opt/genblastg_extension/* .

  genblastG -q ${nb_lrr_fasta} -t ${genome_fasta} -gff -cdna -pro -o ${prefix}_genblastG-output

  awk 'BEGIN{FS="[> ]"} /^>/{val=\$2;next}  {print val,length(\$0);val=""} END{if(val!=""){print val}}' ${prefix}_genblastG-output.pro | tr ' ' \\t > ${prefix}_genblastG-output_FbL_length.txt
  """
}