/*
Here I attempt to replicate the HRP workflow
https://github.com/AndolfoG/HRP/
All modules are in modules/HRP
*/
nextflow.enable.dsl = 2 
params.out = './results'
params.publish_dir_mode = 'copy'

include { AGAT_FILTER_BY_LENGTH } from './modules/HRP/agat/main'
include { AGAT_EXTRACT_PROTEINS } from './modules/HRP/agat/main'
include { MEME } from './modules/HRP/memesuite/main'
include { MAST } from './modules/HRP/memesuite/main'
include { BEDTOOLS_GETFASTA } from './modules/HRP/bedtools/main'
include { BEDTOOLS_CLUSTER } from './modules/HRP/bedtools/main'
include { BEDTOOLS_NR_CLUSTERS } from './modules/HRP/bedtools/main'
include { INTERPROSCAN } from './modules/HRP/interproscan/main'
include { INTERPROSCAN_PFAM } from './modules/HRP/interproscan/main'
include { INTERPROSCAN_SUPERFAMILY } from './modules/HRP/interproscan/main'
include { IPS2FPG } from './modules/HRP/IPS2fpGs/main'
include { SEQTK_SUBSET} from './modules/HRP/seqtk/main'

workflow HRP {
    take: 
      hrp_in // tuple val(meta), path(genome_fasta), path(genome_gff)

    main:
      // Step 1 Extract proteins
      AGAT_EXTRACT_PROTEINS(hrp_in)
      // Step 2 Interproscan
      // This step works with spack module interproscan/5.63-95.0-gcc11-csy
      // It does not seem to give proper results with biocontainers/interproscan:5.59_91.0--hec16e2b_1
      // This seems to be version related, there is no container for 5.63-95.0 and I cannot build it, attempt to use spack for this
      INTERPROSCAN_PFAM(AGAT_EXTRACT_PROTEINS.out)
      // Step 3.1 Bedfile
      bedtools_gf_in = AGAT_EXTRACT_PROTEINS.out.join(INTERPROSCAN_PFAM.out.nb_bed)
      // Step 3.2 Extract
      BEDTOOLS_GETFASTA(bedtools_gf_in)
      // Step 3.3 MEME
      MEME(BEDTOOLS_GETFASTA.out)
      // Step 4 MAST
      MAST(AGAT_EXTRACT_PROTEINS.out.join(MEME.out))
      // Step 5
      // UNCLEAR WHERE SUBSET COMES FROM
      // Subset should be created from
      //to_subset = INTERPROSCAN_PFAM.out.protein_tsv.join(MAST.OUT)
      //SEQTK_SUBSET()
      //INTERPROSCAN_SUPERFAMILY()
      // Step 6
      //IPS2FPG(INTERPROSCAN_SUPERFAMILY.out)
      // Step 7

      // Step 8.1
      // Step 8.2
      // Step 8.3
      // Step 8.4
      // Step 9
}

workflow {
  ch_input = Channel.fromPath(params.samplesheet) 
                      .splitCsv(header:true) 
  HRP(ch_input)
}