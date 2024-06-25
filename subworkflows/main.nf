/*

Include modules

*/

include { AGAT_FILTER_BY_LENGTH } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_PROTEINS } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_MINIPROT_NLR } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_NLR } from '../modules/HRP/agat/main'
include { AGAT_COMPLEMENT } from '../modules/HRP/agat/main'
include { MEME } from '../modules/HRP/memesuite/main'
include { MAST } from '../modules/HRP/memesuite/main'
include { BEDTOOLS_GETFASTA } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_CLUSTER } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_NR_CLUSTERS } from '../modules/HRP/bedtools/main'
include { INTERPROSCAN as INTERPROSCAN_FULL } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_EXTENDED } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_PFAM } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_SUPERFAMILY } from '../modules/HRP/interproscan/main'
include { FILTER_R_GENES as HRP_FILTER_R_GENES } from '../modules/HRP/local/main'
include { FILTER_R_GENES_SINGLE_INPUT as FIND_R_GENES } from '../modules/HRP/local/main'
include { SEQTK_SUBSET_RPS } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_FL } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_CANDIDATES } from '../modules/HRP/seqtk/main'
include { MINIPROT_HRP } from '../modules/HRP/miniprot/main'
include { SEQKIT_GET_LENGTH } from '../modules/HRP/seqkit/main'
include { GET_R_GENE_GFF } from '../modules/HRP/local/main'

log.info """\
  Parameters:
     samplesheet     : ${params.samplesheet}
     exlude_pattern  : ${params.exclude_pattern}
     outdir          : ${params.out}
"""
    .stripIndent(false)

workflow HRP {
    take: 
      hrp_in // tuple val(meta), path(fasta), path(gff)

    main:
      // Step 1 Extract proteins
      hrp_in
        .map { row -> [row.sample, row.fasta] }
        .set { genome }

      hrp_in
        .map { row -> [row.sample, row.gff] }
        .set { ref_gff }

      AGAT_EXTRACT_PROTEINS(hrp_in, params.exclude_pattern)

      AGAT_EXTRACT_PROTEINS
        .out
        .set { proteins }
       INTERPROSCAN_PFAM(proteins)
      //INTERPROSCAN_EXTENDED(proteins)
      // Step 3.1 Bedfile
      proteins
        .join(INTERPROSCAN_PFAM.out.nb_bed)
        .set { bedtools_gf_in }
      // Step 3.2 Extract
      BEDTOOLS_GETFASTA(bedtools_gf_in)
      // Step 3.3 MEME
      MEME(BEDTOOLS_GETFASTA.out)
      // Step 4 MAST
      MAST(proteins
            .join(MEME.out),
            params.mast_gene_pattern)
      // Step 5

      proteins
        .join(INTERPROSCAN_PFAM
                .out
                .protein_tsv
                .join(
                  MAST
                  .out
                  .mast_geneids
                  )
          )
        .set { to_subset }
                    
      SEQTK_SUBSET_RPS(to_subset)

      INTERPROSCAN_SUPERFAMILY(SEQTK_SUBSET_RPS.out)
      // Step 6
      HRP_FILTER_R_GENES(
        INTERPROSCAN_PFAM
        .out
        .protein_tsv
        .join(INTERPROSCAN_SUPERFAMILY.out)
        )
      // Step 7
      SEQTK_SUBSET_FL(proteins
                      .join(
                        HRP_FILTER_R_GENES
                        .out
                        .full_length_tsv)
                      )

      // Replacing genblast with miniprot
      genome
        .join(SEQTK_SUBSET_FL.out)
        .set { miniprot_in }

      miniprot_in.view()

      MINIPROT_HRP(miniprot_in) // emits gff

      AGAT_FILTER_BY_LENGTH(MINIPROT_HRP
                            .out
                            .miniprot_nlrs)

      genome
        .join(MINIPROT_HRP.out.miniprot_nlrs)
        .set { miniprot_nlr_to_extract }
      // Create proteins from gff
      AGAT_EXTRACT_MINIPROT_NLR(miniprot_nlr_to_extract)

      // Get lengths
      SEQKIT_GET_LENGTH(AGAT_EXTRACT_MINIPROT_NLR.out.extracted_nlrs)

      // Step 8.2
      BEDTOOLS_CLUSTER(AGAT_FILTER_BY_LENGTH
                        .out
                        .filtered_bed)

      // Step 8.3
      // Step 8.4
      BEDTOOLS_NR_CLUSTERS(BEDTOOLS_CLUSTER
                            .out
                            .join(SEQKIT_GET_LENGTH.out))

      // Step 9
      //   Extract annotations of non-redundant genes
      GET_R_GENE_GFF(AGAT_FILTER_BY_LENGTH
                      .out
                      .filtered_gff
                      .join(BEDTOOLS_NR_CLUSTERS.out))

      //   Extract protein sequences
      AGAT_EXTRACT_NLR(genome
                        .join(GET_R_GENE_GFF
                                .out
                                .r_genes_merged_gff))

      //   Merge R-Gene gff and input gff
      AGAT_COMPLEMENT(ref_gff
                        .join(GET_R_GENE_GFF
                                .out
                                .r_genes_merged_gff))

      AGAT_COMPLEMENT
        .out
        .merged_gff
        .set { merged_gff }

      AGAT_COMPLEMENT
        .out
        .merged_gtf
        .set { merged_gtf }

    emit:
        merged_gff
        merged_gtf

}