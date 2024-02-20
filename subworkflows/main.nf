/*

Include modules

*/

include { AGAT_FILTER_BY_LENGTH } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_PROTEINS } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_NLR } from '../modules/HRP/agat/main'
include { AGAT_COMPLEMENT } from '../modules/HRP/agat/main'
include { MEME } from '../modules/HRP/memesuite/main'
include { MAST } from '../modules/HRP/memesuite/main'
include { BEDTOOLS_GETFASTA } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_CLUSTER } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_NR_CLUSTERS } from '../modules/HRP/bedtools/main'
include { INTERPROSCAN } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_EXTENDED } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_PFAM } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_SUPERFAMILY } from '../modules/HRP/interproscan/main'
/*
  IPS2FPG no longer included and replaced by R script in FILTER_R_GENES
  include { IPS2FPG } from './modules/HRP/IPS2fpGs/main'
*/
include { FILTER_R_GENES } from '../modules/HRP/local/main'
include { SEQTK_SUBSET_RPS } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_FL } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_CANDIDATES } from '../modules/HRP/seqtk/main'
include { GENBLAST_G } from '../modules/HRP/genblastG/main'
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
      /* This procedure produces a few broken proteins on TAIR11.
         This is somewhat concerning, but.. well.
      */ 
      genome = hrp_in.map(row -> [row.sample, row.fasta])
      ref_gff = hrp_in.map(row -> [row.sample, row.gff])

      AGAT_EXTRACT_PROTEINS(hrp_in, params.exclude_pattern)
      AGAT_EXTRACT_PROTEINS
        .out
        .set { proteins }
      // Step 2 Interproscan
      // This step works with spack module interproscan/5.63-95.0
      // I could not locate a container with this version.
      // Internal container build pipeline failed, the container exceeds storage on our gitlab container storage.
      //
      // It does not seem to give proper results with biocontainers/interproscan:5.59_91.0--hec16e2b_1
      // This seems to be version related.
      // For me interproscan:5.59_91.0 with -dp did not run, subjobs failed and the result was incomplete.
      // To use interproscan without -dp the version needs to be the current release.
      // As of July 2023 this is 5.63-95-0
      // I guess there are work-arounds for this, it should work with an updated container.
      INTERPROSCAN_PFAM(proteins)
      // Step 3.1 Bedfile
      proteins
        .join(INTERPROSCAN_PFAM.out.nb_bed)
        .set { bedtools_gf_in }
      // Step 3.2 Extract
      BEDTOOLS_GETFASTA(bedtools_gf_in)
      // Step 3.3 MEME
      MEME(BEDTOOLS_GETFASTA.out)
      // Step 4 MAST
      MAST(proteins.join(MEME.out))
      // Step 5
      proteins
        .join(INTERPROSCAN_PFAM.out.protein_tsv
        .join(MAST.out.mast_geneids))
        .set { to_subset }

      SEQTK_SUBSET_RPS(to_subset)

      INTERPROSCAN_SUPERFAMILY(SEQTK_SUBSET_RPS.out)
      // Step 6
      //IPS2FPG(INTERPROSCAN_PFAM.out.protein_tsv.join(INTERPROSCAN_SUPERFAMILY.out))
      FILTER_R_GENES(INTERPROSCAN_PFAM.out.protein_tsv.join(INTERPROSCAN_SUPERFAMILY.out))
      // Step 7
      SEQTK_SUBSET_FL(proteins
                      .join(FILTER_R_GENES.out.full_length_tsv))
      // Genblast
      genome
        .join(SEQTK_SUBSET_FL.out)
        .set { genblast_in }
        
      GENBLAST_G(genblast_in)
      SEQKIT_GET_LENGTH(GENBLAST_G.out.genblast_pro)
      // Step 8.1
      AGAT_FILTER_BY_LENGTH(GENBLAST_G.out.genblast_gff)
      // Step 8.2
      BEDTOOLS_CLUSTER(AGAT_FILTER_BY_LENGTH.out.filtered_bed)
      // Step 8.3
      // Step 8.4
      BEDTOOLS_NR_CLUSTERS(BEDTOOLS_CLUSTER.out.join(SEQKIT_GET_LENGTH.out))
      // Step 9
      //   Extract annotations of non-redundant genes
      GET_R_GENE_GFF(AGAT_FILTER_BY_LENGTH.out.filtered_gff.join(BEDTOOLS_NR_CLUSTERS.out))
      //   Extract protein sequences
      AGAT_EXTRACT_NLR(genome.join(GET_R_GENE_GFF.out.r_genes_merged_gff))
      //   Merge R-Gene gff and input gff
      AGAT_COMPLEMENT(ref_gff.join(GET_R_GENE_GFF.out.r_genes_merged_gff))
      //   Interproscan of NLR-Candidates
      // INTERPROSCAN(AGAT_EXTRACT_NLR.out)
}