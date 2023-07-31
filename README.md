# nf-hrp

This pipeline implements the '*H*omology based *R* gene *P*rediction' workflow (https://github.com/AndolfoG/HRP; https://onlinelibrary.wiley.com/doi/10.1111/tpj.15756) in nextflow.

# Usage

The pipeline takes a gff3 file with annotations and a fasta file of the genome as input in a samplesheet and works through the steps of HRP.
The pipeline requires nextflow >= 23.04.2 to work with spack modules, tested on biohpc_gen.

Running this on other infrastructure probably needs some changes to modules.

Clone this repository into ~ (or somewhere else and adapt below) and do

```
nextflow run ~/nf-hrp/ --samplesheet sheet.csv --out './results' -profile biohpc_gen
```

samplesheet layout, header is important:

```
sample,fasta,gff
sample1,genome1.fasta,genome1.gff
```

# Output

The relevant gff files are `results/agat/sampleName_filtered_trascripts.[gff|bed]`, the output of the protein annotations from interproscan are in `results/interproscan/sampleName_NBLRR_gene_candidates.[gff3|tsv]`

# Notes

I ran into some issues during implementation, I am documenting these here.

## Intitial translation

I used this pipeline on automatic annotations, naturally some are broken. 

There is an option to remove a pattern from the initial gff file, e.g. I am removing "ATMG" to get rid of mitochondrial genes.

In addition, internal stop codons are stripped from protein translations, otherwise interproscan will fail. This means that there are probably some bad translations in the protein file.

I do not consider this a problem for the pipeline here, but I would not recommend using the protein.fasta generated here for other workflows.

## IPS2fpGs.sh

This file is bundled in bin/, the respective conf files are in assets/, original here https://github.com/AndolfoG/HRP

IPS2fpgGs.sh has been replaced with an R script that does the same, but is slower.
It is easier for me to maintain R than bash, the overall approach is similar, keeping proteins based on the presence of NB and LRR domains, and adding in extra info based on Coil, RPW8 or TIR domains. The output from this script is a bit more informative, it retains which domains were found in the protein in `_NLR_table.tsv`, `_NLR_genes.tsv` contains only gene-id and classification.

## Changes compared to the original pipeline

Subsetting of fasta files based on lists of gene IDs (input for step 4 and step 9) is done using seqtk.

Step 8.4 (sequence lengths) is done using seqkit.

I have added SSF52540 as an nucleotide binding domain into conf1.tsv (numeric code 7 for IPS2fpG). It appears that SUPERFAMILY does assign many of the NB domains into this category and not into PF00931.

## Interproscan

I was unable to run `interproscan` version 5.59_91.0 (the latest version on biocontainers) with `-dp` successfully: some of the interproscan pfam jobs with `-dp` errored out, and this produced a very short pfam table. 

I resorted to using spack, but I assume a container with the most recent version would suffice, currently that is 5.63-95.0

After I wrote this, a new container for interproscan:5.63-95.0 was released by interpro, I have added it to the charliecloud configuration.

## genblastG

genblastG is not maintained anymore, the official website does not seem to exist anymore.
I obtained geblastG from https://github.com/thecgs/genblastG_extension but I was unable to format a protein DB using formatdb (legacy blast db formatter), it removed valid letters from the input sequences (e.g. "E"?).

I could not figure out how to make formatdb work, instead I am using makeblastdb (recent blast+).
genblastG also tries to run ./formatdb, meaning that adding genblastG to the PATH is not sufficient, formatdb needs to be in the current folder.
This is "solved" by linking the whole genblastG folder into the work directory in the container.

The Dockerfile to build the container used is included in modules/HRP/genblastG.

I have not declared a conda package for this step, it appears the bioconda genblastG suffers similar issues.
