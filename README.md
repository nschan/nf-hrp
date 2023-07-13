# nf-hrp

This pipeline implements the '*H*omology based *R* gene *P*rediction' workflow (https://github.com/AndolfoG/HRP; https://onlinelibrary.wiley.com/doi/10.1111/tpj.15756) in nextflow.

# Usage

The pipeline takes a gff3 file with annotations and a fasta file as input in a samplesheet and works through the steps of HRP.
The pipeline requires nextflow >= 23.04.2 to work with spack modules, tested on biohpc_gen.

Running this on other infrastructure probably needs some changes to modules.

```
nextflow run ~/nf-hrp/ --samplesheet sheet.csv --out './results' -profile biohpc_gen
```

samplesheet layout, header is important:

```
sample,fasta,gff
sample1,genome1.fasta,genome1.gff
```

# Output

The relevant gff files are `results/agat/*filtered_trascripts.[gff|bed]`, the output of the protein annotations from interproscan are in `results/interproscan/*NBLRR_gene_candidates.[gff3|tsv]`

# Notes

I ran into some issues during implementation, I am documenting these here.

## IPS2fpGs.sh

This file is bundled in bin/, the respective conf files are in assets/, original here https://github.com/AndolfoG/HRP

## Changes compared to the original pipeline

Subsetting of fasta files based on lists of gene IDs (input for step 4 and step 9) is done using seqtk.

Step 8.4 (sequence lengths) is done using seqkit.

## Interproscan

I was unable to run `interproscan` version 5.59_91.0 (the latest version on biocontainers) with `-dp` successfully: some of the interproscan pfam jobs with `-dp` errored out, and this produced a very short pfam table. 

I resorted to using spack, but I assume a container with the most recent version would suffice, currently that is 5.63-95.
A dockerfile to build such a container is included in modules/HRP/interproscan.

## genblastG

genblastG is not maintained anymore, the official website does not seem to exist anymore.
I obtained geblastG from https://github.com/thecgs/genblastG_extension but I was unable to format a protein DB using formatdb (legacy blast db formatter), it removed valid letters from the input sequences (e.g. "E"?).

I could not figure out how to make formatdb work, instead I am using makeblastdb (recent blast+).
genblastG also tries to run ./formatdb, meaning that adding genblastG to the PATH is not sufficient, formatdb needs to be in the current folder.
This is "solved" by linking the whole genblastG folder into the work directory in the container.

The Dockerfile to build the container used is included in modules/HRP/genblastG.
