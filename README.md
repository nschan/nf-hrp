# nf-hrp

This pipeline implements the homology based R gene prediction workflow (https://github.com/AndolfoG/HRP) in nextflow.

## Usage

This pipeline takes a gff3 file with annotations and a fasta file as input in a samplesheet and works through the steps of HRP.
This pipeline requires nextflow >= 23.04.2 to work with spack modules, the profile is quite specific for biohpc_gen.

```
nextflow run ~/nf-hrp/ --samplesheet sheet.csv --out './results' -profile biohpc_gen
```

samplesheet layout, header is important:
```
sample,fasta,gff
sample1,genome1.fasta,genome1,gff
```