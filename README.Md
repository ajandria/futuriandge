[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
## Introduction

Welcome to **DE_RNA-Seq_nf** - a bioinformatics best-practice analysis pipeline for up-stream processing of RNA-Seq data with the main focus on differential expression analysis.

<!-- TODO: Base on nf-core template -->

## Run test data as:
```
nextflow run ajandria/futuriandge -r main -latest -profile cluster,conda,test --organism mouse --outDir .
```

## Standard pipeline run:
```
nextflow run ajandria/futuriandge -r main -latest -profile cluster,singularity --organism mouse/human --outDir . --input samplesheet.csv
```

**Currently the pipeline supports only the MUB infrastructure (hardcoded references).**