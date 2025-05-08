# A.calliptera-ATAC-analysis
This repository contains code, results, and documentation for an ATAC-seq analysis project exploring chromatin accessibility across three developmental stages (3 dpf, 7 dpf, and 12 dpf) in Astatotilapia calliptera.


Introduction
Astatotilapia calliptera, a generalist cichlid fish species, serves as a valuable model for studying epigenetic mechanisms underlying phenotypic plasticity and environmental adaptation. This project aims to explore how chromatin accessibility patterns change during early development (3, 7, and 12 days post-fertilization) using ATAC-seq. By identifying regulatory elements and transcription factors involved at each developmental stage, the study seeks to provide insights into the genetic and epigenetic basis of developmental regulation in this species.

Methods
Sample Collection & ATAC-seq Data Acquisition: ATAC-seq data from A. calliptera embryos at 3, 7, and 12 dpf were obtained.

Peak Calling & Differential Accessibility: Quality control, alignment to the reference genome, and peak calling were performed. Differential peak analysis across developmental stages was conducted using edgeR.

Functional Annotation: Peaks were annotated using GAT (Genomic Association Tester) to determine enrichment in genomic features and nearby genes.

GO Enrichment & Pathway Analysis: Genes associated with differentially accessible regions were analyzed using g:Profiler to identify overrepresented biological processes and pathways.

Transcription Factor Analysis: TOBIAS was used for TF footprinting and binding dynamics to identify stage-specific TF activity and regulatory networks.

Clustering: TF footprints were clustered to reveal patterns of regulatory coordination across development.

Expected Results
Identification of dynamic changes in chromatin accessibility that correspond with key developmental transitions.

Enrichment of open chromatin regions in promoter and enhancer elements of developmental and regulatory genes.

Stage-specific activation of transcription factor families involved in early patterning, organogenesis, and neural development.

Insight into the gene regulatory networks that support developmental plasticity in a generalist fish species.
