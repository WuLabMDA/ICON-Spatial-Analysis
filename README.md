# ICON-Spatial-Analysis

**Computational workflows accompanying the study:**

> **Integrated spatial multimodal analysis identifies multicellular niches associated with postsurgical recurrence and ctDNA shedding in early-stage non-small cell lung cancer**

## Overview

**ICON-Spatial** contains the complete computational workflows used to analyze the spatial organization of the tumor microenvironment in early-stage non-small cell lung cancer (NSCLC). The repository integrates multimodal spatial and molecular profiling to identify multicellular niches associated with recurrence, adaptive immune responses, genomic alterations, and circulating tumor DNA (ctDNA) shedding.

The analytical framework combines multiple orthogonal data modalities, including:

- Imaging Mass Cytometry (IMC)
- GeoMx Spatial Transcriptomics
- Bulk RNA Sequencing
- Whole-Exome Sequencing (WES)
- T-cell Receptor (TCR) Sequencing
- Circulating Tumor DNA (ctDNA)

Together, these analyses define clinically relevant multicellular spatial niches and their relationships with immune function, tumor genomics, T-cell clonal dynamics, and systemic biomarkers.

---

## Repository Structure

```
ICON-Spatial/
│
├── CellAnnotation/            # Cell phenotyping and annotation
├── CellularNeighborhoods/     # Cellular neighborhood (CN) identification
├── SpatialInteractions/       # Cell-cell interaction analyses
├── TLS/                       # Tertiary lymphoid structure analyses
├── NeutrophilHubs/            # Graph-based neutrophil-centered hub analysis
├── GeoMx/                     # Spatial transcriptomic analyses
├── BulkRNA/                   # Bulk RNA integration
├── Genomics/                  # WES and genomic association analyses
├── TCR/                       # T-cell receptor repertoire analyses
├── ctDNA/                     # ctDNA association analyses
├── Figures/                   # Figure generation scripts
└── README.md
```

---

## Major Analyses

The repository includes workflows for:

- IMC preprocessing and cell type annotation
- Identification of multicellular cellular neighborhoods (CNs)
- Spatial interaction analysis
- CellChat analysis
- Graph-based neutrophil-centered multicellular hub identification
- GeoMx spatial transcriptomic integration
- Bulk RNA-seq integration
- Whole-exome sequencing association analyses
- T-cell receptor repertoire analyses
- ctDNA association analyses
- Statistical analyses
- Figure generation

---

## Software

The analyses were primarily implemented in **R**.

Major packages include:

- SpatialExperiment
- cytomapper
- igraph
- CellChat
- Seurat
- edgeR
- limma
- data.table
- ggplot2
- ComplexHeatmap

---

## Citation

If you use this repository, please cite:

> **Integrated spatial multimodal analysis identifies multicellular niches associated with postsurgical recurrence and ctDNA shedding in early-stage non-small cell lung cancer.**

(*Nature Cancer*, in review)

---

## Contact

**Muhammad Aminu, Ph.D.**

Department of Imaging Physics

The University of Texas MD Anderson Cancer Center

Email: muhammadaminu47@gmail.com
