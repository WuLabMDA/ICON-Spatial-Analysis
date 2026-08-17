# ICON-Spatial-Analysis

**Codes accompanying the study:**

> **Integrated spatial multimodal analysis identifies multicellular niches associated with postsurgical recurrence and ctDNA shedding in early-stage non-small cell lung cancer**

## Overview

**ICON-Spatial-Analysis** contains the computational workflows used to analyze the spatial organization of the tumor microenvironment in early-stage non-small cell lung cancer (NSCLC). The repository integrates multimodal spatial and molecular profiling to identify multicellular niches associated with recurrence, adaptive immune responses, genomic alterations, and circulating tumor DNA (ctDNA) shedding.

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
├── Codes/
│   ├── Landscape/                      # Cell-cell interaction analyses
│   ├── Transcriptomic validation/      # IMC, GeoMx, Bulk RNA-seq integration
│   ├── Multicellular neighborhoods/    # Cellular neighborhood (CN) identification, TLS analysis
│   ├── TCR analysis/                   # T-cell receptor repertoire analyses
│   ├── ctDNA analysis/                 # Circulating tumor DNA analyses, Neutrophil-centered hub analysis
└── README.md
```

---

## Major Analyses

The repository includes workflows for:

- Spatial interaction analysis
- Identification of multicellular cellular neighborhoods (CNs)
- Tertiary lymphoid structure analyses
- GeoMx spatial transcriptomic integration
- Bulk RNA-seq integration
- Whole-exome sequencing association analyses
- T-cell receptor repertoire analyses
- ctDNA analyses
- Neutrophil-centered multicellular hub identification

---

## Software

The analyses were primarily implemented in **R**.

Major packages include:

- SpatialExperiment
- cytomapper
- igraph
- Seurat
- edgeR
- imcRtools

---

## Citation

If you use this repository, please cite:

> **Integrated spatial multimodal analysis identifies multicellular niches associated with postsurgical recurrence and ctDNA shedding in early-stage non-small cell lung cancer.**

(*TBD*, in review)


Email: muhammadaminu47@gmail.com
