# osteosarcoma-Mn-programs

Analysis code for:

**Cell-type-specific distribution of manganese-related gene programs in the osteosarcoma microenvironment: an integrative single-cell, spatial and Mendelian randomization study**

Yizhe He, Meiyun Zhao, Zhi Zeng, Xiaoxu Wang. The Second Affiliated Hospital of the University of South China, Hengyang, China.

The study is purely computational and uses only publicly available, de-identified data. No data files are stored in this repository; see *Data* below.

## Scripts

Run from the project root, in this order.

| Script | Content |
|---|---|
| `scripts/00_main_all_lessons.R` | Main analysis: data download and processing, curated Mn-related tiers (T1–T4), two-sample and multivariable Mendelian randomization, single-cell QC, integration, annotation and tier scoring (discovery and validation cohorts), patient-level pseudobulk tests, CellChat, bulk cohorts (TARGET-OS, GSE21257), Visium spatial analysis, drug–gene interactions and clinical-trial search |
| `scripts/19_revision_checks.R` | Revision analyses and final main figures (Figs 1–8): figure theme, patient-level re-analysis of the spatial data (19 patients), MR with FinnGen R13, MVMR diagnostics, survival analyses and random gene-set controls |
| `scripts/20_supplementary_figures_and_spatial_revision.R` | Supplementary Figures S1–S4, Additional file 1, the lesion-level spatial analyses with 25 lesions (MR10A and MR10B treated as separate lesions; MR11A and MR11B combined) and the final Figure 5 |

The scripts were written step by step during the analysis and contain intermediate checks. Lines that only acted on the authors' local file system (renaming, opening files, setting the working directory) are commented out and marked `[local-only, disabled for public release]`.

## Data

| Data | Source |
|---|---|
| Discovery single-cell RNA-seq (11 patients) | GEO GSE152048 |
| Validation single-cell RNA-seq (6 treatment-naive patients) | GEO GSE162454 |
| Visium spatial transcriptomics (26 sections, 19 patients) | GEO GSE293065 |
| Bulk expression and clinical data | GEO GSE21257; TARGET-OS via the UCSC Xena GDC hub |
| Whole-blood Mn GWAS (HUNT, MoBa, PIVUS) | DataverseNO, doi:10.18710/UYPCL0 |
| Osteosarcoma and bone cancer GWAS | FinnGen release R13 (C3_OSTEOSARCOMA_WIDE, C3_BONE_CARTILAGE) |
| Serum iron GWAS | IEU OpenGWAS, ieu-a-1049 |

Raw data are expected under `data_raw/`, intermediate objects are written to `data_proc/`, results to `results/` and final figures to `results/figures/final/` and `final/`.

Access to IEU OpenGWAS requires a personal token. Store it in your `.Renviron` file as `OPENGWAS_JWT=...`; it is read with `Sys.getenv("OPENGWAS_JWT")` and is never written in the scripts.

## Main software

R with Seurat, harmony, scDblFinder, CellChat, TwoSampleMR (v0.7.14), MRPRESSO, MVMR, ieugwasr, nlme, survival, data.table, ggplot2, patchwork, ggrastr, magick and qpdf. Package versions are printed in `scripts/19_revision_checks.R` (step 494).

## License

MIT (see `LICENSE`).
