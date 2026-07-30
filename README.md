# pySCENIC-analysis-of-Memory-B-and-Plasma-Cells
A pySCENIC analysis of Memory B and Plasma Cell Transcriptional Regulatory Networks in Healthy Bone Marrow and its Implications on the Myeloma Propagating Cell Debate


The project characterises the gene regulatory networks of normal Memory B cells and Plasma cells across a cohort of eight healthy bone marrow donors using single-cell RNA sequencing data from the Human Cell Atlas. Key analyses include gene regulatory network inference (pySCENIC), pseudobulk differential expression (PyDESeq2), Lasso regression modelling, pathway enrichment analysis (g:Profiler), and single-cell UMAP visualisation. Findings are interpreted in the context of the unresolved myeloma propagating cell (MPC) debate.

**Data Availability**

Input data (raw FASTQ files) were obtained from the Human Cell Atlas data portal: https://www.humancellatlas.org/
The dataset comprises 8 healthy bone marrow donors (4 female, 4 male, ages 26–52) sequenced using 10x Genomics Chromium 3' GEX chemistry. Raw data access requires registration with the Human Cell Atlas portal.

Reference files used:

Genome reference: GRCh38-2024-A (10x Genomics)


TF list: hs_hgnc_curated_tfs.txt (pySCENIC resources)


cisTarget databases: hg38 feather ranking databases (pySCENIC resources, available at https://resources.aertslab.org/cistarget/)


Motif table: motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl


Scripts should be run in the following order. Each step depends on the outputs of the previous step.


Step 0:  00_run_cellranger.sh           - filtered_feature_bc_matrix.h5 + BAM per donor


Step 1:  01_run_qc.sh                   - QC-filtered Seurat objects per donor


Step 2a: 02a_run_azimuth.sh             - per-donor Azimuth label CSVs


Step 2b: 02b_run_merge_labels.sh        - combined_annotated.h5ad + per-donor/cell-type h5ads


Step 3:  03_run_pyscenic.sh             - adjacencies.tsv + regulons.csv + auc_mtx.csv (×16)


Step 4:  04_run_grn_analysis.sh         - most_common_pairs.csv + violin plots


Step 5:  05_run_pseudobulk_de.sh        - balanced_cells.h5ad + de_memoryB_vs_plasma.csv


Step 6:  06_run_lasso.sh BACH2          - lasso_out/BACH2/*
         06_run_lasso.sh ATF4           - lasso_out/ATF4/*
         06_run_lasso.sh RORA           - lasso_out/RORA/*

         
Step 7:  07_run_gprofiler.sh            - gprofiler_out/ (CSVs + dot plots)


Step 8:  08_run_network_focused.sh      - network_out/ (BACH2 + ATF4 network PNGs)


Step 9: 10_run_overlap.sh               - overlap_out/ (overlap CSVs + figures)


Step 10: Run visualisation scripts      - figure PNGs for dissertation
