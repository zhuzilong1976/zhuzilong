# Input data

This directory is not tracked by version control. Download the following files
before running the pipeline.

## Primary dataset

**GSE279180** — *Cell type mapping reveals tissue niches and interactions in
subcortical multiple sclerosis lesions [snRNA-seq]*, Homo sapiens.

```bash
# microglia subset (127.7 MB), raw integer counts
curl -L -o data/GSE279180_ctype_MG.h5ad \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279180/suppl/GSE279180_ctype_MG.h5ad
```

The file is an AnnData object with 23,627 genes × 9,239 microglia. `X` holds raw
counts (verified: all integers, range 1–1217). `obs` contains `condition`
(Control/MS), `lesion_type` (CA/CI/Ctrl), `patient_id`, `sample_id`, `subtype`,
`age`, `sex`, `batch_sn`.

Derived files created by `code/01_data/01_prepare_microglia_data.R`:

| File | Contents |
|---|---|
| `work/data/counts.mtx` | 23,627 genes × 9,239 cells, raw counts (Matrix Market) |
| `work/data/genes.txt` | gene symbols (rows) |
| `work/data/barcodes.txt` | cell barcodes (columns) |
| `work/data/metadata.csv` | barcode, celltype, condition, donor, patient, lesion_type, subtype, sample_id, age, sex, batch_sn |

## Replication networks

Ten network objects published by the developers of scTenifoldKnk, fetched
automatically by `code/04_replication/30_fetch_published_networks.mjs` into
`work/author_networks/` (108 MB total):

| File | Knockout | System | Network genes |
|---|---|---|---|
| `GSM3477499.RData` | Hnf4a + Hnf4g | mouse liver | 2,591 |
| `GSM3716703.RData` | Nkx2-1 | human lung | 8,647 |
| `GSM4116571.RData` | Dmd | human muscle | 9,783 |
| `MALAT1.RData` | Malat1 | human | 9,970 |
| `Preenterocytes.RData` | Ahr | mouse intestine | 8,478 |
| `SRS3059998.RData`, `SRS3059999.RData` | Mecp2 | mouse brain | 8,652 / 8,555 |
| `SRS3161261.RData`, `SRS4245406.RData` | Cftr | human lung AT2 | 6,585 / 7,107 |
| `SRS3161261_...Akap7.RData` | Akap7 (negative control) | human lung AT2 | 6,605 |

Source: `https://github.com/cailab-tamu/scTenifoldKnk/tree/master/inst/manuscript`.

## Gene-set libraries (for enrichment)

Fetched by `code/utils/fetch_gmt.mjs` (Enrichr REST API) into `work/gmt/`:
KEGG 2021 Human, Reactome 2022, GO Biological Process 2021, MSigDB Hallmark 2020.

## Base GRN (CellOracle)

`hg19_gimmemotifsv5_fpr2` (37,003 peaks × 1,094 TFs), 5.26 MB, downloaded
automatically by `celloracle.data.load_human_promoter_base_GRN()`:

```
https://raw.githubusercontent.com/morris-lab/CellOracle/master/celloracle/data/promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet
```

## Note on identifiers

Gene symbols are used as provided by each source (human symbols are upper-case,
mouse symbols are title-case). No cross-species harmonisation was performed.
