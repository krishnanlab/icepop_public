# ICePop: Metacell-informative Cell Population
This repository contains source code for ICePop ([DOI](TBD)). 

The data used in this study are available on Zenodo: https://doi.org/10.5281/zenodo.19238928

The code used to reproduce the analyses in the paper is available at: https://github.com/krishnanlab/icepop_analysis

## Dependencies
`python>=3.11,<3.12`

## Installation
ICePop can be installed easily via pip from PyPI: `pip install icepop`

## Run ICePop
Before running the analysis, we recommend downloading the processed data from [Zenodo](https://github.com/krishnanlab/icepop_analysis).

Expand and place the downloaded files under `../data`, then run the following commands.

A more detailed tutorial is available at [`notebook/ICePop_tutorial.ipynb`](https://github.com/krishnanlab/icepop_public/blob/main/notebook/ICePop_tutorial.ipynb)

### Step 1: Extract metacells
```
icepop metacell \
    --h5ad ../data/TM_FACS/TM_FACS_cnt.h5ad \
    --outdir ../results/TM_FACS \
    --save_name TM_FACS
```

#### Input options

1. `--h5ad` (str) Path to input AnnData (.h5ad) file containing single-cell expression ***count*** data
2. `--outdir` (str) Output directory where MetaQ results will be written
3. `--save_name` (str; default='metaq_res') prefix of metaq output under `./save/*`, do not write a path
4. `--ncell_per_mc` (int; default=75) Target number of cells per metacell. The total number of metacells is \n determined as approximately `n_cells / ncell_per_mc`
5. `--ct_key` (str; default='cell_type') Column name in `adata.obs` specifying cell-type annotations. Used to evaluate metacell purity
6. `--device` (str; default='cuda') Compute device to use. Options include 'cuda' or 'cpu'

this step need gpu for faster speed

#### Outputs
1. metacell assignment: `outdir/mc_assign.csv`
2. metacell statistics: `outdir/mc_stats.csv`

### Step 2: Get association, mixture and influence diagnoistics
```
icepop association \
    --h5ad ../data/TM_FACS/TM_FACS_cnt.h5ad \
    --mc_assign ../results/TM_FACS_mc/mc_assign.csv \
    --magmaz ../data/magmaz/asd.genes.out \
    --sp mmusculus \
    --outdir ../results/TM_FACS
```

#### Input options

1. `--h5ad` (str) Input AnnData file containing single-cell expression ***count*** data
2. `--mc_assign` (str) CSV file mapping cells to metacell assignments (output from step 1: `outdir/mc_assign.csv`)
3. `--magmaz` (str) [magmaz](https://doi.org/10.1371/journal.pcbi.1004219) MAGMA gene-level association file (*.genes.out) of a trait of interest
4. `--spec_score` (str; default=None) Precomputed specificity scores; will be calculated if not provided
5. `--outdir` (str) Output directory for association results
6. `--n_jobs` (int; default=20) Number of parallel workers
7. `--sp` (str; default='mmusculus') Species identifier for gene ID conversion
8. `--ct_key` (str; default='cell_type') Column in `adata.obs` defining cell types
9. `--trait_name` (str; optional) Trait name used for output file naming
10. `--n_perm` (int; default=1000) Number of permutations for null distribution estimation
11. `--q_thres` (float; default=0.1) FDR threshold for significance
12. `--output_dfbs` (boolean; default=True) If output influential testing results

#### Outputs

1. `outdir/celltype__trait-*.csv`: Disease-cell type association table
2. `outdir/dfbs__trait-*.npz`: Gene-level influence scores (DFBETAS) for each disease–cell type association
3. `outdir/metacell__trait-*.csv`: Disease-metacell type association table
4. `outdir/mc_spec_score.npz`: (if not provided as part of run arguments)
5. `outdir/mcfdr__trait-*.csv`: Cell type × metacell matrix indicating significant disease-associated metacells within each cell type

where `*` is trait name we assume magmaz file name is `*.genes.out`

### Step3: Enrichment Analysis and Interactive output
```
icepop interactive \
  --outdir ../results/TM_FACS \
  --geneset_collections KEGG \
  --adata_path ../data/TM_FACS/TM_FACS_cnt.h5ad

or 

icepop interactive \
  --outdir ../results/TM_FACS \
  --geneset_collections none \
  --geneset_path custom.gmt \
  --adata_path ../data/TM_FACS/TM_FACS_cnt.h5ad
```

#### Input options
1. `--outdir` (str) Output directory for association results and metacell results
2. `--geneset_collections` (str) All, 'BIOCARTA', 'KEGG', 'REACTOME', 'WIKIPATHWAYS', 'MIR', 'TF', 'GOBP', 'GOCC', 'GOMF', 'HP'
3. `--geneset_path` (str) path to custom gmt file
4. `--adata_path` (str) path to AnnData file containing single-cell expression ***count*** data

#### Outputs
1. `outdir/icepop-report.ipynb`: Interactive Jupyter notebook containing all results
2. `outdir/icepop-report.html`: Rendered HTML version of the notebook for easy viewing
3. `outdir/enrichment`: Directory containing gene set enrichment analysis results

> **Note:** We recommend using the same `--outdir` for Step 1 and Step 2, as the interactive step expects results from both steps to be located in the same output directory.