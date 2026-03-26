# ICePop: Informative Cell Population
This repository contains source code for ICePop ([DOI](TBD))

## Dependencies
`python>=3.11,<3.12`
outdir/mc_spec_score.npz (if not provided as part of run arguments)

## Installation
ICePop can be installed easily via pip from PyPI:`pip install icepop`

## Run ICePop
### Step 1: Extract metacells
```
icepop metacell \
    --h5ad ../data/mouse_colon/mouse_colon_cnt.h5ad \
    --outdir ../results/mouse_colon_mc \
    --save_name mouse_colon
```

#### Input options

1. `--h5ad` (str) Path to input AnnData (.h5ad) file containing single-cell expression data
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
    --magmaz ../data/TM_FACS/magmaz/asd.genes.out \
    --spec_score ../results/TM_FACS_mc/mc_spec_score.npz \
    --sp mmusculus \
    --outdir ./test
```

#### Input options

1. `--h5ad` (str) Input AnnData file containing single-cell expression data
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
