# ICePop: Informative Cell Population
This repository contains source code for ICePop ([DOI](TBD))

## Installation
### Dependencies
ICePop dependencies are handled by [poetry](https://python-poetry.org/) with `python>=3.11,<3.12`
To install poetry, please follow the [instructions on poetry's home page](https://python-poetry.org/docs/#installation).
Then run the following command to 
- create a virtual environment
- install all dependencies from `poetry.lock`
- install icepop in editable mode:
```
poetry install
```

### [Optional] Using icepop in another local project 
If you want to use this repository as a local editable dependency in a different Poetry project:
```
poetry add --editable /path/to/icepop
```
where the path must point to the folder containing `pyproject.toml`

## Run ICePop
### Step 1: Extract metacells
```
poetry run icepop metacell \
    --h5ad ../data/mouse_colon/mouse_colon_cnt.h5ad \
    --outdir ../results/mouse_colon_mc \
    --save_name mouse_colon
```

#### input options

1. `--h5ad` (str) Path to input AnnData (.h5ad) file containing single-cell expression data
2. `--outdir` (str) Output directory where MetaQ results will be written
3. `--save_name` (str; default='metaq_res') prefix of metaq output under `./save/*`, do not write a path
4. `--ncell_per_mc` (int; default=75) Target number of cells per metacell. The total number of metacells is \n determined as approximately `n_cells / ncell_per_mc`
5. `--ct_key` (str; default='cell_type') Column name in `adata.obs` specifying cell-type annotations. Used to evaluate metacell purity
6. `--device` (str; default='cuda') Compute device to use. Options include 'cuda' or 'cpu'

this step need gpu for faster speed

#### output
metacell assignment

### Step 2: Get association, mixture and influence diagnoistics
```
poetry run icepop association \
    --h5ad ../data/TM_FACS/TM_FACS_cnt.h5ad \
    --mc_assign ../results/TM_FACS_mc/mc_assign.csv \
    --magmaz ../data/TM_FACS/magmaz/asd.genes.out \
    --spec_score ../results/TM_FACS_mc/mc_spec_score.npz \
    --outdir ./test
```

#### input options

1. `--h5ad` (str) Input AnnData file containing single-cell expression data
2. `--mc_assign` (str) CSV file mapping cells to metacell assignments
3. `--magmaz` (str) [magmaz](https://doi.org/10.1371/journal.pcbi.1004219) summary statistics
4. `--spec_score` (str) MAGMA gene-level association file (*.genes.out)
5. `--outdir` (str) Output directory for association results
6. `--n_jobs` (int; default=20) Number of parallel workers
7. `--sp` (str; default='mmusculus') Species identifier for gene ID conversion
8. `--ct_key` (str; default='cell_type') Column in `adata.obs` defining cell types
9. `--trait_name` (str; optional) Trait name used for output file naming
10. `--n_perm` (int; default=1000) Number of permutations for null distribution estimation
11. `--q_thres` (float; default=0.1) FDR threshold for significance

#### output

1. `outdir/celltype__trait-*.csv`
2. `outdir/dfbs__trait-*.npz`
3. `outdir/metacell__trait-*.csv`

where `*` is trait name we assume magmaz file name is `*.genes.out`
