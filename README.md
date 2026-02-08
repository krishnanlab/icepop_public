# ICePop: Informative Cell Population
This repository contains source code for ICePop ([DOI](TBD))

# Installation
## Dependencies
- ICePop dependencies are handled by [poetry](https://python-poetry.org/) with `python>=3.11,<3.12`
To install poetry, please follow the [instructions on poetry's home page](https://python-poetry.org/docs/#installation).
- Install dependencies
```
poetry install poetry.lock
```

## Install ICePop
```
poetry add --editable /path/to/icepop
```

# Run ICePop
## Step 1: Extract metacells
```
poetry run icepop metacell \
    --h5ad ../data/mouse_colon/mouse_colon_cnt.h5ad \
    --outdir ../results/mouse_colon_mc \
    --save_name mouse_colon
```

### input options

1. `--h5ad` single-cell expression count in .X of h5ad
2. `--outdir` output directory path
3. `--save_name` prefix of metaq output under `./save/*`, do not write a path

this step need gpu for faster speed

### output
metacell assignment

## Step 2: Get association, mixture and influence diagnoistics
```
poetry run icepop association \
    --h5ad ../data/TM_FACS/TM_FACS_cnt.h5ad \
    --mc_assign ../results/TM_FACS_mc/mc_assign.csv \
    --magmaz ../data/TM_FACS/magmaz/asd.genes.out \
    --spec_score ../results/TM_FACS_mc/mc_spec_score.npz \
    --outdir ./test
```

### input options

1. `--h5ad` single-cell expression count in .X of h5ad
2. `--mc_assign` metacell assignment output from last step
3. `--magmaz` [magmaz](https://doi.org/10.1371/journal.pcbi.1004219) summary statistics
4. `--spec_score` specificity score of metacell, will auto generated in Step 1
5. `--outdir` output directory path

### output

1. `outdir/celltype__trait-*.csv`
2. `outdir/dfbs__trait-*.npz`
3. `outdir/metacell__trait-*.csv`

where `*` is trait name we assume magmaz file name is `*.genes.out`
