# icepop_public
public repo for icepop

# extract metacells
poetry run icepop metacell \
    --h5ad ../data/mouse_colon/mouse_colon_cnt.h5ad \
    --outdir ../results/mouse_colon_mc \
    --save_name mouse_colon

## input options
--h5ad count in .X of h5ad
--outdir outdir
--save_name prefix of metaq output under ./save/*, do not write a path

this step need gpu for faster speed

## output
metacell assignment

# get association, mixture and influence diagnoistics
poetry run icepop association \
    --h5ad ../data/TM_FACS/TM_FACS_cnt.h5ad \
    --mc_assign ../results/TM_FACS_mc/mc_assign.csv \
    --magmaz ../data/TM_FACS/magmaz/asd.genes.out \
    --spec_score ../results/TM_FACS_mc/mc_spec_score.npz \
    --outdir ./test

## input options
--h5ad count in .X of h5ad
--mc_assign metacell assignment from last step
--magmaz magmaz sum stat
--spec_score spec score of metacell, will auto generated in the first run
--outdir outdir

## output
outdir/celltype__trait-*.csv
outdir/dfbs__trait-*.npz
outdir/metacell__trait-*.csv

* is trait name we assume magmaz file name is *.genes.out
