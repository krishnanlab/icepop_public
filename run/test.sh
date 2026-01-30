poetry run icepop metacell \
--h5ad ../data/TM_FACS_cnt.h5ad \
--outdir ../data/test_run \
--save_name TM_FACS_cnt

poetry run icepop association \
--h5ad ../data/TM_FACS_cnt.h5ad \
--mc_assign ../results/test_run/mc_assign.csv \
--magmaz ../data/magmaz/IPF.genes.out \
--spec_score ../results/test_run/mc_spec_score.npz \
--outdir ../results/test_run
