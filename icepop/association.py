import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
from icepop.specificity_score import specificity_score
from icepop.convert_score import CrossSpeciesScoreConverter
from icepop.model import MetacellAssoc


def association(
    h5ad: str,
    mc_assign: str,
    magmaz: str,
    outdir: str,
    spec_score: str = None,
    n_jobs: int = 20,
    sp: str = 'mmusculus',
    ct_key: str = 'cell_type',
    trait_name: str = None,
    n_perm: int = 1000,
    min_purity: float = 0.2,
    q_thres: float = 0.1
):
    # load data and mc assignment
    adata = sc.read_h5ad(h5ad)
    adata.obs['metacell'] = pd.read_csv(mc_assign, header=None)[0].values

    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # mkdir
    Path(outdir).mkdir(exist_ok=True)

    # check if metacell columns in adata
    # get spec score
    if spec_score is not None:
        f = np.load(spec_score, allow_pickle=True)
        adata.uns['spec_score'] = pd.DataFrame(f['score'], index=f['mc'], columns=f['genes'])
    else:
        spec = specificity_score(adata, n_jobs=n_jobs)
        spec.get_metacell_spec_score()
        np.savez_compressed(
            f'{outdir}/mc_spec_score.npz',
            score=np.asarray(adata.uns['spec_score'], np.float32),
            mc=adata.uns['spec_score'].index.values,
            genes=adata.uns['spec_score'].columns.values,
        )

    # translate ID
    score_converter = CrossSpeciesScoreConverter(adata, sp=sp)
    mc_spec_score = score_converter.convert_score_across_species(adata.uns['spec_score'], normed=False)
    metacells = np.asarray(mc_spec_score.index)

    # get metacell freq
    freq_df = pd.crosstab(adata.obs[ct_key], adata.obs['metacell'])
    freq_df = freq_df.loc[:, metacells]
    freq_df = freq_df.div(freq_df.sum(0))

    # output name
    if trait_name is None:
        trait_name = Path(magmaz).name.replace('.genes.out', '')
    metacell_res = f'{outdir}/metacell__trait-{trait_name}.csv'
    celltype_res = f'{outdir}/celltype__trait-{trait_name}.csv'

    # align gene set between magma z and sc data
    magmaz_df = pd.read_csv(magmaz, header=0, index_col=0, sep=r'\s+')
    magmaz_df.index = magmaz_df.index.astype(str)
    shared_genes = magmaz_df.index.intersection(mc_spec_score.columns)
    y = np.asarray(magmaz_df.loc[shared_genes, 'ZSTAT'])
    X = mc_spec_score.loc[:, shared_genes].to_numpy()

    # model fitting
    assoc = MetacellAssoc(
        n_perm=n_perm, n_jobs=n_jobs,
        min_purity=min_purity, q_thres=q_thres,
        ct_key=ct_key
    )
    ct_df, mc_df = assoc.fit(
        X, y, freq_df,
        adata.obs.loc[:, ['metacell', ct_key]].copy()
    )

    # save output
    mc_df.to_csv(metacell_res, header=True, index=False)
    ct_df.to_csv(celltype_res, header=True, index=False)
