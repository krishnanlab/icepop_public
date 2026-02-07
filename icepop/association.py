import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
from icepop.specificity_score import specificity_score
from icepop.convert_score import CrossSpeciesScoreConverter
from icepop.model import MetacellAssoc
import logging
from time import time

logger = logging.getLogger("icepop.association")


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
    q_thres: float = 0.1
):
    """
    Run metacell-based gene–trait association analysis.

    This command integrates single-cell metacell specificity scores
    with MAGMA gene-level statistics to infer disease- or trait-associated
    cell types.

    Parameters
    ----------
    h5ad : str
        Input AnnData file containing single-cell expression data.
    mc_assign : str
        CSV file mapping cells to metacell assignments.
    magmaz : str
        MAGMA gene-level association file (*.genes.out).
    outdir : str
        Output directory for association results.
    spec_score : str, optional
        Precomputed metacell specificity score file (.npz).
        If not provided, specificity scores are computed.
    n_jobs : int, default=20
        Number of parallel workers.
    sp : str, default='mmusculus'
        Species identifier for gene ID conversion.
    ct_key : str, default='cell_type'
        Column in `adata.obs` defining cell types.
    trait_name : str, optional
        Trait name used for output file naming.
    n_perm : int, default=1000
        Number of permutations for null distribution estimation.
    q_thres : float, default=0.1
        FDR threshold for significance.

    Outputs
    -------
    - metacell-level association results (CSV)
    - cell-type–level association results (CSV)
    - influence diagnostics (NPZ)
    """

    t0 = time()
    logger.info("Starting association analysis")
    logger.info(
        f"Inputs: h5ad={h5ad}, mc_assign={mc_assign}, magmaz={magmaz}"
    )
    logger.info(
        f"Params: ct_key={ct_key}, sp={sp}, n_perm={n_perm}, "
        f"q_thres={q_thres}, n_jobs={n_jobs}"
    )

    # load data and mc assignment
    logger.info("Loading h5ad and metacell assignment")
    adata = sc.read_h5ad(h5ad)
    logger.info(f"Loaded adata: n_cells={adata.n_obs}, n_genes={adata.n_vars}")
    logger.info("Number of %s: %d" % (ct_key, adata.obs[ct_key].unique().size))

    adata.obs['metacell'] = pd.read_csv(mc_assign, header=None)[0].values
    logger.info("Number of metacells %d" % (adata.obs['metacell'].unique().size))
    logger.info(
        f"Metacell sizes: "
        f"min={adata.obs['metacell'].value_counts().min()}, "
        f"median={adata.obs['metacell'].value_counts().median()}, "
        f"max={adata.obs['metacell'].value_counts().max()}"
    )

    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # mkdir
    Path(outdir).mkdir(exist_ok=True)

    # check if metacell columns in adata
    # get spec score
    if spec_score is not None:
        logger.info("Loading precomputed metacell specificity scores")
        f = np.load(spec_score, allow_pickle=True)
        adata.uns['spec_score'] = pd.DataFrame(f['score'], index=f['mc'], columns=f['genes'])
    else:
        logger.info("Computing metacell specificity scores")
        t_spec = time()
        spec = specificity_score(adata, n_jobs=n_jobs)
        spec.get_metacell_spec_score()
        np.savez_compressed(
            f'{outdir}/mc_spec_score.npz',
            score=np.asarray(adata.uns['spec_score'], np.float32),
            mc=adata.uns['spec_score'].index.values,
            genes=adata.uns['spec_score'].columns.values,
        )
        logger.info(f"Finished in {(time() - t_spec) / 60:.2f} min")

    # translate ID
    logger.info(f"Converting specificity scores across species ({sp})")
    score_converter = CrossSpeciesScoreConverter(adata, sp=sp)
    mc_spec_score = score_converter.convert_score_across_species(adata.uns['spec_score'], normed=False)
    if mc_spec_score.shape[1] == 0:
        logger.error(
            "Converted spec-score matrix is empty. "
            "shape=%s. "
            "Likely causes: gene ID mismatch, wrong spec_score input",
            mc_spec_score.shape,
        )
        raise ValueError("Empty spec-score matrix after conversion")
    metacells = np.asarray(mc_spec_score.index)
    logger.info(f"Converted score matrix: {mc_spec_score.shape}")

    # get metacell freq
    freq_df = pd.crosstab(adata.obs[ct_key], adata.obs['metacell'])
    freq_df = freq_df.loc[:, metacells]
    freq_df = freq_df.div(freq_df.sum(0))

    # check stats
    purity = freq_df.max(axis=0)
    logger.info(
        f"Metacell purity: "
        f"{(purity >= 0.2).sum()}/{len(purity)} pass min_purity=0.2"
    )

    # output name
    if trait_name is None:
        trait_name = Path(magmaz).name.replace('.genes.out', '')
    metacell_out = f'{outdir}/metacell__trait-{trait_name}.csv'
    celltype_out = f'{outdir}/celltype__trait-{trait_name}.csv'
    dfbs_out = f'{outdir}/dfbs__trait-{trait_name}.npz'

    # align gene set between magma z and sc data
    magmaz_df = pd.read_csv(magmaz, header=0, index_col=0, sep=r'\s+')
    magmaz_df.index = magmaz_df.index.astype(str)
    shared_genes = magmaz_df.index.intersection(mc_spec_score.columns)
    y = np.asarray(magmaz_df.loc[shared_genes, 'ZSTAT'])
    X = mc_spec_score.loc[:, shared_genes].to_numpy()
    logger.info(f"Shared genes between MAGMA and scRNA: {len(shared_genes)}")

    t_fit = time()
    logger.info(
        f"Running MetacellAssoc "
        f"(n_perm={n_perm}, n_jobs={n_jobs})"
    )
    # model fitting
    assoc = MetacellAssoc(
        n_perm=n_perm, n_jobs=n_jobs,
        q_thres=q_thres, ct_key=ct_key
    )
    ct_df, mc_df, ctdfbs = assoc.fit(
        X, y, freq_df,
        adata.obs.loc[:, ['metacell', ct_key]].copy()
    )
    logger.info(f"Finished in {(time() - t_fit) / 60:.2f} min")

    # save output
    mc_df.to_csv(metacell_out, header=True, index=False)
    ct_df.to_csv(celltype_out, header=True, index=False)
    np.savez_compressed(
        dfbs_out,
        dfbs=ctdfbs,
        celltypes=freq_df.index.values,
        genes=shared_genes
    )
    logger.info(f"Assoication analysis finished in {(time() - t0) / 60:.2f} min")
