from pathlib import Path
import scanpy as sc
import pandas as pd
from MetaQ_sc import run_metaq
import logging

logger = logging.getLogger("icepop.metacell")


def metacell(
    h5ad: str,
    outdir: str,
    ncell_per_mc: int = 75,
    save_name: str = 'metaq_res',
    ct_key: str = 'cell_type',
    device: str = 'cuda',
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
    min_purity : float, default=0.2
        Minimum metacell purity threshold.
    q_thres : float, default=0.1
        FDR threshold for significance.

    Outputs
    -------
    - metacell-level association results (CSV)
    - cell-type–level association results (CSV)
    - influence diagnostics (NPZ)
    """

    adata = sc.read_h5ad(h5ad)

    # mkdir
    Path(outdir).mkdir(exist_ok=True)

    # check ct key
    if ct_key not in adata.obs:
        raise KeyError(f"Cell-type key '{ct_key}' not found in adata.obs")

    # check count/log normed
    max_val = adata.X.max()
    logger.info(f"Max expression value: {max_val}")

    if max_val <= 15:
        logger.warning(
            "Max expression value (%s) is low; verify that input is raw counts.",
            max_val,
        )

    # check stats
    tot_n_cell = adata.shape[0]
    target_num = int(round(tot_n_cell / ncell_per_mc))
    logger.info(f'Totol number of cells: {tot_n_cell}, target metacell number: {target_num}')

    # metacell id file
    metacellids_h5ad = f'./save/{save_name}_{target_num}metacell_ids.h5ad'

    # run metaq
    if not Path(metacellids_h5ad).exists():
        logger.info("Running MetaQ")
        try:
            run_metaq(
                data_path=[
                    h5ad,
                ],  # the path to the input h5ad data
                data_type=[
                    "RNA",
                ],  # the type of the input data
                metacell_num=target_num,  # the target number of metacells
                save_name=save_name,  # the file name prefix when saving the results
                type_key=ct_key,
                device=device
            )
        except KeyError as e:
            logger.warning(
                f"[INFO] MetaQ evaluation failed (celltype mapping), ignored: {e}"
            )

    # save metacell assignment
    metacellid_adata = sc.read(metacellids_h5ad)
    with open(f'{outdir}/mc_assign.csv', 'w') as f:
        for i in metacellid_adata.obs['metacell']:
            f.write(f'metacell-{i}\n')

    # purity stats
    adata.obs['metacell'] = [f'metacell-{i}' for i in metacellid_adata.obs['metacell']]
    freq_df = pd.crosstab(adata.obs[ct_key], adata.obs['metacell'])
    freq_df = freq_df.div(freq_df.sum(0))

    mc2ct_df = pd.DataFrame({
        "cell_type": freq_df.idxmax(axis=0),
        "purity": freq_df.max(axis=0),
    })
    mc2ct_df.to_csv(f'{outdir}/mc_stats.csv', header=True, index=True)
