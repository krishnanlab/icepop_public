from pathlib import Path
import scanpy as sc
import pandas as pd
import SEACells
from MetaQ_sc import run_metaq
from icepop.logging_config import logger


def celltype_frac(x, col_name):
    val_counts = x[col_name].value_counts()
    return val_counts.values[0] / val_counts.values.sum()


def calc_stats(
    adata,
    build_kernel_on='X_pca', ct_key='cell_type', SEACells_label='metacell',
    n_top_genes=2000, n_comp=50
):
    if build_kernel_on not in adata.obsm.keys():
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        sc.tl.pca(adata, n_comps=n_comp, use_highly_variable=True)

    # compactness, normalized and reverted
    compactness = SEACells.evaluate.compactness(
        adata, low_dim_embedding=build_kernel_on, SEACells_label=SEACells_label
    )
    # separation, normalized and reverted
    separation = SEACells.evaluate.separation(
        adata,
        low_dim_embedding=build_kernel_on, nth_nbr=1, SEACells_label=SEACells_label
    )
    # cell type purity
    celltype_fraction = adata.obs.groupby(SEACells_label).apply(lambda x: celltype_frac(x, ct_key))
    celltype = adata.obs.groupby(SEACells_label).apply(lambda x: x[ct_key].value_counts().index[0])
    ct_purity = pd.concat([celltype, celltype_fraction], axis=1).rename(columns={0: ct_key, 1: 'purity'})
    return pd.concat([compactness, separation, ct_purity], axis=1)


def metacell(
    h5ad: str,
    outdir: str,
    ncell_per_mc: int = 75,
    save_name: str = 'metaq_res',
    ct_key: str = 'cell_type',
    device: str = 'cuda',
    report_stat: bool = True
):
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

    # add metaq
    metacellid_adata = sc.read(metacellids_h5ad)

    # add metacell into adata
    adata.obs['metacell'] = [f'metacell-{i}' for i in metacellid_adata.obs['metacell']]

    # save metacell assignment
    with open(f'{outdir}/mc_assign.csv', 'w') as f:
        for i in metacellid_adata.obs['metacell']:
            f.write(f'metacell-{i}\n')

    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # summary stats using SEACell
    stats_df = calc_stats(adata)
    stats_df.to_csv(f'{outdir}/mc_stats.csv', header=True, index=True)

    # report mc stats
    if report_stat:
        # report stats
        logger.info('Compactness stats:')
        desc = stats_df['compactness'].describe()
        logger.info(
            "count=%d mean=%.4f std=%.4f min=%.4f "
            "25%%=%.4f median=%.4f 75%%=%.4f max=%.4f",
            desc['count'], desc['mean'], desc['std'], desc['min'],
            desc['25%'], desc['50%'], desc['75%'], desc['max']
        )

        logger.info('Separation stats:')
        desc = stats_df['separation'].describe()
        logger.info(
            "count=%d mean=%.4f std=%.4f min=%.4f "
            "25%%=%.4f median=%.4f 75%%=%.4f max=%.4f",
            desc['count'], desc['mean'], desc['std'], desc['min'],
            desc['25%'], desc['50%'], desc['75%'], desc['max']
        )

        logger.info('Cell type purity stats:')
        desc = stats_df['purity'].describe()
        logger.info(
            "count=%d mean=%.4f std=%.4f min=%.4f "
            "25%%=%.4f median=%.4f 75%%=%.4f max=%.4f",
            desc['count'], desc['mean'], desc['std'], desc['min'],
            desc['25%'], desc['50%'], desc['75%'], desc['max']
        )
