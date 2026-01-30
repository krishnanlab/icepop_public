from pathlib import Path
import scanpy as sc
from MetaQ_sc import run_metaq
from icepop.logging_config import logger


def celltype_frac(x, col_name):
    val_counts = x[col_name].value_counts()
    return val_counts.values[0] / val_counts.values.sum()


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
