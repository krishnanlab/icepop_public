from typing import Tuple
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
from multiprocessing.shared_memory import SharedMemory
from scipy.sparse import csr_matrix
from time import time
from icepop.logging_config import logger
from scipy.stats import norm


# ============================================================
# Core linear regression
# ============================================================

def _linear_reg(
    X: np.ndarray,
    y: np.ndarray,
    eps: float = 1e-12,
    dfb: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Perform OLS regression of y on each row of X (implicit intercept).

    Parameters
    ----------
    X : np.ndarray
        Shape (n_metacell, n_gene)
    y : np.ndarray
        Shape (n_gene,)
    eps : float
        Numerical stability constant
    dfb : bool
        Whether to compute DFBETAS

    Returns
    -------
    beta : np.ndarray
        Regression coefficients (n_metacell,)
    se : np.ndarray
        Standard errors (n_metacell,)
    dfb_mat : np.ndarray or None
        DFBETAS (n_metacell, n_gene) if dfb=True
    """
    X = np.nan_to_num(X, nan=0.0)
    y = np.nan_to_num(y, nan=0.0)

    if X.ndim == 1:
        X = X[None, :]

    y = y - y.mean()
    X = X - X.mean(axis=1, keepdims=True)

    xx = np.sum(X ** 2, axis=1) + eps
    beta = (X @ y) / xx

    resid = y[None, :] - beta[:, None] * X
    rss = np.sum(resid ** 2, axis=1)
    s2 = rss / (X.shape[1] - 1)
    se = np.sqrt(s2 / xx)

    if not dfb:
        return beta, se, None

    h = (X ** 2) / (xx[:, None] + eps)
    dfb_mat = ((X / (xx[:, None] + eps)) * resid) / (1 - h + eps)
    return beta, se, dfb_mat


# ============================================================
# Multiprocessing worker (top-level required)
# ============================================================
def _lr_util(args):
    shm_data_name, data_dtype, data_shape, y, eps = args

    shm_data = SharedMemory(name=shm_data_name)
    X = np.ndarray(data_shape, dtype=data_dtype, buffer=shm_data.buf)

    X = np.nan_to_num(X, nan=0.0)
    y = np.nan_to_num(y, nan=0.0)

    if X.ndim == 1:
        X = X[None, :]

    k, n = X.shape

    # mean-center (implicit intercept)
    y = y - y.mean()
    X = X - X.mean(axis=1, keepdims=True)

    xx = np.sum(X**2, axis=1) + eps
    beta = (X @ y) / xx

    resid = y[None, :] - beta[:, None] * X
    rss = np.sum(resid**2, axis=1)
    s2 = rss / (n - 1)
    se = np.sqrt(s2 / xx)

    return beta, se


# ============================================================
# 2. Permutation-estimated covariance (core)
# ============================================================
def _run_parallel_lr(
    X, y,
    n_perm=1000,
    random_state=42,
    n_jobs=20,
    eps=1e-12
):
    """
    run linear reg for permed y in parallel
    """
    rng = np.random.default_rng(random_state)

    if n_jobs > 1:
        # Create shared memory blocks
        shm_data = SharedMemory(create=True, size=X.nbytes)

        # Copy data to shared memory
        np_data_shm = np.ndarray(X.shape, dtype=X.dtype, buffer=shm_data.buf)
        np_data_shm[:] = X[:]

        # prepare permutation bank
        y_perms = [rng.permutation(y) for _ in range(n_perm)]

        # Prepare task parameters
        tasks = [
            (
                shm_data.name, X.dtype, X.shape,
                y_perms[i],
                eps
            )
            for i in range(n_perm)
        ]

        logger.info("[pool] launching workers...")
        t0 = time()
        try:
            with Pool(n_jobs) as pool:
                res = pool.map(_lr_util, tasks)
        finally:
            # Cleanup shared memory in main process
            shm_data.close()
            shm_data.unlink()
        logger.info(f"[pool] finished in {(time() - t0) / 60:.2f} min")

        beta_perm, se_perm = zip(*res)
        beta_perm = np.vstack(beta_perm).astype(np.float32)
        se_perm = np.vstack(se_perm).astype(np.float32)
        return beta_perm, se_perm
    else:
        logger.info("[serial] running without multiprocessing")
        t0 = time()
        beta_perm = np.zeros((n_perm, X.shape[0]), dtype=np.float32)
        se_perm = np.zeros((n_perm, X.shape[0]), dtype=np.float32)
        for i in range(n_perm):
            beta_perm[i], se_perm[i], _ = _linear_reg(X, rng.permutation(y), eps=eps)
        logger.info(f"[serial] finished in {(time() - t0) / 60:.2f} min")
        return beta_perm, se_perm


# ============================================================
# 4. Cell-type aggregation using covariance
# ============================================================

def _celltype_from_metacell(beta_hat, cov_beta, f, null):
    """
    Aggregate metacell effects to cell-type level.

    Parameters
    ----------
    f : (n_metacell,)  weights (e.g. metacell proportions)

    Returns
    -------
    beta_ct, se_ct, z_ct, p_ct
    """
    beta_ct = f @ beta_hat
    var_ct = f @ cov_beta @ f
    se_ct = np.sqrt(var_ct + 1e-12)
    z_ct = beta_ct / se_ct

    mu, sigma = norm.fit(null)
    if null.sum() == 0:
        mu = 0
        sigma = 1e-12
    p_ct = norm.sf(z_ct, loc=mu, scale=sigma)

    return beta_ct, se_ct, z_ct, p_ct


# ============================================================
# Main association model
# ============================================================

class MetacellAssoc:
    """
    Metacell → cell-type association model using permutation GLS.

    Notes
    -----
    - sep and comp must be precomputed
    - DFBETAS are returned for post-hoc analysis
    """

    def __init__(
        self,
        n_perm: int = 1000,
        n_jobs: int = 20,
        eps: float = 1e-12,
        random_state: int = 42,
        q_thres: float = 0.1,
        ct_key: str = 'cell_type'
    ):
        self.n_perm = n_perm
        self.n_jobs = n_jobs
        self.eps = eps
        self.random_state = random_state
        self.q_thres = q_thres
        self.ct_key = ct_key

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        freq_df: pd.DataFrame,
        c2mc2ct: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run metacell → cell-type association.

        Parameters
        ----------
        X : np.ndarray
            Metacell × gene matrix
        y : np.ndarray
            Gene-level score
        freq_df : pd.DataFrame
            Celltype × metacell weight matrix
        comp : np.ndarray
            Metacell compactness score
        sep : np.ndarray
            Metacell specificity score

        Returns
        -------
        ct_df : pd.DataFrame
            Cell-type association results
        mc_df : pd.DataFrame
            Metacell-level results
        """
        X = X.astype(np.float32)
        y = y.astype(np.float32)

        # --------------------------------
        # Metacell-level association
        # --------------------------------
        # get mc beta and se
        beta_hat, se_hat, dfb = _linear_reg(X, y, eps=self.eps, dfb=True)

        # estimate cov of metacell beta
        beta_perm, se_perm = _run_parallel_lr(
            X, y,
            n_perm=self.n_perm,
            n_jobs=self.n_jobs,
            eps=self.eps,
            random_state=self.random_state
        )
        cov_beta = np.cov(beta_perm, rowvar=False, dtype=np.float32)

        z = beta_hat / se_hat
        p = norm.sf(z)

        mc_df = pd.DataFrame({
            "metacell": freq_df.columns,
            "beta": beta_hat,
            "se": se_hat,
            "z": z,
            "p": p,
        })

        # --------------------------------
        # Build metacell weights per cell type
        # --------------------------------
        sig_w = norm.cdf(z)

        metacell_weight = {}
        for ct, f in freq_df.iterrows():
            w = sig_w * np.asarray(f, dtype=np.float32)
            tot = w.sum()
            metacell_weight[ct] = w / (tot if tot > 0 else 1e-12)

        # --------------------------------
        # null t distribution for each cell type
        # --------------------------------
        t0 = time()
        logger.info("[perm-null] start building null distributions of association")
        metacell_perm_t = {}
        sig_w_perm = norm.cdf(beta_perm / se_perm)
        for celltype, f in freq_df.iterrows():
            f = np.asarray(f, dtype=np.float32)

            # perm ct beta
            f_perm = sig_w_perm * f[None, :]
            tot = f_perm.sum(1)
            tot = np.where(tot == 0, 1e-12, tot)
            f_perm = f_perm / tot[:, None]
            beta_ct_perm = (beta_perm * f_perm).sum(1)

            # perm ct se
            f_perm = csr_matrix(f_perm)
            out = np.empty(f_perm.shape[0], dtype=np.float32)
            for i in range(f_perm.shape[0]):
                idx = f_perm.indices[f_perm.indptr[i]:f_perm.indptr[i+1]]
                val = f_perm.data[f_perm.indptr[i]:f_perm.indptr[i+1]]
                out[i] = val @ cov_beta[np.ix_(idx, idx)] @ val
            se_ct_perm = np.sqrt(out + 1e-12)

            # get ct null distribution
            t_ct_perm = beta_ct_perm / se_ct_perm
            metacell_perm_t[celltype] = t_ct_perm
        logger.info(f"[perm-null] finished in {(time() - t0) / 60:.2f} min")

        # --------------------------------
        # Aggregate to cell-type level
        # --------------------------------
        t0 = time()
        logger.info(f"[Aggregate] Aggregate metacell association to {self.ct_key} level")
        ct_res = []
        for ct, f in metacell_weight.items():
            beta_ct, se_ct, z_ct, p_ct = _celltype_from_metacell(
                beta_hat, cov_beta, f, metacell_perm_t[ct]
            )
            ct_res.append((ct, beta_ct, se_ct, z_ct, p_ct))
        ct_df = pd.DataFrame(
            ct_res,
            columns=[self.ct_key, "beta", "se", "z", "p"]
        )

        # multiple testing on cell type
        ct_df['q'] = multipletests(ct_df['p'], method="fdr_bh")[1]
        logger.info(f"[mixture] finished in {(time() - t0) / 60:.2f} min")

        # --------------------------------
        # Estimate mixture of association
        # --------------------------------
        t0 = time()
        logger.info(f"[mixture] Estimate significance of metacells within {self.ct_key}")
        mc_p = mc_df['p'].values
        metacells = freq_df.columns.values
        sig_pct_list = []
        for ct, q in zip(ct_df[self.ct_key], ct_df['q']):
            sig_pct = 0.0
            if q <= self.q_thres:
                purity = freq_df.loc[ct, :].values
                if (purity > 0).sum() >= 2:
                    p = mc_p[purity > 0]
                    pos_mcs = metacells[np.argwhere(purity > 0).ravel()]
                    w = purity[purity > 0]
                    w = w / w.sum() * w.size
                    q = multipletests(p / w, method='fdr_bh')[1]
                    mc2q = dict(zip(pos_mcs, q))
                    c_q = np.asarray([
                        mc2q[mc] if mc in mc2q else 1.0
                        for mc in c2mc2ct[c2mc2ct[self.ct_key] == ct]['metacell']
                    ])
                    sig_pct = np.mean(c_q <= self.q_thres)
            sig_pct_list.append(sig_pct)
        ct_df['sig_pct'] = sig_pct_list
        logger.info(f"[mixture] finished in {(time() - t0) / 60:.2f} min")

        # --------------------------------
        # Calculate influence per cell type
        # --------------------------------
        t0 = time()
        logger.info(f"[influence] Get influence score for given {self.ct_key}")
        ctdfbs = np.vstack([i for i in metacell_weight.values()]) @ dfb / ct_df['se'].values[:, None]
        logger.info(f"[influence] finished in {(time() - t0) / 60:.2f} min")

        return ct_df.sort_values("q", ignore_index=True), mc_df, ctdfbs
