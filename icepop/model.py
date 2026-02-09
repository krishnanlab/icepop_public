from typing import Tuple
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
from multiprocessing.shared_memory import SharedMemory
from scipy.sparse import csr_matrix
from time import time
from scipy.stats import norm
import logging

logger = logging.getLogger(__name__)


def _linear_reg(
    X: np.ndarray,
    y: np.ndarray,
    eps: float = 1e-12,
    dfb: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Perform row-wise ordinary least squares (OLS) regression of ``y`` on each
    row of ``X`` using an implicit intercept (mean-centering).

    Each row of ``X`` is treated as an independent predictor vector for the
    same response vector ``y``. This is equivalent to fitting multiple simple
    linear regressions:

        y = beta_i * X_i + error

    after subtracting row-wise means from ``X`` and the mean from ``y``.

    Parameters
    ----------
    X : np.ndarray of shape (n_metacell, n_gene)
        Predictor matrix where each row represents a metacell profile across genes.
    y : np.ndarray of shape (n_gene,)
        Gene-level response vector.
    eps : float, default=1e-12
        Small numerical constant added to denominators for stability.
    dfb : bool, default=False
        Whether to compute DFBETAS influence statistics for each observation.

    Returns
    -------
    beta : np.ndarray of shape (n_metacell,)
        Estimated regression coefficient for each row of ``X``.
    se : np.ndarray of shape (n_metacell,)
        Standard error of each regression coefficient.
    dfb_mat : np.ndarray of shape (n_metacell, n_gene) or None
        DFBETAS influence values if ``dfb=True``; otherwise ``None``.
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


def _lr_util(args):
    """
    Worker utility for permutation-based row-wise OLS regression using shared memory.

    This function is designed to be executed inside a multiprocessing pool.
    It reconstructs the predictor matrix ``X`` from a shared memory buffer,
    performs mean-centered simple linear regression of a (possibly permuted)
    response vector ``y`` on each row of ``X``, and returns regression
    coefficients and standard errors.

    The regression model for each row ``i`` is:

        y = beta_i * X_i + error

    after subtracting the mean from both ``y`` and each row of ``X``
    (equivalent to fitting an intercept).

    Parameters
    ----------
    args : tuple
        Packed arguments required for multiprocessing execution:

        - shm_data_name : str  
          Name of the shared memory block containing ``X``.
        - data_dtype : numpy.dtype  
          Data type of the shared array.
        - data_shape : tuple of int  
          Shape of ``X`` as ``(n_metacell, n_gene)``.
        - y : np.ndarray of shape (n_gene,)  
          Response vector (typically a permutation of the original response).
        - eps : float  
          Small numerical stability constant added to denominators.

    Returns
    -------
    beta : np.ndarray of shape (n_metacell,)
        Estimated regression coefficient for each row of ``X``.
    se : np.ndarray of shape (n_metacell,)
        Standard error of each regression coefficient.

    Notes
    -----
    - Uses shared memory to avoid duplicating large arrays across processes.
    - Assumes rows of ``X`` are independent predictors evaluated against
      the same response vector ``y``.
    - Does **not** compute influence statistics (e.g., DFBETAS); see
      :func:`_linear_reg` for the full regression implementation.
    """
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


def _run_parallel_lr(
    X, y,
    n_perm=1000,
    random_state=42,
    n_jobs=20,
    eps=1e-12
):
    """
    Compute permutation-based linear regressions in parallel.

    For each permutation of ``y``, performs row-wise OLS regression of the
    permuted response on ``X``. Shared memory is used to avoid duplicating
    large arrays across worker processes.

    Parameters
    ----------
    X : np.ndarray of shape (n_metacell, n_gene)
        Predictor matrix.
    y : np.ndarray of shape (n_gene,)
        Response vector to be permuted.
    n_perm : int, default=1000
        Number of random permutations of ``y``.
    random_state : int, default=42
        Seed for reproducible permutations.
    n_jobs : int, default=20
        Number of parallel worker processes. If ``1``, runs serially.
    eps : float, default=1e-12
        Numerical stability constant used in regression.

    Returns
    -------
    beta_perm : np.ndarray of shape (n_perm, n_metacell)
        Regression coefficients from each permutation.
    se_perm : np.ndarray of shape (n_perm, n_metacell)
        Standard errors from each permutation.
    """
    rng = np.random.default_rng(random_state)

    if n_jobs > 1:
        # create shared memory blocks
        shm_data = SharedMemory(create=True, size=X.nbytes)

        # copy data to shared memory
        np_data_shm = np.ndarray(X.shape, dtype=X.dtype, buffer=shm_data.buf)
        np_data_shm[:] = X[:]

        # prepare permutation bank
        y_perms = [rng.permutation(y) for _ in range(n_perm)]

        # prepare task parameters
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
            # clean up shared memory in main process
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


def _celltype_from_metacell(beta_hat, cov_beta, f, null):
    """
    Aggregate metacell-level regression effects to a cell-type–level statistic.

    The cell-type effect is computed as a weighted linear combination of
    metacell coefficients:

        beta_ct = f^T beta_hat
        var_ct  = f^T cov_beta f

    A z-score is derived and evaluated against an empirical null distribution
    estimated from permutation statistics.

    Parameters
    ----------
    beta_hat : np.ndarray of shape (n_metacell,)
        Observed metacell regression coefficients.
    cov_beta : np.ndarray of shape (n_metacell, n_metacell)
        Covariance matrix of metacell coefficients estimated from permutations.
    f : np.ndarray of shape (n_metacell,)
        Normalized metacell weights for the given cell type.
    null : np.ndarray of shape (n_perm,)
        Permutation-derived null distribution of cell-type z-scores.

    Returns
    -------
    beta_ct : float
        Aggregated cell-type regression coefficient.
    se_ct : float
        Standard error of the aggregated coefficient.
    z_ct : float
        Z-score of the cell-type association.
    p_ct : float
        One-sided p-value computed from the fitted normal null distribution.
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


class MetacellAssoc:
    """
    Permutation-based metacell → cell-type association model.

    This class estimates gene-level associations at the metacell level using
    row-wise OLS regression, derives covariance from permutation testing,
    and aggregates effects to the cell-type level using generalized least
    squares (GLS)-style weighting.

    The workflow includes:

    1. Metacell-level regression and standard errors
    2. Permutation-derived covariance estimation
    3. Cell-type aggregation with empirical null calibration
    4. Multiple-testing correction across cell types
    5. Within–cell-type metacell significance mixture estimation
    6. Influence scoring via DFBETAS

    Parameters
    ----------
    n_perm : int, default=1000
        Number of permutations used to estimate null distributions.
    n_jobs : int, default=20
        Number of parallel worker processes.
    eps : float, default=1e-12
        Numerical stability constant.
    random_state : int, default=42
        Seed for permutation reproducibility.
    q_thres : float, default=0.1
        False discovery rate (FDR) threshold used for significance summaries.
    ct_key : str, default="cell_type"
        Column name used to label cell types in output tables.
    output_dfbs: bool: default=True
        Whether output dfbetas
    """

    def __init__(
        self,
        n_perm: int = 1000,
        n_jobs: int = 20,
        eps: float = 1e-12,
        random_state: int = 42,
        q_thres: float = 0.1,
        ct_key: str = 'cell_type',
        output_dfbs: bool = True
    ):
        self.n_perm = n_perm
        self.n_jobs = n_jobs
        self.eps = eps
        self.random_state = random_state
        self.q_thres = q_thres
        self.ct_key = ct_key
        self.output_dfbs = output_dfbs

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        freq_df: pd.DataFrame,
        c2mc2ct: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run the full metacell → cell-type association pipeline.

        Parameters
        ----------
        X : np.ndarray of shape (n_metacell, n_gene)
            Metacell-by-gene expression or feature matrix.
        y : np.ndarray of shape (n_gene,)
            Gene-level score or phenotype association vector.
        freq_df : pd.DataFrame of shape (n_celltype, n_metacell)
            Cell-type–specific metacell weights (e.g., proportions or purity).
        c2mc2ct : pd.DataFrame
            Mapping table linking cells, metacells, and cell types. Must contain
            columns ``["metacell", ct_key]``.

        Returns
        -------
        ct_df : pd.DataFrame
            Cell-type–level association statistics including beta, standard error,
            z-score, p-value, FDR q-value, and significant metacell percentage.
        mc_df : pd.DataFrame
            Metacell-level regression statistics.
        ctdfbs : np.ndarray of shape (n_celltype, n_gene)
            Cell-type–level influence scores derived from DFBETAS.
        """
        X = X.astype(np.float32)
        y = y.astype(np.float32)

        # metacell-level association
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

        # build metacell weights per cell type
        sig_w = norm.cdf(z)

        metacell_weight = {}
        for ct, f in freq_df.iterrows():
            w = sig_w * np.asarray(f, dtype=np.float32)
            tot = w.sum()
            metacell_weight[ct] = w / (tot if tot > 0 else 1e-12)

        # null t distribution for each cell type
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

        # aggregate to cell-type level
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

        # estimate mixture of association
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

        ctdfbs = None
        if self.output_dfbs:
            # calculate influence per cell type
            t0 = time()
            logger.info(f"[influence] Get influence score for given {self.ct_key}")
            ctdfbs = np.vstack([i for i in metacell_weight.values()]) @ dfb / ct_df['se'].values[:, None]
            logger.info(f"[influence] finished in {(time() - t0) / 60:.2f} min")

        return ct_df.sort_values("q", ignore_index=True), mc_df, ctdfbs
