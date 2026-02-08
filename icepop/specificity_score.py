import numpy as np
import pandas as pd
from scipy import stats
from time import time
from scipy.sparse import csr_matrix
from collections import defaultdict
from multiprocessing import Pool
from multiprocessing.shared_memory import SharedMemory
import logging

logger = logging.getLogger(__name__)


class specificity_score:
    """
    Compute gene expression specificity scores across metacells or cell types.

    This class calculates, for each gene, the probability that expression
    within a given group (metacell or cell type) is higher than expression
    in the remaining cells. The computation is optimized for sparse matrices
    and parallel execution using shared memory.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset containing:

        - ``adata.X`` → gene expression matrix (cells × genes)
        - ``adata.obs['metacell']`` → metacell labels (for metacell scoring)
        - ``adata.obs['cell_type']`` → cell type labels (for cell-type scoring)
    n_jobs : int, default=1
        Number of parallel worker processes used for score computation.

    Notes
    -----
    - Uses sparse CSR representation for memory efficiency.
    - Parallel workers attach to shared memory instead of copying matrices.
    - Final specificity score combines:

        **Normal CDF(z-statistic) × expression coverage**

      followed by column-wise normalization.
    """

    def __init__(
        self,
        adata,
        n_jobs=1,
    ):
        """Initialize specificity score calculator."""
        self.adata = adata
        self.n_jobs = n_jobs

    @staticmethod
    def calculate_statistics(data):
        """
        Compute per-gene mean and variance from a sparse matrix.

        Parameters
        ----------
        data : scipy.sparse matrix (genes × cells)
            Sparse expression matrix subset.

        Returns
        -------
        mean : ndarray of shape (n_genes,)
            Mean expression per gene.
        var : ndarray of shape (n_genes,)
            Variance per gene.
        count : int
            Number of cells used in the calculation.

        Notes
        -----
        Variance is set to zero when only one cell is present.
        """
        count = data.shape[1]
        sum_mat = data.sum(axis=1)
        square_sum_mat = data.power(2).sum(axis=1)
        mean_mat = sum_mat / count
        if count > 1:
            var_mat = (square_sum_mat - (np.power(sum_mat, 2) / count)) / (count - 1)
        else:
            # for metacell with only one cell, var is 0
            var_mat = np.zeros_like(square_sum_mat, dtype=square_sum_mat.dtype)
        return np.asarray(mean_mat).flatten(), np.asarray(var_mat).flatten(), count

    @staticmethod
    def compute_score_for_metacell(args):
        """
        Worker function to compute specificity statistics for one metacell.

        This function reconstructs the shared sparse matrix, extracts the
        subset of cells belonging to a metacell, and computes:

        - mean/variance within the metacell
        - inferred mean/variance of remaining cells
        - z-statistic comparing the two
        - expression coverage within the metacell

        Returns
        -------
        stat : ndarray
            Z-statistic per gene.
        rc : ndarray
            Expression coverage ratio per gene.
        """
        # unpack task parameters
        (global_mean, global_var, global_count, global_mean_sq, idx_list,
         calculate_statistics, shm_data_name,
         shm_indices_name, shm_indptr_name, data_dtype, indices_dtype,
         indptr_dtype, data_shape, indices_shape, indptr_shape, csr_shape) = args

        # attach to shared memory
        shm_data = SharedMemory(name=shm_data_name)
        shm_indices = SharedMemory(name=shm_indices_name)
        shm_indptr = SharedMemory(name=shm_indptr_name)

        try:
            # reconstruct arrays from shared memory
            data = np.ndarray(data_shape, dtype=data_dtype, buffer=shm_data.buf)
            indices = np.ndarray(indices_shape, dtype=indices_dtype, buffer=shm_indices.buf)
            indptr = np.ndarray(indptr_shape, dtype=indptr_dtype, buffer=shm_indptr.buf)

            # reconstruct CSR matrix
            data_mat = csr_matrix((data, indices, indptr), shape=csr_shape)

            # slice the matrix using idx_list (columns)
            sub_data = data_mat[:, idx_list]

            # calculate stats for metacells
            sub_mean, sub_var, sub_count = calculate_statistics(sub_data)

            # infer mean and var of rest cells based on global mean and var
            rest_count = global_count - sub_count
            rest_mean = (global_mean * global_count - sub_mean * sub_count) / rest_count
            rest_var = (global_count * global_var - sub_count * (sub_var + np.power(sub_mean, 2) - 2 * sub_mean * global_mean + global_mean_sq)) / rest_count - ((np.power(rest_mean, 2) - 2 * rest_mean * global_mean + global_mean_sq))

            # calculate prob of metacell expressed higher than the rest
            denom = np.sqrt(sub_var / sub_count + rest_var / rest_count)
            denom = np.where(denom == 0, 1e-10, denom)  # avoid divide-by-zero
            stat = (sub_mean - rest_mean) / denom

            # expression coverage across current cell types
            rc = np.asarray((sub_data > 0).sum(axis=1)).flatten() / sub_count
            return stat, rc

        finally:
            # close shared memory handles in worker
            shm_data.close()
            shm_indices.close()
            shm_indptr.close()

    def get_metacell_spec_score(self):
        """
        Compute specificity scores for all genes across metacells.

        Results are stored in:

        ``adata.uns["spec_score"]`` → DataFrame (metacells × genes)

        Notes
        -----
        - Uses multiprocessing with shared memory.
        - Scores are:

          ``NormalCDF(z) × coverage``

          followed by column-wise normalization.
        """
        logger.info('Calculate metacell specificity score')

        t0 = time()
        # transpose of input matrix and convert to sparse matrix
        data = self.adata.X.transpose().tocsr()
        # extract CSR components
        data_array = data.data
        indices_array = data.indices
        indptr_array = data.indptr
        csr_shape = data.shape

        # create shared memory blocks
        shm_data = SharedMemory(create=True, size=data_array.nbytes)
        shm_indices = SharedMemory(create=True, size=indices_array.nbytes)
        shm_indptr = SharedMemory(create=True, size=indptr_array.nbytes)

        # copy data to shared memory
        np_data_shm = np.ndarray(data_array.shape, dtype=data_array.dtype, buffer=shm_data.buf)
        np_data_shm[:] = data_array[:]

        np_indices_shm = np.ndarray(indices_array.shape, dtype=indices_array.dtype, buffer=shm_indices.buf)
        np_indices_shm[:] = indices_array[:]

        np_indptr_shm = np.ndarray(indptr_array.shape, dtype=indptr_array.dtype, buffer=shm_indptr.buf)
        np_indptr_shm[:] = indptr_array[:]

        # calculate var, mean and mean square for entire matrix
        global_mean, global_var, global_count = self.calculate_statistics(data)
        global_mean_sq = np.power(global_mean, 2)

        # get metacell and corresponding cell idx
        metacell2idx = defaultdict(list)
        for metacell, idx in zip(self.adata.obs["metacell"], np.arange(self.adata.shape[0])):
            metacell2idx[metacell].append(idx)

        # prepare task parameters
        tasks = [
            (
                global_mean, global_var, global_count, global_mean_sq,
                idx_list,
                self.calculate_statistics,
                shm_data.name,
                shm_indices.name,
                shm_indptr.name,
                data_array.dtype,
                indices_array.dtype,
                indptr_array.dtype,
                data_array.shape,
                indices_array.shape,
                indptr_array.shape,
                csr_shape
            )
            for idx_list in metacell2idx.values()
        ]
        logger.info('Took %.2f min to parepare input' % ((time() - t0) / 60))

        t0 = time()
        try:
            with Pool(self.n_jobs) as pool:
                res = pool.map(self.compute_score_for_metacell, tasks)
        finally:
            # clean up shared memory in main process
            shm_data.close()
            shm_data.unlink()
            shm_indices.close()
            shm_indices.unlink()
            shm_indptr.close()
            shm_indptr.unlink()
        logger.info('Took %.2f min to calculate metacell specificity score' % ((time() - t0) / 60))

        # aggregate into a matrix, then normalized score across cells
        stat, rc = zip(*res)
        stat = np.vstack(stat).astype(np.float32)
        rc = np.vstack(rc).astype(np.float32)

        s_scores = stats.norm.cdf(stat) * rc
        self.adata.uns["spec_score"] = pd.DataFrame(
            s_scores / (np.sum(s_scores, axis=0) + 1e-12), index=list(metacell2idx.keys()), columns=list(self.adata.var_names)
        )

    @staticmethod
    def compute_score_for_celltype(args):
        """
        Worker function computing specificity statistics for one cell type.

        Identical statistical formulation to metacell scoring, but applied
        to groups defined by ``adata.obs['cell_type']``.
        """
        # unpack task parameters
        (global_mean, global_var, global_count, global_mean_sq, idx_list,
         calculate_statistics, shm_data_name,
         shm_indices_name, shm_indptr_name, data_dtype, indices_dtype,
         indptr_dtype, data_shape, indices_shape, indptr_shape, csr_shape) = args

        # attach to shared memory
        shm_data = SharedMemory(name=shm_data_name)
        shm_indices = SharedMemory(name=shm_indices_name)
        shm_indptr = SharedMemory(name=shm_indptr_name)

        try:
            # reconstruct arrays from shared memory
            data = np.ndarray(data_shape, dtype=data_dtype, buffer=shm_data.buf)
            indices = np.ndarray(indices_shape, dtype=indices_dtype, buffer=shm_indices.buf)
            indptr = np.ndarray(indptr_shape, dtype=indptr_dtype, buffer=shm_indptr.buf)

            # reconstruct CSR matrix
            data_mat = csr_matrix((data, indices, indptr), shape=csr_shape)

            # slice the matrix using idx_list (columns)
            sub_data = data_mat[:, idx_list]

            # calculate stats for cell types
            sub_mean, sub_var, sub_count = calculate_statistics(sub_data)

            # infer mean and var of rest cells based on global mean and var
            rest_count = global_count - sub_count
            rest_mean = (global_mean * global_count - sub_mean * sub_count) / rest_count
            rest_var = (global_count * global_var - sub_count * (sub_var + np.power(sub_mean, 2) - 2 * sub_mean * global_mean + global_mean_sq)) / rest_count - ((np.power(rest_mean, 2) - 2 * rest_mean * global_mean + global_mean_sq))

            # calculate prob of cell types expressed higher than the rest
            denom = np.sqrt(sub_var / sub_count + rest_var / rest_count)
            denom = np.where(denom == 0, 1e-10, denom)  # avoid divide-by-zero
            stat = (sub_mean - rest_mean) / denom

            # expression coverage across current cell types
            rc = np.asarray((sub_data > 0).sum(axis=1)).flatten() / sub_count
            return stat, rc

        finally:
            # close shared memory handles in worker
            shm_data.close()
            shm_indices.close()
            shm_indptr.close()

    def get_celltype_spec_score(self):
        """
        Compute specificity scores for all genes across cell types.

        Results are stored in:

        ``adata.uns["cell_type_spec_score"]`` → DataFrame (cell types × genes)
        """
        logger.info('Calculate cell type specificity score')

        t0 = time()
        # transpose of input matrix and convert to sparse matrix
        data = self.adata.X.transpose().tocsr()
        # extract CSR components
        data_array = data.data
        indices_array = data.indices
        indptr_array = data.indptr
        csr_shape = data.shape

        # create shared memory blocks
        shm_data = SharedMemory(create=True, size=data_array.nbytes)
        shm_indices = SharedMemory(create=True, size=indices_array.nbytes)
        shm_indptr = SharedMemory(create=True, size=indptr_array.nbytes)

        # copy data to shared memory
        np_data_shm = np.ndarray(data_array.shape, dtype=data_array.dtype, buffer=shm_data.buf)
        np_data_shm[:] = data_array[:]

        np_indices_shm = np.ndarray(indices_array.shape, dtype=indices_array.dtype, buffer=shm_indices.buf)
        np_indices_shm[:] = indices_array[:]

        np_indptr_shm = np.ndarray(indptr_array.shape, dtype=indptr_array.dtype, buffer=shm_indptr.buf)
        np_indptr_shm[:] = indptr_array[:]

        # calculate var, mean and mean square for entire matrix
        global_mean, global_var, global_count = self.calculate_statistics(data)
        global_mean_sq = np.power(global_mean, 2)

        # get cell type and corresponding cell idx
        celltype2idx = defaultdict(list)
        for cell_type, idx in zip(self.adata.obs["cell_type"], np.arange(self.adata.shape[0])):
            celltype2idx[cell_type].append(idx)

        # prepare task parameters
        tasks = [
            (
                global_mean, global_var, global_count, global_mean_sq,
                idx_list,
                self.calculate_statistics,
                shm_data.name,
                shm_indices.name,
                shm_indptr.name,
                data_array.dtype,
                indices_array.dtype,
                indptr_array.dtype,
                data_array.shape,
                indices_array.shape,
                indptr_array.shape,
                csr_shape
            )
            for idx_list in celltype2idx.values()
        ]
        logger.info('Took %.2f min to parepare input' % ((time() - t0) / 60))

        t0 = time()
        try:
            with Pool(self.n_jobs) as pool:
                res = pool.map(self.compute_score_for_celltype, tasks)
        finally:
            # clean up shared memory in main process
            shm_data.close()
            shm_data.unlink()
            shm_indices.close()
            shm_indices.unlink()
            shm_indptr.close()
            shm_indptr.unlink()
        logger.info('Took %.2f min to calculate cell type specificity score' % ((time() - t0) / 60))

        # aggregate into a matrix, then normalized score across cells
        stat, rc = zip(*res)
        s_scores = stats.norm.cdf(np.vstack(stat)) * np.vstack(rc)
        self.adata.uns["cell_type_spec_score"] = pd.DataFrame(
            s_scores / np.sum(s_scores, axis=0) + 1e-12, index=list(celltype2idx.keys()), columns=list(self.adata.var_names)
        )
