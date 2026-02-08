import numpy as np
import pandas as pd
from icepop.data import HomologyData


class CrossSpeciesScoreConverter:
    """
    Convert gene-level scores between species using ortholog mappings.

    This class builds a binary ortholog conversion matrix that maps
    model-organism genes (e.g., mouse) to human genes, enabling
    cross-species projection of gene-level scores such as expression
    specificity, enrichment, or association statistics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing gene expression or specificity
        scores. Gene identifiers must be stored in ``adata.var_names``.
    sp : str, default='mmusculus'
        Species identifier for the source organism. If ``'hsapiens'``,
        no conversion is performed.

    Notes
    -----
    - Ortholog relationships are loaded using :class:`icepop.data.HomologyData`.
    - Conversion is performed via matrix multiplication between the
      input score matrix and a binary ortholog mapping matrix.
    - Optional normalization accounts for the number of orthologs per
      human gene.
    """

    def __init__(
        self,
        adata,
        sp='mmusculus',
    ):
        """
        Initialize the cross-species score converter.

        Parameters
        ----------
        adata : AnnData
            Annotated data object containing gene-level information.
        sp : str, default='mmusculus'
            Source species name used for ortholog lookup.
        """
        self.adata = adata
        self.sp = sp

    def generate_cross_sp_matrix(self):
        """
        Construct the ortholog conversion matrix from source species to human.

        This method:

        1. Loads human-to-species ortholog mappings.
        2. Filters mappings to genes present in ``adata.var_names``.
        3. Builds a binary matrix of shape:

           ``(n_species_genes, n_human_genes)``

           where
           - Rows correspond to source-species genes.
           - Columns correspond to human genes.
           - Entries are 1 if an ortholog relationship exists.

        Attributes Created
        ------------------
        human_genes_sorted : list of str
            Sorted list of human genes included in the mapping.
        ortho_mat : np.ndarray
            Binary ortholog mapping matrix.
        n_gene_mat : np.ndarray
            Number of orthologs per human gene (column-wise sum of
            ``ortho_mat``).
        """
        # get human to model sp genes map
        ortho_map = HomologyData(sp=self.sp).load()

        # get all species genes
        all_sp_genes = set(self.adata.var_names)

        # filter the ortholog map
        filtered_ortho_map = {
            human_gene: list(set(sp_genes) & all_sp_genes)
            for human_gene, sp_genes in ortho_map.items() if set(sp_genes) & all_sp_genes
        }

        # generate a 0/1 conversion matrix orthologs
        sp_gene_idx_dict = {gene: idx for idx, gene in enumerate(self.adata.var_names)}

        # precompute indices for matrix assignment
        self.human_genes_sorted = sorted(filtered_ortho_map.keys())
        num_sp_genes = len(all_sp_genes)
        num_human_genes = len(self.human_genes_sorted)
        self.ortho_mat = np.zeros((num_sp_genes, num_human_genes), dtype=np.float32)

        col_indices = []
        row_indices = []
        for col_idx, human_gene in enumerate(self.human_genes_sorted):
            for sp_gene in filtered_ortho_map[human_gene]:
                row_indices.append(sp_gene_idx_dict[sp_gene])
                col_indices.append(col_idx)
        self.ortho_mat[row_indices, col_indices] = 1

        self.n_gene_mat = np.sum(self.ortho_mat, axis=0)

    def convert_score_across_species(self, score, normed=True):
        """
        Project gene-level scores from the source species to human genes.

        Parameters
        ----------
        score : pandas.DataFrame
            Gene-level score matrix with rows representing samples or
            cell types and columns corresponding to source-species genes
            aligned with ``adata.var_names``.
        normed : bool, default=True
            Whether to normalize projected scores by the square root of
            the number of orthologs per human gene.

            - ``True`` → divide by ``sqrt(n_orthologs)`` (recommended for
              variance stabilization).
            - ``False`` → divide by ``n_orthologs`` (simple averaging).

        Returns
        -------
        pandas.DataFrame
            Score matrix in **human gene space** with:

            - Same row index as ``score``.
            - Columns equal to mapped human genes.

            If ``sp == 'hsapiens'``, the input score is returned unchanged.

        Notes
        -----
        Conversion is computed efficiently using matrix multiplication:

        ``score @ ortholog_matrix``

        followed by normalization based on ortholog counts.
        """
        if self.sp != 'hsapiens':
            self.generate_cross_sp_matrix()

            # vectorized matrix multiplication for score conversion
            # average double norm score
            if normed:
                return pd.DataFrame(
                    score.to_numpy() @ self.ortho_mat / np.sqrt(self.n_gene_mat),
                    index=score.index,
                    columns=self.human_genes_sorted,
                )
            else:
                return pd.DataFrame(
                    score.to_numpy() @ self.ortho_mat / self.n_gene_mat,
                    index=score.index,
                    columns=self.human_genes_sorted,
                )
        else:
            return score
