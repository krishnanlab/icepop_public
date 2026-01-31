import numpy as np
import pandas as pd
from icepop.data import HomologyData


class CrossSpeciesScoreConverter:
    """
    class to calculate expression specificity score
    """

    def __init__(
        self,
        adata,
        sp='mmusculus',
    ):
        self.adata = adata
        self.sp = sp

    def generate_cross_sp_matrix(self):
        # get human to model sp genes map
        ortho_map = HomologyData(sp=self.sp).load()

        # Get all species genes
        all_sp_genes = set(self.adata.var_names)

        # Filter the ortholog map
        filtered_ortho_map = {
            human_gene: list(set(sp_genes) & all_sp_genes)
            for human_gene, sp_genes in ortho_map.items() if set(sp_genes) & all_sp_genes
        }

        # generate a 0/1 conversion matrix orthologs
        sp_gene_idx_dict = {gene: idx for idx, gene in enumerate(self.adata.var_names)}

        # Precompute indices for matrix assignment
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
        if self.sp != 'hsapiens':
            self.generate_cross_sp_matrix()

            # Vectorized matrix multiplication for score conversion
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
