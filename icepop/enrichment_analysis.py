import os
import pickle
import logging
import pandas as pd
import numpy as np
from math import sqrt
from scipy.stats import hypergeom
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests


class EnrichmentPipeline:
    def __init__(self, outdir: str, geneset_collections: str, geneset_path: str = None):
        self.outdir = outdir
        self.geneset_collections = geneset_collections
        self.geneset_path = geneset_path

    @staticmethod
    def gene_contribution(gene_contribution_path: str, factor: float):
        """Load and filter gene contribution data."""

        gene_contribution = np.load(gene_contribution_path, allow_pickle=True)

        gene_contribution_table = pd.DataFrame(
            gene_contribution["dfbs"],
            index=gene_contribution["celltypes"],
            columns=gene_contribution["genes"],
        )

        gc_threshold = factor / sqrt(gene_contribution_table.shape[1])
        is_gene = (gene_contribution_table > gc_threshold).any(axis=0)
        gene_contribution_filtered = gene_contribution_table.loc[:, is_gene]

        return gene_contribution_table, gc_threshold, gene_contribution_filtered

    def parse_geneset_args(self):
        """Parse arguments for geneset collections"""

        collections_arg = self.geneset_collections

        if collections_arg.lower() == "none":
            if not self.geneset_path:
                raise ValueError(
                    "When --geneset_collections=None, you must provide --geneset_path."
                )
            return None
        if collections_arg.lower() == "all":
            return ["All"]
        return [c.strip() for c in collections_arg.split(",")]

    @staticmethod
    def universe_genes(
        gene_contribution_table: pd.DataFrame,
        geneset_diction_filtered: dict,
        prefix: str,
    ) -> set[str]:
        """Get universe of genes for geneset collections."""

        geneset_all_genes = [
            gene
            for genes in geneset_diction_filtered[prefix].values()
            for gene in genes
        ]
        universe = set(gene_contribution_table.columns).intersection(
            set(geneset_all_genes)
        )
        return universe

    @staticmethod
    def run_hypergeometric(
        cell_gene: pd.DataFrame,
        celltype: str,
        universe: set[str],
        geneset_dict: dict,
        gc_threshold: float,
    ) -> pd.DataFrame:
        """Run hypergeometric test for a given cell type against a give geneset collection"""

        geneset_results = []
        universe_len = len(universe)

        is_cell = cell_gene.loc[celltype] > gc_threshold
        cell_genelist = set(cell_gene.columns[is_cell]).intersection(universe)
        cell_genelist_len = len(cell_genelist)

        for gset_name, gene_set in geneset_dict.items():
            gene_set = set(gene_set).intersection(universe)
            geneset_len = len(gene_set)
            overlap_genes = cell_genelist.intersection(gene_set)
            overlap_len = len(overlap_genes)

            union_len = len(cell_genelist.union(gene_set))
            normalized_overlap = overlap_len / universe_len

            hypergeom_test = hypergeom(universe_len, geneset_len, cell_genelist_len)
            pval = hypergeom_test.sf(overlap_len - 1)
            geneset_results.append(
                {
                    "celltype": celltype,
                    "geneset": gset_name,
                    "geneset_size": geneset_len,
                    "celltype_genelist_size": cell_genelist_len,
                    "overlap": overlap_len,
                    "universe": universe_len,
                    "normalized_overlap": normalized_overlap,
                    "pvalue": pval,
                }
            )
        return pd.DataFrame(geneset_results)

    @staticmethod
    def extract_trait(file_path: str):
        fname = os.path.basename(file_path)
        return fname.split("__")[1].split(".")[0]

    @staticmethod
    def load_msigdb(path: str):
        with open(path, "rb") as f:
            return pickle.load(f)

    @staticmethod
    def find_gene_contribution_file(outdir: str):
        return next(
            (
                os.path.join(outdir, f)
                for f in os.listdir(outdir)
                if f.startswith("dfbs") and f.endswith(".npz")
            ),
            None,
        )

    @staticmethod
    def significant_celltypes(outdir: str):
        celltype_assoc_file = next(
            (
                os.path.join(outdir, f)
                for f in os.listdir(outdir)
                if f.startswith("celltype__trait") and f.endswith(".csv")
            ),
            None,
        )
        if celltype_assoc_file is None:
            raise FileNotFoundError("No celltype association file found.")

        celltype_assoc = pd.read_csv(celltype_assoc_file)
        celltype_assoc_filtered = celltype_assoc[
            (celltype_assoc["q"] < 0.1) & (celltype_assoc["sig_pct"] > 0.2)
        ].reset_index(drop=True)
        return celltype_assoc_filtered["cell_type"].tolist()

    @staticmethod
    def load_gmt(gmt_path: str) -> dict:
        """loading user provided gmt file"""
        genesets = {}
        with open(gmt_path, "r") as f:
            for line in f:
                geneset_name = line.strip().split("\t")[0]
                genes = line.strip().split("\t")[1:]
                genesets[geneset_name] = genes
        return genesets

    def write_enrichment(
        self,
        trait: str,
        celltypes: list,
        gene_contribution_table: pd.DataFrame,
        gene_contribution_filtered: pd.DataFrame,
        geneset_diction_filtered: dict,
        gc_threshold: float,
    ):

        for prefix, genesets in geneset_diction_filtered.items():
            if not genesets:
                continue

            universe = self.universe_genes(
                gene_contribution_table, geneset_diction_filtered, prefix
            )
            results_list = Parallel(n_jobs=-1, backend="loky")(
                delayed(self.run_hypergeometric)(
                    cell_gene=gene_contribution_filtered,
                    celltype=ct,
                    universe=universe,
                    geneset_dict=genesets,
                    gc_threshold=gc_threshold,
                )
                for ct in celltypes
            )

            for hypergeom_results in results_list:
                hypergeom_results["prefix"] = prefix
                hypergeom_results["fdr"] = multipletests(
                    hypergeom_results["pvalue"], method="fdr_bh"
                )[1]

            genesets_all = (
                pd.concat(results_list, ignore_index=True)
                .sort_values("fdr", ascending=True)
                .reset_index(drop=True)
            )
            enrich_outdir = os.path.join(self.outdir, "enrichment")
            os.makedirs(enrich_outdir, exist_ok=True)
            output_path = os.path.join(
                enrich_outdir, f"enrichment_anlaysis_{trait}_{prefix}.txt"
            )
            genesets_all.to_csv(output_path, sep="\t", index=False)

    def run(self):
        gene_file = self.find_gene_contribution_file(self.outdir)
        celltypes = self.significant_celltypes(self.outdir)
        if gene_file is None:
            raise FileNotFoundError("No gene contribution file found.")
        trait = self.extract_trait(gene_file)
        gene_df, gc_threshold, gene_df_filtered = self.gene_contribution(
            gene_file, factor=2
        )
        collections = self.parse_geneset_args()

        if collections is None:
            custom_genesets = self.load_gmt(self.geneset_path)
            self.write_enrichment(
                trait,
                celltypes,
                gene_df,
                gene_df_filtered,
                {"CUSTOM": custom_genesets},
                gc_threshold,
            )
            return

        geneset_pkl_path = "/mnt/research/FishEvoDevoGeno/Hao/icepop_public/genesets/msigdb_genesets.pkl"
        msigdb = self.load_msigdb(geneset_pkl_path)
        if collections == ["All"]:
            collections = list(msigdb.keys())
        for gc in collections:
            self.write_enrichment(
                trait,
                celltypes,
                gene_df,
                gene_df_filtered,
                {gc: msigdb[gc]},
                gc_threshold,
            )
