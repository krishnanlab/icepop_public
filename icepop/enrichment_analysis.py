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
    def __init__(self, 
                results_path: str,
                geneset_collections: str,
                geneset_path: str=None):
        self.results_path = results_path
        self.geneset_collections = geneset_collections
        self.geneset_path = geneset_path

    @staticmethod
    def gene_contribution(gene_contribution_path: str,
                          factor: float):
        """Load and filter gene contribution data."""
        
        gene_contribution = np.load(
            gene_contribution_path,
            allow_pickle=True)
        
        gene_contribution_table = pd.DataFrame(
            gene_contribution['dfbs'],
            index=gene_contribution['celltypes'],
            columns=gene_contribution['genes'])
        
        gc_threshold = factor/sqrt(gene_contribution_table.shape[1])
        is_gene = (gene_contribution_table > gc_threshold).any(axis=0)
        gene_contribution_filtered = gene_contribution_table.loc[:,is_gene]

        return gene_contribution_table, gc_threshold, gene_contribution_filtered


    def parse_geneset_args(self):
        """Parse arguments for geneset collections"""
        
        collections_arg = self.geneset_collections
        
        if collections_arg.lower() == "none":
            if not args['geneset_path']:
                raise ValueError(
                    "When --geneset_collections=None, you must provide --geneset_path."
                    )
            return None
        if collections_arg.lower() == "all":
            return ["All"]    
        return [c.strip() for c in collections_arg.split(",")]


    @staticmethod
    def universe_genes(gene_contribution_table: pd.DataFrame,
                       geneset_diction_filtered: dict,
                       prefix: str) -> set[str]:
        """Get universe of genes for geneset collections."""

        geneset_all_genes = [gene 
                        for genes in geneset_diction_filtered[prefix].values()
                        for gene in genes]
        universe = set(gene_contribution_table.columns).intersection(set(geneset_all_genes))
        return universe


    @staticmethod
    def run_hypergeometric(cell_gene: pd.DataFrame,
                           celltype: str,
                           universe: set[str],
                           geneset_dict: dict,
                           gc_threshold: float) -> pd.DataFrame:
        """Run hypergeometric test for a given cell type against a give geneset collection"""

        geneset_results = []
        universe_len = len(universe)
        
        is_cell = cell_gene.loc[celltype] > gc_threshold
        cell_genelist =  set(cell_gene.columns[is_cell]).intersection(universe)
        cell_genelist_len = len(cell_genelist)

        for gset_name, gene_set in geneset_dict.items():
            gene_set = set(gene_set).intersection(universe)
            geneset_len = len(gene_set)
            overlap_genes = cell_genelist.intersection(gene_set)
            overlap_len = len(overlap_genes)
            hypergeom_test = hypergeom(universe_len, 
                                        geneset_len, cell_genelist_len)
            pval = hypergeom_test.sf(overlap_len-1)
            geneset_results.append({
                    "celltype": celltype,
                    "geneset": gset_name,
                    "geneset_size": geneset_len,
                    "celltype_genelist_size": cell_genelist_len,
                    "overlap": overlap_len,
                    "universe": universe_len,
                    "pvalue": pval,
                })
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
    def find_gene_contribution_file(results_path: str):
        return next(
            (
                os.path.join(results_path, f)
                for f in os.listdir(results_path)
                if f.startswith("dfbs") and f.endswith(".npz")
            ),
            None,
        )
    

    def write_enrichment(self,
                        trait: str,
                        gene_contribution_table: pd.DataFrame,
                        gene_contribution_filtered: pd.DataFrame,
                        geneset_diction_filtered: pd.DataFrame,
                        gc_threshold: float):
        for prefix, genesets in geneset_diction_filtered.items():
            if not genesets:
                continue
            
            universe = self.universe_genes(
                gene_contribution_table, 
                geneset_diction_filtered, 
                prefix)
            celltypes = gene_contribution_filtered.index.tolist()
            results_list = Parallel(n_jobs=-1, backend="loky")(
                delayed(self.run_hypergeometric)(
                    cell_gene=gene_contribution_filtered,
                    celltype=ct,
                    universe=universe,
                    geneset_dict=genesets,
                    gc_threshold=gc_threshold,)
                for ct in celltypes)

            for hypergeom_results in results_list:
                hypergeom_results["prefix"] = prefix
                hypergeom_results["fdr"] = multipletests(hypergeom_results["pvalue"], method = "fdr_bh")[1]

                genesets_all = (pd.concat(results_list, ignore_index=True)
                        .sort_values("fdr", ascending=True)
                        .reset_index(drop=True))
                output_path = os.path.join(self.results_path, 
                                    f"enrichment_anlaysis_{trait}_{prefix}.txt")
                genesets_all.to_csv(output_path, sep="\t", index=False)


    def run(self):
        gene_file = self.find_gene_contribution_file(self.results_path)
        if gene_file is None:
            raise FileNotFoundError("No gene contribution file found.")

        trait = self.extract_trait(gene_file)

        gene_df, gc_threshold, gene_df_filtered = self.gene_contribution(
            gene_file, factor=2)

        geneset_pkl_path = "./data/msigdb_genesets.pkl"
        msigdb = self.load_msigdb(geneset_pkl_path)

        collections = self.parse_geneset_args()

        if collections is None:
            raise NotImplementedError("Custom geneset not implemented yet")

        if collections == ["All"]:
            collections = list(msigdb.keys())

        for gc in collections:
            self.write_enrichment(
                trait,
                gene_df,
                gene_df_filtered,
                {gc: msigdb[gc]},
                gc_threshold,)                    
