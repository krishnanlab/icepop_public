#!/usr/bin/env python

import os
import time
import logging
import subprocess
import papermill
from enrichment_analysis import EnrichmentPipeline

logger = logging.getLogger("icepop.interactive")

def interactive(outdir, geneset_collections, notebook, geneset_path=None, adata_path=None):
    """
    Run enrichment analysis and execute notebook
    outdir : str
    geneset_collections : str
    notebook : str
    geneset_path : str, optional
    """


    if not outdir:
        raise ValueError("outdir is required")

    if not notebook:
        raise ValueError("notebook path is required")

    if not os.path.exists(notebook):
        raise FileNotFoundError(f"Notebook not found: {notebook}")

    if geneset_collections.lower() == "none" and not geneset_path:
        raise ValueError("--geneset_path required when geneset_collections=None")

    if geneset_path and not os.path.exists(geneset_path):
        raise FileNotFoundError(f"geneset_path not found: {geneset_path}")


    os.makedirs(outdir, exist_ok=True)
    logger.info(f"Enrichment output directory: {outdir}")
    logger.info("Starting enrichment analysis...")
    start = time.perf_counter()
    try:
        pipeline = EnrichmentPipeline(
            outdir=outdir,
            geneset_collections=geneset_collections,
            geneset_path=geneset_path,)
        pipeline.run()
    except Exception as e:
        logger.error(f"Enrichment analysis failed: {e}")
        raise

    elapsed = time.perf_counter() - start
    logger.info(f"Enrichment completed in {elapsed:.2f} seconds")
    logger.info("Generating summary report...")
    try:
        executed_nb = os.path.join(outdir, "icepop-report.ipynb")
        subprocess.run([
            "papermill",
            notebook,
            executed_nb,
            "-p","outdir", outdir,
            "-p","adata_path",adata_path], check=True)
        subprocess.run([
            "jupyter", "nbconvert",
            "--to", "html",
            executed_nb,
            "--output-dir", outdir].check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running notebook: {e}")
        raise

    logger.info(f"Notebook executed and saved to {outdir}")