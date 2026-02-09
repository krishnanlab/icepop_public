import time
import logging
import fire
from enrichment_analysis import EnrichmentPipeline

logging.basicConfig(level=logging.INFO)

def run_enrichment(results_path: str,
                    geneset_collections: str,
                    geneset_path: str | None = None,):
    """
    Run enrichment analysis using hypergeometric test.

    Args:
        results_path: Path to results directory.
        geneset_collections: Geneset collections: All, KEGG, REACTOME, etc.
        geneset_path: Optional path to custom geneset GMT file.
    """
    start = time.perf_counter()

    pipeline = EnrichmentPipeline(
        results_path=results_path,
        geneset_collections=geneset_collections,
        geneset_path=geneset_path,
    )

    logging.info("Starting enrichment analysis")
    pipeline.run()

    elapsed = time.perf_counter() - start
    logging.info(f"Execution time: {elapsed:.2f} seconds")

    return elapsed  

def main():
    fire.Fire(run_enrichment)

if __name__ == "__main__":
    main()
