import os
import time
import argparse
import logging
from enrichment_analysis import EnrichmentPipeline

def args_parser():
    parser = argparse.ArgumentParser(
        description="Run enrichment analysis using hypergeometric test.")
    
    parser.add_argument("--results_path",
                        type=str,
                        required=True,
                        help="Path to results directory.",)

    parser.add_argument("--geneset_collections",
                        type=str,
                        required=True,
                        help="Geneset collections: All, KEGG, REACTOME, etc.",)

    parser.add_argument("--geneset_path",
                        type=str,
                        required=False,
                        help="Path to custom geneset GMT file.",)
    
    return parser.parse_args()


def main():
    start = time.perf_counter()
    args = args_parser()

    pipeline = EnrichmentPipeline(results_path=args.results_path,
                                geneset_collections=args.geneset_collections,
                                geneset_path=args.geneset_path,)

    logging.info("Starting enrichment analysis")
    pipeline.run()

    elapsed = time.perf_counter() - start
    logging.info(f"Execution time: {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()