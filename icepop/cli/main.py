#!/usr/bin/env python
"""
ICePop command-line interface.

This module defines the top-level CLI entry point for ICePop,
exposing subcommands for metacell construction, association testing,
and downstream analyses.

Example
-------
icepop metacell --h5ad data.h5ad --outfile results/metaq
icepop association --h5ad data.h5ad --magmaz trait.genes.out
"""

import fire
import logging
from icepop.metacell import metacell
from icepop.association import association
# from icepop.interactive import interactive


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def main(verbose: bool = False):
    setup_logging(verbose)
    fire.Fire({
        "metacell": metacell,
        "association": association,
        # "interactive": interactive,
    })


if __name__ == "__main__":
    main()
