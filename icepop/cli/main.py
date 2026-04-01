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
from icepop.interactive import interactive


def setup_logging(verbose=False):
    """
    Configure global logging behavior for the ICePop CLI.

    Parameters
    ----------
    verbose : bool, default=False
        If ``True``, sets logging level to ``DEBUG`` for detailed output.
        Otherwise, uses ``INFO`` level for standard progress messages.

    Notes
    -----
    This function initializes the root logger using ``logging.basicConfig``.
    It should be called before invoking any CLI subcommands.
    """
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def main(verbose: bool = False):
    """
    Entry point for the ICePop command-line interface.

    This function initializes logging and dispatches CLI subcommands
    using the ``fire`` library.

    Parameters
    ----------
    verbose : bool, default=False
        Enable debug-level logging output.

    Available Subcommands
    ---------------------
    metacell
        Construct metacells from single-cell expression data.
    association
        Perform gene-level and cell-type–level association testing.
    interactive
        Generate interactive output to explore affected cells and influential genes
    """
    setup_logging(verbose)
    fire.Fire({
        "metacell": metacell,
        "association": association,
        "interactive": interactive,
    })


if __name__ == "__main__":
    main()
