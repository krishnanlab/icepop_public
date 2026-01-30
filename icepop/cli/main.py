#!/usr/bin/env python
import fire
from icepop.metacell import metacell
from icepop.association import association
# from icepop.influence import influence
# from icepop.interactive import interactive


def main():
    fire.Fire({
        "metacell": metacell,
        "association": association,
        # "influence": influence,
        # "interactive": interactive,
    })


if __name__ == "__main__":
    main()
