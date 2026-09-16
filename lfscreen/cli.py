"""Command-line interface: ``lfscreen candidates.csv -o ranked.csv --mode indoor``."""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from .screen import BANDGAP_WINDOWS, rank_candidates


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="lfscreen",
        description="Rank lead-free perovskite / perovskite-inspired PV absorbers.",
    )
    p.add_argument("csv", help="input CSV with 'formula' and optional 'band_gap' (eV) columns")
    p.add_argument("-o", "--output", help="output CSV (default: print to stdout)")
    p.add_argument("--mode", choices=sorted(BANDGAP_WINDOWS), default="outdoor")
    p.add_argument("--window", nargs=2, type=float, metavar=("MIN", "MAX"),
                   help="custom band-gap window in eV (overrides --mode window)")
    p.add_argument("--radius", action="append", default=[], metavar="ION=R:CHARGE",
                   help="override/add an ionic radius, e.g. --radius Sn=1.18:2")
    args = p.parse_args(argv)

    radii = {}
    for item in args.radius:
        try:
            ion, rest = item.split("=")
            r, q = rest.split(":")
            radii[ion] = (float(r), int(q))
        except ValueError:
            p.error(f"bad --radius value {item!r}; expected ION=R:CHARGE")

    df = pd.read_csv(args.csv, comment="#", skipinitialspace=True)
    ranked = rank_candidates(df, mode=args.mode,
                             window=tuple(args.window) if args.window else None,
                             radii=radii or None)
    if args.output:
        ranked.to_csv(args.output, index=False)
        print(f"wrote {len(ranked)} rows to {args.output}")
    else:
        ranked.to_csv(sys.stdout, index=False)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
