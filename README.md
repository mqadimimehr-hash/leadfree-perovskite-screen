# leadfree-perovskite-screen

A small, dependency-light Python toolkit (`lfscreen`) for screening lead-free
perovskite and perovskite-inspired absorbers for photovoltaics, including
indoor PV.

## Features

- **Composition parsing**: `CsSnI3`, `Cs2AgBiBr6`, `FA3Bi2Br9`, and mixed
  halides such as `CsSnI2Br`. `MA` (methylammonium) and `FA` (formamidinium)
  are treated as organic A-site cations.
- **Structural descriptors from ionic radii**
  - Goldschmidt tolerance factor, `t = (rA + rX) / (sqrt(2) (rB + rX))` (Goldschmidt, 1926)
  - octahedral factor, `mu = rB / rX`
  - tau factor, `tau = rX/rB - nA (nA - (rA/rB) / ln(rA/rB))`, where
    `tau < 4.18` predicts a perovskite (Bartel et al., *Science Advances*, 2019)
  - ABX3 and A2BB'X6 double perovskites (the average B-site radius is used)
  - A3B2X9 compositions are flagged as **non-perovskite**. Only composition
    information is returned for them.
- **Toxicity flag**: any composition containing Pb, Cd, Hg, Tl or As fails.
- **Band-gap window filter**: you supply the band gaps in a CSV.
- **Simple heuristic score** that combines the filters. It is *not* a
  detailed-balance efficiency estimate.
- **Optional Materials Project fetcher** (`mp-api`, loaded only when used).
- **CLI**: `lfscreen`.

## Installation

```bash
pip install -e .            # core (numpy, pandas)
pip install -e ".[mp]"      # optional Materials Project support
pip install -e ".[test]"    # pytest
```

## Usage

### CLI

```bash
lfscreen examples/candidates.csv -o ranked.csv --mode indoor
lfscreen examples/candidates.csv --mode outdoor
lfscreen examples/candidates.csv --window 1.8 2.1 --radius Sn=1.18:2
```

The input CSV needs a `formula` column and, optionally, a `band_gap` column
(eV). Lines that start with `#` are ignored.

### Python

```python
import pandas as pd
from lfscreen import describe, rank_candidates

describe("Cs2AgBiBr6")
# {'family': "A2BB'X6", 'tolerance_factor': ..., 'octahedral_factor': ..., 'tau': ..., ...}

df = pd.read_csv("examples/candidates.csv", comment="#")
ranked = rank_candidates(df, mode="indoor", radii={"Sn": (1.18, 2)})
```

### Materials Project (optional)

```python
import os
from lfscreen import fetch_materials_project
os.environ["MP_API_KEY"] = "..."   # or pass api_key=
fetch_materials_project(["CsSnI3", "Cs2AgBiBr6"])
```

DFT (GGA/PBE) band gaps are usually underestimated. Correct them or check them
before you use them for screening.

## Band-gap windows (approximate)

| mode      | window (eV) | rationale                                      |
|-----------|-------------|------------------------------------------------|
| `outdoor` | 1.1 – 1.6   | near the single-junction optimum under AM1.5G  |
| `indoor`  | 1.7 – 2.0   | typical white-LED / fluorescent indoor spectra |

These windows are rough guides, not strict limits. Use `--window` or
`window=` to set your own.

## Scoring

Each criterion gives a value between 0 and 1:

- **Band gap** (weight 0.4): 1 inside the window, falling linearly to 0 at
  0.3 eV outside it.
- **tau < 4.18** (weight 0.3).
- **0.8 <= t <= 1.0** (weight 0.15).
- **0.414 <= mu <= 0.9** (weight 0.15).

The score is the weighted mean of the criteria that apply. For A3B2X9
entries, the structural criteria are skipped rather than counted as failures.
Toxic compositions score 0. `passes` is true when a composition is non-toxic,
its band gap is inside the window, and tau does not predict a non-perovskite.
Treat the score as a ranking aid only.

## Ionic radii: please verify

The built-in table (`lfscreen.IONIC_RADII`, in angstrom) holds a small set of
commonly used values. Most are Shannon radii (CN 12 for A-site, CN 6 for
B-site and X-site ions). The table also has commonly used effective radii for
MA+ and FA+. Some values (for example Sn2+ in 6-fold coordination) differ
between sources. **Check every radius for the coordination number that fits
your system.** You can override or extend the table with
`radii={"Ion": (radius, charge)}` or with `--radius Ion=R:CHARGE` on the
command line.

## Example data

`examples/candidates.csv` lists a few well-known compositions (Sn-, Ge-, Bi-
and Sb-based, plus Pb references for comparison). **Its `band_gap` values are
rough, approximate literature figures for demonstration only. Check them
against primary sources before any real use.**

## Tests

```bash
pytest
```

## References

- V. M. Goldschmidt, 1926 (tolerance factor).
- C. J. Bartel et al., *Science Advances*, 2019 (tau factor).

## License

MIT © 2026 Mohammad Ghadimimehr ([mqadimimehr-hash](https://github.com/mqadimimehr-hash))
