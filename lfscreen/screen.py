"""Screening filters, a simple heuristic score, and an optional MP fetcher.

The score is a transparent heuristic ranking aid, NOT a detailed-balance
efficiency estimate.
"""

from __future__ import annotations

import math
import os
from typing import Iterable, Mapping, Optional, Tuple

import pandas as pd

from .descriptors import describe, parse_composition

# Approximate band-gap windows (eV). Outdoor: AM1.5G single junction optimum
# region; indoor: white LED / fluorescent spectra favour wider gaps.
BANDGAP_WINDOWS = {
    "outdoor": (1.1, 1.6),
    "indoor": (1.7, 2.0),
}

TOXIC_ELEMENTS = frozenset({"Pb", "Cd", "Hg", "Tl", "As"})

# Heuristic windows for the geometric descriptors.
TOLERANCE_RANGE = (0.8, 1.0)
OCTAHEDRAL_RANGE = (0.414, 0.9)

WEIGHTS = {"band_gap": 0.4, "tau": 0.3, "tolerance": 0.15, "octahedral": 0.15}
BANDGAP_SOFTNESS = 0.3  # eV; linear score decay outside the window


def is_toxic(formula: str, toxic: Iterable[str] = TOXIC_ELEMENTS) -> bool:
    """True if the composition contains any element from ``toxic``."""
    return bool(parse_composition(formula).elements & set(toxic))


def bandgap_score(eg: float, window: Tuple[float, float],
                  softness: float = BANDGAP_SOFTNESS) -> float:
    """1 inside the window, decaying linearly to 0 at ``softness`` eV outside."""
    if eg is None or (isinstance(eg, float) and math.isnan(eg)):
        return math.nan
    lo, hi = window
    d = max(lo - eg, 0.0, eg - hi)
    return max(0.0, 1.0 - d / softness) if softness > 0 else float(d == 0)


def _in(x: float, rng: Tuple[float, float]) -> bool:
    return not math.isnan(x) and rng[0] <= x <= rng[1]


def screen_candidates(df: pd.DataFrame, mode: str = "outdoor",
                      window: Optional[Tuple[float, float]] = None,
                      radii: Optional[Mapping[str, tuple]] = None,
                      formula_col: str = "formula",
                      bandgap_col: str = "band_gap") -> pd.DataFrame:
    """Add descriptor, filter and score columns to a candidate table.

    Score = weighted mean of the applicable criteria (band gap, tau < 4.18,
    tolerance factor in range, octahedral factor in range). Structural
    criteria are skipped (not penalised) for non-perovskite A3B2X9 entries,
    which are flagged via ``family``. Toxic compositions score 0.
    """
    if window is None:
        if mode not in BANDGAP_WINDOWS:
            raise ValueError(f"mode must be one of {sorted(BANDGAP_WINDOWS)}")
        window = BANDGAP_WINDOWS[mode]
    if formula_col not in df.columns:
        raise KeyError(f"missing column {formula_col!r}")

    rows = []
    for _, rec in df.iterrows():
        formula = str(rec[formula_col]).strip()
        d = describe(formula, radii=radii)
        eg = float(rec[bandgap_col]) if bandgap_col in df.columns and pd.notna(rec[bandgap_col]) else math.nan
        d["toxic"] = is_toxic(formula)
        d["band_gap_ok"] = _in(eg, window)

        parts = {}
        bg = bandgap_score(eg, window)
        if not math.isnan(bg):
            parts["band_gap"] = bg
        if d["tau_perovskite"] is not None:
            parts["tau"] = float(d["tau_perovskite"])
        if not math.isnan(d["tolerance_factor"]):
            parts["tolerance"] = float(_in(d["tolerance_factor"], TOLERANCE_RANGE))
            parts["octahedral"] = float(_in(d["octahedral_factor"], OCTAHEDRAL_RANGE))
        wsum = sum(WEIGHTS[k] for k in parts)
        score = sum(WEIGHTS[k] * v for k, v in parts.items()) / wsum if wsum else math.nan
        d["score"] = 0.0 if d["toxic"] else round(score, 4)
        d["passes"] = bool(not d["toxic"] and d["band_gap_ok"]
                           and d["tau_perovskite"] is not False)
        d.pop("formula")
        rows.append(d)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    out["mode"] = mode
    return out


def rank_candidates(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Screen and sort by (passes, score) descending."""
    out = screen_candidates(df, **kwargs)
    return out.sort_values(["passes", "score"], ascending=False,
                           na_position="last").reset_index(drop=True)


def fetch_materials_project(formulas: Iterable[str],
                            api_key: Optional[str] = None) -> pd.DataFrame:
    """Fetch lowest-energy-above-hull band gaps from the Materials Project.

    Requires the optional ``mp-api`` package and an API key (argument or the
    ``MP_API_KEY`` environment variable). Note that DFT (GGA/PBE) band gaps
    are typically underestimated.
    """
    try:
        from mp_api.client import MPRester  # lazy optional import
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError("install the optional dependency: pip install mp-api") from exc
    key = api_key or os.environ.get("MP_API_KEY")
    if not key:
        raise RuntimeError("set MP_API_KEY or pass api_key=")

    records = []
    with MPRester(key) as mpr:  # pragma: no cover - network
        for f in formulas:
            reduced = f.replace("MA", "CH6N").replace("FA", "CH5N2")
            docs = mpr.materials.summary.search(
                formula=reduced,
                fields=["material_id", "formula_pretty", "band_gap", "energy_above_hull"],
            )
            if not docs:
                records.append({"formula": f, "material_id": None,
                                "band_gap": math.nan, "energy_above_hull": math.nan})
                continue
            best = min(docs, key=lambda d: d.energy_above_hull)
            records.append({"formula": f, "material_id": str(best.material_id),
                            "band_gap": best.band_gap,
                            "energy_above_hull": best.energy_above_hull})
    return pd.DataFrame(records)
