import math
from pathlib import Path

import pandas as pd
import pytest

from lfscreen import (
    bartel_tau, describe, is_toxic, octahedral_factor, parse_composition,
    rank_candidates, screen_candidates, tolerance_factor,
)
from lfscreen.cli import main

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "candidates.csv"


def test_cspbi3_tolerance_in_expected_range():
    t = describe("CsPbI3")["tolerance_factor"]
    assert 0.8 < t < 0.9
    assert t == pytest.approx((1.88 + 2.20) / (math.sqrt(2) * (1.19 + 2.20)))


def test_tau_matches_hand_calculation():
    ra, rb, rx, na = 1.88, 1.19, 2.20, 1
    # hand calculation: rX/rB = 1.848739..., rA/rB = 1.579832..., ln = 0.457318...
    expected = 2.20 / 1.19 - 1 * (1 - (1.88 / 1.19) / math.log(1.88 / 1.19))
    assert bartel_tau(ra, rb, rx, na) == pytest.approx(expected)
    assert bartel_tau(ra, rb, rx, na) == pytest.approx(4.3032, abs=1e-3)
    assert describe("CsPbI3")["tau"] == pytest.approx(expected)


def test_tau_requires_ra_gt_rb():
    with pytest.raises(ValueError):
        bartel_tau(1.0, 1.1, 2.2, 1)


def test_octahedral_factor():
    assert octahedral_factor(1.10, 2.20) == pytest.approx(0.5)
    assert tolerance_factor(1.0, 1.0, 1.0) == pytest.approx(1 / math.sqrt(2))


@pytest.mark.parametrize("formula,family,organic", [
    ("CsSnI3", "ABX3", False),
    ("MAPbI3", "ABX3", True),
    ("FASnI3", "ABX3", True),
    ("Cs2AgBiBr6", "A2BB'X6", False),
    ("FA3Bi2Br9", "A3B2X9", True),
    ("CsSnI2Br", "ABX3", False),
])
def test_parse_families(formula, family, organic):
    c = parse_composition(formula)
    assert c.family == family
    assert c.has_organic is organic


def test_parse_counts_and_errors():
    c = parse_composition("FA3Bi2Br9")
    assert c.counts == {"FA": 3, "Bi": 2, "Br": 9}
    assert "N" in c.elements and "FA" not in c.elements
    with pytest.raises(ValueError):
        parse_composition("Cs-SnI3")
    with pytest.raises(ValueError):
        parse_composition("")


def test_double_perovskite_uses_average_b_radius():
    d = describe("Cs2AgBiBr6")
    rb = (1.15 + 1.03) / 2
    assert d["octahedral_factor"] == pytest.approx(rb / 1.96)
    assert d["tolerance_factor"] == pytest.approx((1.88 + 1.96) / (math.sqrt(2) * (rb + 1.96)))


def test_a3b2x9_flagged_non_perovskite():
    d = describe("FA3Bi2Br9")
    assert d["family"] == "A3B2X9"
    assert math.isnan(d["tolerance_factor"]) and math.isnan(d["tau"])
    assert "non-perovskite" in d["note"]


def test_radius_override_and_missing():
    base = describe("CsSnI3")["tolerance_factor"]
    over = describe("CsSnI3", radii={"Sn": (1.18, 2)})["tolerance_factor"]
    assert over < base
    d = describe("CsCuI3")
    assert math.isnan(d["tolerance_factor"]) and "Cu" in d["note"]


def test_toxicity():
    assert is_toxic("MAPbI3")
    assert is_toxic("Cs2TlBiBr6")
    assert not is_toxic("Cs2AgBiBr6")


def test_screen_modes_and_toxic_zero():
    df = pd.DataFrame({"formula": ["CsSnI3", "Cs2AgBiBr6", "CsPbI3", "Cs3Bi2I9"],
                       "band_gap": [1.3, 1.9, 1.7, 2.1]})
    out = screen_candidates(df, mode="indoor").set_index("formula")
    assert out.loc["CsPbI3", "score"] == 0 and not out.loc["CsPbI3", "passes"]
    assert out.loc["Cs2AgBiBr6", "band_gap_ok"]
    assert not out.loc["CsSnI3", "band_gap_ok"]
    assert out.loc["Cs3Bi2I9", "family"] == "A3B2X9"
    outdoor = screen_candidates(df, mode="outdoor").set_index("formula")
    assert outdoor.loc["CsSnI3", "band_gap_ok"]
    with pytest.raises(ValueError):
        screen_candidates(df, mode="space")


def test_rank_and_cli(tmp_path):
    df = pd.read_csv(EXAMPLE, comment="#")
    ranked = rank_candidates(df, mode="outdoor")
    assert ranked.iloc[0]["passes"]
    assert not ranked[ranked.formula == "MAPbI3"].iloc[0]["passes"]
    out = tmp_path / "ranked.csv"
    assert main([str(EXAMPLE), "-o", str(out), "--mode", "indoor"]) == 0
    res = pd.read_csv(out)
    assert len(res) == len(df) and "score" in res.columns
    assert (res["mode"] == "indoor").all()
