"""Composition parsing and ionic-radius structural descriptors.

Radii are in angstrom. A-site values use 12-fold coordination, B-site values
6-fold coordination and halides 6-fold coordination (Shannon, 1976), except the
molecular cations MA+ and FA+, for which commonly used *effective* radii are
given. Always verify values for the coordination numbers relevant to your
system; pass ``radii=`` to override or extend the table.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional

# ion -> (radius in angstrom, formal charge)
IONIC_RADII: Dict[str, tuple] = {
    # A-site cations (CN 12; MA/FA are effective radii for molecular cations)
    "Cs": (1.88, 1),
    "Rb": (1.72, 1),
    "K": (1.64, 1),
    "MA": (2.17, 1),  # CH3NH3+ effective radius
    "FA": (2.53, 1),  # CH(NH2)2+ effective radius
    # B-site cations (CN 6)
    "Pb": (1.19, 2),  # for comparison only
    "Sn": (1.10, 2),  # value widely used in the perovskite literature
    "Ge": (0.73, 2),
    "Bi": (1.03, 3),
    "Sb": (0.76, 3),
    "Ag": (1.15, 1),
    # X-site anions (CN 6)
    "Cl": (1.81, -1),
    "Br": (1.96, -1),
    "I": (2.20, -1),
}

A_SITE = {"Cs", "Rb", "K", "MA", "FA"}
X_SITE = {"F", "Cl", "Br", "I"}
ORGANIC = {"MA", "FA"}

TAU_THRESHOLD = 4.18  # Bartel et al. (2019): tau < 4.18 -> perovskite predicted

_TOKEN = re.compile(r"(MA|FA|[A-Z][a-z]?)(\d*\.?\d*)")


@dataclass
class Composition:
    """Parsed composition with site assignment."""

    formula: str
    counts: Dict[str, float]
    a_site: Dict[str, float] = field(default_factory=dict)
    b_site: Dict[str, float] = field(default_factory=dict)
    x_site: Dict[str, float] = field(default_factory=dict)
    family: str = "unknown"  # ABX3 | A2BB'X6 | A3B2X9 | unknown

    @property
    def elements(self) -> set:
        """Chemical elements present (organic cations expanded to C, H, N)."""
        out = {k for k in self.counts if k not in ORGANIC}
        if any(k in ORGANIC for k in self.counts):
            out |= {"C", "H", "N"}
        return out

    @property
    def is_perovskite_family(self) -> bool:
        return self.family in ("ABX3", "A2BB'X6")

    @property
    def has_organic(self) -> bool:
        return any(k in ORGANIC for k in self.counts)


def parse_composition(formula: str) -> Composition:
    """Parse formulas such as ``CsSnI3``, ``Cs2AgBiBr6`` or ``FA3Bi2Br9``.

    ``MA`` and ``FA`` are treated as single organic A-site cations. Mixed
    halides (e.g. ``CsSnI2Br``) are supported.
    """
    formula = formula.strip()
    if not formula:
        raise ValueError("empty formula")
    counts: Dict[str, float] = {}
    pos = 0
    for m in _TOKEN.finditer(formula):
        if m.start() != pos:
            raise ValueError(f"cannot parse formula {formula!r} at position {pos}")
        sym, num = m.group(1), m.group(2)
        n = float(num) if num else 1.0
        counts[sym] = counts.get(sym, 0.0) + n
        pos = m.end()
    if pos != len(formula):
        raise ValueError(f"cannot parse formula {formula!r} at position {pos}")

    comp = Composition(formula=formula, counts=counts)
    comp.a_site = {k: v for k, v in counts.items() if k in A_SITE}
    comp.x_site = {k: v for k, v in counts.items() if k in X_SITE}
    comp.b_site = {k: v for k, v in counts.items() if k not in A_SITE and k not in X_SITE}

    na, nb, nx = (sum(d.values()) for d in (comp.a_site, comp.b_site, comp.x_site))
    if na and nb and nx:
        if math.isclose(na / nx, 1 / 3) and math.isclose(nb / nx, 1 / 3):
            # Same site ratios for ABX3 and A2BB'X6: written as A2 B B' X6 with
            # two distinct B species -> ordered double perovskite; otherwise
            # treated as a (possibly alloyed) ABX3 with averaged B radius.
            if len(comp.b_site) == 2 and math.isclose(nx, 6.0) and math.isclose(na, 2.0):
                comp.family = "A2BB'X6"
            else:
                comp.family = "ABX3"
        elif math.isclose(na / nx, 3 / 9) and math.isclose(nb / nx, 2 / 9):
            comp.family = "A3B2X9"
    return comp


def _radius(ion: str, radii: Mapping[str, tuple]) -> float:
    try:
        return float(radii[ion][0])
    except KeyError:
        raise KeyError(
            f"no ionic radius for {ion!r}; pass radii={{'{ion}': (r, charge)}}"
        ) from None


def _charge(ion: str, radii: Mapping[str, tuple]) -> int:
    try:
        return int(radii[ion][1])
    except KeyError:
        raise KeyError(f"no charge for {ion!r}; pass radii={{'{ion}': (r, charge)}}") from None


def _weighted(site: Mapping[str, float], radii: Mapping[str, tuple]) -> float:
    total = sum(site.values())
    return sum(_radius(k, radii) * n for k, n in site.items()) / total


def tolerance_factor(r_a: float, r_b: float, r_x: float) -> float:
    """Goldschmidt (1926) tolerance factor t = (rA + rX) / (sqrt(2) (rB + rX))."""
    return (r_a + r_x) / (math.sqrt(2.0) * (r_b + r_x))


def octahedral_factor(r_b: float, r_x: float) -> float:
    """Octahedral factor mu = rB / rX."""
    return r_b / r_x


def bartel_tau(r_a: float, r_b: float, r_x: float, n_a: float) -> float:
    """Bartel et al. (2019) tau = rX/rB - nA (nA - (rA/rB) / ln(rA/rB)).

    Requires rA > rB (as in the original formulation). tau < 4.18 predicts a
    perovskite.
    """
    ratio = r_a / r_b
    if ratio <= 1.0:
        raise ValueError("tau is defined for rA > rB")
    return r_x / r_b - n_a * (n_a - ratio / math.log(ratio))


def describe(formula: str, radii: Optional[Mapping[str, tuple]] = None) -> dict:
    """Return composition info and (where applicable) structural descriptors."""
    table = dict(IONIC_RADII)
    if radii:
        table.update(radii)
    comp = parse_composition(formula)
    out = {
        "formula": formula,
        "family": comp.family,
        "organic": comp.has_organic,
        "tolerance_factor": math.nan,
        "octahedral_factor": math.nan,
        "tau": math.nan,
        "tau_perovskite": None,
        "note": "",
    }
    if not comp.is_perovskite_family:
        out["note"] = (
            "non-perovskite (A3B2X9) structure; descriptors not applied"
            if comp.family == "A3B2X9" else "unrecognised stoichiometry"
        )
        return out
    try:
        r_a = _weighted(comp.a_site, table)
        r_b = _weighted(comp.b_site, table)  # average B radius for double perovskites
        r_x = _weighted(comp.x_site, table)
        n_a = sum(_charge(k, table) * n for k, n in comp.a_site.items()) / sum(comp.a_site.values())
    except KeyError as exc:
        out["note"] = str(exc.args[0])
        return out
    out["tolerance_factor"] = tolerance_factor(r_a, r_b, r_x)
    out["octahedral_factor"] = octahedral_factor(r_b, r_x)
    try:
        tau = bartel_tau(r_a, r_b, r_x, n_a)
        out["tau"] = tau
        out["tau_perovskite"] = bool(tau < TAU_THRESHOLD)
    except ValueError as exc:
        out["note"] = str(exc)
    return out
