"""lfscreen: screening of lead-free perovskite and perovskite-inspired PV absorbers."""

from .descriptors import (
    IONIC_RADII,
    Composition,
    bartel_tau,
    describe,
    octahedral_factor,
    parse_composition,
    tolerance_factor,
)
from .screen import (
    BANDGAP_WINDOWS,
    fetch_materials_project,
    is_toxic,
    rank_candidates,
    screen_candidates,
)

__version__ = "0.1.0"

__all__ = [
    "IONIC_RADII", "Composition", "bartel_tau", "describe", "octahedral_factor",
    "parse_composition", "tolerance_factor", "BANDGAP_WINDOWS",
    "fetch_materials_project", "is_toxic", "rank_candidates", "screen_candidates",
]
