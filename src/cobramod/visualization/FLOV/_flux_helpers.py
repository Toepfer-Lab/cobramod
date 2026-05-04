"""Stateless helpers used by the FLoV widget.

Includes model classification (compartment, currency, proton, reaction kind),
edge/flux/colour mapping, spread parsing, hull-path rendering, station
metabolite-class detection, pair-key resolution, and CSV loaders.

Split out of ``flux_network.py`` so the main file holds the widget class only.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Optional, Union

import cobra
import numpy as np
import pandas as pd

from ._flux_config import (
    ACTIVE_EDGE_EPS,
    EDGE_BINS,
    EDGE_LOG_STRETCH,
    EDGE_MIN_VIS_NORM,
    _AA_CODES,
    _COFAC_KW,
    _INORG_KW,
    _SUGAR_KW,
)


# ==========================================================================
# MODEL CLASSIFICATION HELPERS (config-aware, stateless)
# ==========================================================================

def _met_comp(met: cobra.Metabolite, cfg: dict) -> str:
    if met.compartment:
        key = met.compartment.lower()
        if key in cfg["compartments"]:
            return key
    return cfg["exchange_key"]


def _is_currency(met: cobra.Metabolite, cfg: dict) -> bool:
    t = (met.name or "").lower() + " " + met.id.lower()
    return any(kw in t for kw in cfg["currency_kw"])


def _normalise_met_token(raw: str, cfg: dict) -> str:
    """Strip compartment suffix/bracket, collapse repeated underscores, and
    remove any residual trailing token.  Used by _is_proton on both the
    metabolite id and its name (TAIR10 stores names as 'H_plus__c0')."""
    s = cfg["bracket_re"].sub("", raw).lower().strip()
    s = cfg["suffix_re"].sub("", s).strip("_")
    s = re.sub(r"_+", "_", s)
    s = re.sub(r"_[a-z0-9]+$", "", s)
    if s.startswith("m_"):
        s = s[2:]
    return s


def _is_proton(met: cobra.Metabolite, cfg: dict) -> bool:
    name_lower = (met.name or "").lower().strip()
    id_base   = _normalise_met_token(met.id, cfg)
    name_base = _normalise_met_token(name_lower, cfg)
    _proton_tokens = {"h", "h_plus", "hplus"}
    return (
        id_base   in _proton_tokens
        or name_base in _proton_tokens
        or name_lower in {"h+", "proton", "h"}
        or "proton" in name_lower
    )


def _rxn_kind(rxn: cobra.Reaction, cfg: dict) -> str:
    rxn_id_lower = rxn.id.lower()
    if any(rxn_id_lower.startswith(p.lower()) for p in cfg["import_prefix"]):
        return "import"
    if any(rxn_id_lower.startswith(p.lower()) for p in cfg["export_prefix"]):
        return "export"
    if rxn_id_lower.startswith("biomass"):
        return "biomass"
    if any(rxn_id_lower.startswith(p.lower()) for p in cfg["biomass_prefix"]):
        return "biomass"
    all_comps = {_met_comp(m, cfg) for m in rxn.metabolites}
    non_exch = all_comps - {cfg["exchange_key"]}
    if non_exch and len(all_comps) > 1:
        return "transport"
    if len(rxn.metabolites) == 1 and non_exch:
        return "transport"
    return "regular"


# ==========================================================================
# MISC SCALAR HELPERS
# ==========================================================================

def _stable_unit(met_id: str) -> float:
    digest = hashlib.md5(met_id.encode("utf-8")).hexdigest()
    v = int(digest[:8], 16) / float(0xFFFFFFFF)
    return 2.0 * v - 1.0


def _has_flux_data(flux: float) -> bool:
    return (not np.isnan(flux)) and (abs(float(flux)) > ACTIVE_EDGE_EPS)


def _emit_progress(message: str) -> None:
    print(f"FLoV: {message}", flush=True)


def _parse_spread(
    value: Union[float, int, str],
) -> tuple[float, float, float, float]:
    """Parse a spread specification into ``(u, d, l, r)`` scale factors.

    Accepts:

    * A plain number – applies uniformly in all four directions.
    * A string of one or more ``<number><direction>`` tokens (case-insensitive,
      optional whitespace between tokens).  Direction letters are
      ``u`` (up / +y), ``d`` (down / -y), ``l`` (left / -x),
      ``r`` (right / +x).  Any direction not mentioned defaults to ``1.0``
      (no change).  Examples::

          "2u"          # double distance upward, keep everything else
          "2u1.5r"      # spread up ×2, right ×1.5
          "2u 0.8d 1.5r 1.2l"  # all four directions independently
    """
    if isinstance(value, (int, float)):
        s = float(value)
        return s, s, s, s
    s = str(value).strip()
    try:
        v = float(s)
        return v, v, v, v
    except ValueError:
        pass
    tokens = re.findall(r'([0-9]*\.?[0-9]+)\s*([udlrUDLR])', s)
    if not tokens:
        raise ValueError(
            f"Invalid spread value {value!r}. "
            "Use a number or tokens like '2u', '1.5r', '2u1.5l'."
        )
    result: dict[str, float] = {"u": 1.0, "d": 1.0, "l": 1.0, "r": 1.0}
    for num_str, direction in tokens:
        result[direction.lower()] = float(num_str)
    return result["u"], result["d"], result["l"], result["r"]


# ==========================================================================
# EDGE / FLUX STYLE HELPERS
# ==========================================================================

def _edge_bucket(flux: float, abs_max: float) -> tuple[int, float, int]:
    abs_flux = abs(flux)
    if abs_max > 1e-9:
        rel = abs_flux / abs_max
        mag_norm = np.log1p(EDGE_LOG_STRETCH * rel) / np.log1p(EDGE_LOG_STRETCH)
        if abs_flux > 0.0:
            mag_norm = EDGE_MIN_VIS_NORM + (1.0 - EDGE_MIN_VIS_NORM) * mag_norm
    else:
        mag_norm = 0.0
    b = int(np.clip(np.floor(mag_norm * EDGE_BINS), 0, EDGE_BINS - 1))
    sign = 1 if flux >= 0 else -1
    return b, float(mag_norm), sign


def _edge_rgba(mag_norm: float, sign: int, has_flux: bool, density_scale: float = 1.0, color_modifier: float = 1.0) -> str:
    if not has_flux:
        return "rgba(110,118,129,0.20)"
    alpha = (0.55 + 0.45 * mag_norm) * (0.60 + 0.40 * density_scale)
    eff_norm = float(np.clip(mag_norm * color_modifier, 0.0, 1.0))
    if sign > 0:
        r = int(255)
        g = int(180 - 180 * eff_norm)
        b_ = int(70 - 70 * eff_norm)
    else:
        r = int(100 - 50 * eff_norm)
        g = int(190 - 90 * eff_norm)
        b_ = int(255)
    return f"rgba({r},{g},{b_},{alpha:.3f})"


def _flux_to_hex(flux: float, abs_max: float) -> str:
    """Map a flux value to a dark-mode hex colour using a diverging red/blue scale."""
    if abs_max < 1e-9 or not _has_flux_data(flux):
        return "#6e7681"
    t = float(np.clip(abs(flux) / abs_max, 0.0, 1.0))
    if flux >= 0:
        r = int(239 + (198 - 239) * t)
        g = int(154 + (40 - 154) * t)
        b_ = int(154 + (40 - 154) * t)
    else:
        r = int(144 + (21 - 144) * t)
        g = int(202 + (101 - 202) * t)
        b_ = int(249 + (192 - 249) * t)
    return f"#{r:02x}{g:02x}{b_:02x}"


# ==========================================================================
# COLOUR HELPERS
# ==========================================================================

def _tint_hex(hex_color: str, alpha: float = 0.15) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _modulate_route_color(base_hex: str, mag_norm: float, has_flux: bool) -> str:
    """Return an rgba color keeping the base hue; opacity scaled by flux magnitude."""
    if not has_flux:
        return "rgba(110,118,129,0.18)"
    h = base_hex.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    alpha = 0.42 + 0.58 * float(mag_norm)
    return f"rgba({r},{g},{b},{alpha:.3f})"


# ==========================================================================
# SHAPE HELPERS
# ==========================================================================

def _smooth_hull_path(verts: np.ndarray, tension: float = 0.4) -> str:
    n = len(verts)
    parts = [f"M {verts[0][0]:.4f},{verts[0][1]:.4f}"]
    for i in range(n):
        p0 = verts[(i - 1) % n]
        p1 = verts[i % n]
        p2 = verts[(i + 1) % n]
        p3 = verts[(i + 2) % n]
        cp1 = p1 + tension * (p2 - p0) / 3.0
        cp2 = p2 - tension * (p3 - p1) / 3.0
        parts.append(
            f"C {cp1[0]:.4f},{cp1[1]:.4f} "
            f"{cp2[0]:.4f},{cp2[1]:.4f} "
            f"{p2[0]:.4f},{p2[1]:.4f}"
        )
    parts.append("Z")
    return " ".join(parts)


# ==========================================================================
# STATION METABOLITE-CLASS / PAIR-KEY HELPERS
# ==========================================================================

def _met_class(rxn: cobra.Reaction, cfg: dict) -> str:
    """Classify a transport/exchange reaction by its dominant non-proton
    metabolite. Cheap, static, no DB lookups — keep it fast.
    """
    cands = [
        (abs(s), m) for m, s in rxn.metabolites.items()
        if not _is_proton(m, cfg)
    ]
    if not cands:
        return "other"
    cands.sort(key=lambda t: t[0], reverse=True)
    primary = cands[0][1]
    name = (primary.name or "").lower()
    base = _normalise_met_token(primary.id, cfg)
    for aa in _AA_CODES:
        if base == aa or base.startswith(aa + "_") or base.startswith(aa + "-"):
            return "amino_acid"
    name_padded = f" {name} "
    for aa in _AA_CODES:
        if f" {aa} " in name_padded:
            return "amino_acid"
    if any(kw in base or kw in name for kw in _SUGAR_KW):
        return "sugar"
    if any(kw in base for kw in _COFAC_KW):
        return "cofactor"
    for kw in _INORG_KW:
        if base == kw or base.startswith(kw + "_"):
            return "inorganic"
    return "other"


def _pair_key(
    kind: str, non_exch_comps: set[str], cfg: dict,
) -> Optional[tuple[str, str]]:
    """Canonical (compA, compB) pair for a non-regular reaction.

    For ``transport``: the two non-exchange compartments it touches
    (lex-sorted). For ``import``/``export``/``biomass``: ``(comp, exchange)``.
    Returns ``None`` if no compartment can be resolved (the user's
    "ignore reactions that lead to no compartment" rule).
    """
    exch = cfg["exchange_key"]
    comps = sorted(non_exch_comps - {exch})
    if kind in ("import", "export", "biomass"):
        return (comps[0], exch) if comps else None
    if kind == "transport":
        if len(comps) >= 2:
            return (comps[0], comps[1])
        return (comps[0], exch) if comps else None
    return None


# ==========================================================================
# CSV LOADERS
# ==========================================================================

def _load_ignore_csv(path: Union[str, Path]) -> tuple[set[str], set[str]]:
    """Parse a CSV listing reaction and metabolite IDs to exclude from the graph.

    Expected format — two columns, order does not matter, extra columns ignored::

        reactions,metabolites
        RXN_001,cpd00067_c0
        RXN_002,cpd00009_m0
        ,cpd00012_d0

    Empty cells are skipped.  Column headers are matched case-insensitively;
    accepted names are ``reactions``/``rxn``/``rxn_id`` and
    ``metabolites``/``met``/``met_id``.

    Returns
    -------
    rxn_ids : set[str]
        Reaction IDs to ignore.
    met_ids : set[str]
        Metabolite IDs to ignore.
    """
    df = pd.read_csv(Path(path), dtype=str)
    df.columns = df.columns.str.strip().str.lower()

    _rxn_kw = ["reactions", "rxn", "rxn_id", "reaction"]
    _met_kw = ["metabolites", "met", "met_id", "metabolite"]

    rxn_col = next((c for c in _rxn_kw if c in df.columns), None)
    met_col = next((c for c in _met_kw if c in df.columns), None)

    rxn_ids: set[str] = set()
    met_ids: set[str] = set()

    if rxn_col:
        rxn_ids = {v.strip() for v in df[rxn_col].dropna() if v.strip()}
    if met_col:
        met_ids = {v.strip() for v in df[met_col].dropna() if v.strip()}

    return rxn_ids, met_ids


def _load_single_csv(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    _id_kw = ["reaction_id", "rxn_id", "rxn", "id", "reaction"]
    _val_kw = ["flux", "value", "fba_flux", "v"]
    id_col = next((c for c in _id_kw if c in df.columns), df.columns[0])
    val_col = next((c for c in _val_kw if c in df.columns), df.columns[1])
    return pd.Series(
        pd.to_numeric(df[val_col], errors="coerce").values,
        index=df[id_col].astype(str).str.strip(),
    )
