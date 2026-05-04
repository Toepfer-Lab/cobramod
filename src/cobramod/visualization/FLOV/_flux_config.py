"""Constants, dataclasses, defaults, and TOML config loaders for the FLoV widget.

Split out of ``flux_network.py`` so the main file holds the widget class only.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ImportError:
        tomllib = None  # type: ignore[assignment]


# ==========================================================================
# CONSTANTS
# ==========================================================================

ACTIVE_EDGE_EPS = 1e-6
EDGE_BINS = 6
_MET_COLOR = "#58a6ff"
EDGE_LOG_STRETCH = 99.0
EDGE_MIN_VIS_NORM = 0.22


# ==========================================================================
# DATA CLASSES
# ==========================================================================

@dataclass
class ViewSpec:
    """One flux view shown in the visualisation."""
    label: str
    flux_dict: dict[str, float]
    std_dict: dict[str, float] = field(default_factory=dict)
    hover_extra: dict[str, str] = field(default_factory=dict)
    pipeline_tag: str = ""


# ==========================================================================
# DEFAULTS  (AraCore layout — overridden by TOML config)
# ==========================================================================

_DEFAULT_COMPARTMENTS: dict[str, dict] = {
    "h": dict(region=(-11.5, -2.5, -4.5, 6.0), color="#3fb950",
              fill="rgba(63,185,80,0.15)", label="Chloroplast"),
    "c": dict(region=(-2.0, 2.0, -4.5, 6.0), color="#58a6ff",
              fill="rgba(88,166,255,0.15)", label="Cytosol"),
    "m": dict(region=(2.5, 7.5, 0.0, 6.0), color="#ffa657",
              fill="rgba(255,166,87,0.15)", label="Mitochondria"),
    "p": dict(region=(2.5, 6.5, -4.5, -0.5), color="#d2a8ff",
              fill="rgba(210,168,255,0.15)", label="Peroxisome"),
    "u": dict(region=(-11.5, 7.5, -8.5, -5.0), color="#8b949e",
              fill="rgba(139,148,158,0.10)", label="Exchange / Transport"),
}

_DEFAULT_CURRENCY_KW: set[str] = {
    "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
    "water", "h2o", "proton", "h+", "orthophosphate", "pyrophosphate",
    "phosphate", "co2", "carbon dioxide", "oxygen", " o2",
    "coa", "coenzyme a", "fad", "fadh2",
}

_DEFAULT_IMPORT_PREFIX: tuple[str, ...] = ("Im_",)
_DEFAULT_EXPORT_PREFIX: tuple[str, ...] = ("Ex_", "Si_")
_DEFAULT_BIOMASS_PREFIX: tuple[str, ...] = ("Bio_", "BIOMASS")
_DEFAULT_EXCHANGE_KEY = "u"


# ==========================================================================
# RAILWAY / STATION CONSTANTS
# ==========================================================================

# Standard 20 amino-acid 3-letter codes (lowercase). Used by _met_class.
_AA_CODES: frozenset[str] = frozenset({
    "ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile",
    "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val",
})
_SUGAR_KW: tuple[str, ...] = (
    "glc", "fru", "gal", "man", "lac", "suc", "rib", "xyl", "ara", "tre",
    "malt", "starch", "sorb",
)
_INORG_KW: tuple[str, ...] = (
    "pi", "ppi", "o2", "co2", "h2o", "nh3", "nh4", "no3", "no2",
    "so3", "so4", "hco3", "hco2",
)
_COFAC_KW: tuple[str, ...] = (
    "atp", "adp", "amp", "gtp", "gdp", "ctp", "utp", "ttp",
    "nad", "nadp", "fad", "fmn", "coa", "ubq", "q8", "q10",
)
_LINE_KLASSES: tuple[str, ...] = (
    "amino_acid", "sugar", "cofactor", "inorganic", "other",
)
_LINE_KLASS_COLORS: dict[str, str] = {
    "amino_acid": "#22c55e",
    "sugar":      "#f59e0b",
    "cofactor":   "#a855f7",
    "inorganic":  "#06b6d4",
    "other":      "#B9B84E",
}
_LINE_KLASS_INACTIVE = "#6e7681"
MAX_LINES_PER_PAIR = 10
MAX_RXNS_PER_LINE = 10

_RAIL_ROUTE_PALETTE: tuple[str, ...] = (
    "#f97316", "#22c55e", "#a855f7", "#06b6d4", "#f43f5e",
    "#eab308", "#3b82f6", "#84cc16", "#ec4899", "#14b8a6",
    "#8b5cf6", "#fb923c",
)


# ==========================================================================
# CONFIG LOADING
# ==========================================================================

def _load_toml(path: Path) -> dict:
    if tomllib is None:
        raise ImportError(
            "TOML config requires Python >= 3.11 or 'pip install tomli'."
        )
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def _parse_config(cfg: dict) -> dict:
    """Parse a TOML config dict into internal config format."""
    compartments: dict[str, dict] = {}
    comp_raw = cfg.get("compartment", {})
    if comp_raw:
        for key, entry in comp_raw.items():
            x0, x1, y0, y1 = [float(v) for v in entry["region"]]
            compartments[key] = dict(
                region=(x0, x1, y0, y1),
                color=entry.get("color", "#888888"),
                fill=entry.get("fill", "rgba(128,128,128,0.14)"),
                label=entry.get("label", key),
            )
    else:
        compartments = dict(_DEFAULT_COMPARTMENTS)

    lay = cfg.get("layout", {})
    exchange_key = lay.get("exchange_key", _DEFAULT_EXCHANGE_KEY)
    currency_kw = set(lay.get("currency_keywords", list(_DEFAULT_CURRENCY_KW)))
    import_prefix = tuple(lay.get("import_prefixes", list(_DEFAULT_IMPORT_PREFIX)))
    export_prefix = tuple(lay.get("export_prefixes", list(_DEFAULT_EXPORT_PREFIX)))
    biomass_prefix = tuple(lay.get("biomass_prefixes", list(_DEFAULT_BIOMASS_PREFIX)))

    _alts = "|".join(
        re.escape(k) for k in sorted(compartments.keys(), key=len, reverse=True)
    )
    suffix_re = re.compile(rf"_({_alts})(?:_[a-z])?$", re.IGNORECASE)
    bracket_re = re.compile(rf"\[({_alts})\]$", re.IGNORECASE)
    return dict(
        compartments=compartments,
        exchange_key=exchange_key,
        currency_kw=currency_kw,
        import_prefix=import_prefix,
        export_prefix=export_prefix,
        biomass_prefix=biomass_prefix,
        suffix_re=suffix_re,
        bracket_re=bracket_re,
    )
