# -*- coding: utf-8 -*-
"""Color themes and visual styling palettes for PySkewT v2."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Union, List


@dataclass
class Theme:
    """Encapsulates all visual colors and styling for a Skew-T diagram."""
    name: str
    description: str
    # Canvas and typography
    bg_color: str
    fg_color: str
    # Thermodynamic background lines
    grid_isobars: str
    grid_isotherms: str
    grid_zero_isotherm: str
    grid_dry_adiabats: str
    grid_moist_adiabats: str
    grid_mixing_ratio: str
    # Sounding profiles
    temp_profile: str
    dewpoint_profile: str
    parcel_profile: str
    # Convective shading
    cape_fill: str
    cin_fill: str
    lcl_marker: str
    # Wind barbs
    wind_barbs: str
    wind_flags: str
    # Diagnostic HUD box
    hud_bg: str
    hud_border: str
    hud_text: str


# -----------------------------------------------------------------------------
# Curated Palettes
# -----------------------------------------------------------------------------

THEME_MODERN_LIGHT = Theme(
    name="modern_light",
    description="Paleta moderna de alto contraste y claridad científica (predeterminada)",
    bg_color="#ffffff",
    fg_color="#1e293b",
    # Grid lines: subtle slate tones so sounding curves pop
    grid_isobars="#94a3b8",
    grid_isotherms="#cbd5e1",
    grid_zero_isotherm="#0284c7",       # Vivid azure
    grid_dry_adiabats="#ea580c",        # Rich copper/terracotta
    grid_moist_adiabats="#059669",      # Forest emerald
    grid_mixing_ratio="#0891b2",        # Cyan teal
    # Profiles
    temp_profile="#e11d48",             # Crisp vivid crimson
    dewpoint_profile="#15803d",         # Deep emerald green
    parcel_profile="#374151",           # Charcoal slate
    # Convective
    cape_fill="#f43f5e",                # Vibrant coral rose
    cin_fill="#38bdf8",                 # Sky blue
    lcl_marker="#2563eb",               # Royal blue
    # Winds
    wind_barbs="#1e293b",
    wind_flags="#e11d48",
    # HUD
    hud_bg="#ffffff",
    hud_border="#cbd5e1",
    hud_text="#0f172a",
)

THEME_DARK_MODE = Theme(
    name="dark_mode",
    description="Modo oscuro de alta luminiscencia para monitores y guardias operacionales",
    bg_color="#0f172a",                 # Deep dark navy/slate
    fg_color="#f8fafc",                 # Clean light text
    # Grid lines
    grid_isobars="#334155",
    grid_isotherms="#475569",
    grid_zero_isotherm="#38bdf8",       # Glowing cyan
    grid_dry_adiabats="#f59e0b",        # Warm amber
    grid_moist_adiabats="#34d399",      # Mint green
    grid_mixing_ratio="#2dd4bf",        # Turquoise
    # Profiles
    temp_profile="#fb7185",             # Neon coral
    dewpoint_profile="#4ade80",         # Vibrant lime
    parcel_profile="#e2e8f0",           # Off-white
    # Convective
    cape_fill="#f43f5e",                # Translucent neon rose
    cin_fill="#06b6d4",                 # Translucent cyan
    lcl_marker="#60a5fa",               # Electric light blue
    # Winds
    wind_barbs="#f1f5f9",
    wind_flags="#fb7185",
    # HUD
    hud_bg="#1e293b",
    hud_border="#475569",
    hud_text="#f8fafc",
)

THEME_CLASSIC_METEO = Theme(
    name="classic_meteo",
    description="Estilo meteorológico tradicional de papel y publicaciones de la OMM",
    bg_color="#faf8f0",                 # Warm paper / parchment
    fg_color="#292524",                 # Warm dark brown
    # Grid lines
    grid_isobars="#a89f91",
    grid_isotherms="#cbbfae",
    grid_zero_isotherm="#1e3a8a",       # Deep navy
    grid_dry_adiabats="#c2410c",        # Burnt sienna
    grid_moist_adiabats="#166534",      # Traditional olive forest
    grid_mixing_ratio="#0369a1",        # Slate blue
    # Profiles
    temp_profile="#b91c1c",             # Classic sounding red
    dewpoint_profile="#15803d",         # Classic sounding green
    parcel_profile="#44403c",           # Dark taupe
    # Convective
    cape_fill="#ef4444",
    cin_fill="#7dd3fc",
    lcl_marker="#1d4ed8",
    # Winds
    wind_barbs="#292524",
    wind_flags="#b91c1c",
    # HUD
    hud_bg="#fefcf8",
    hud_border="#d6cfc7",
    hud_text="#1c1917",
)

THEME_ACCESSIBLE = Theme(
    name="accessible",
    description="Paleta accesible (CUD / Okabe-Ito) segura para daltonismo (protanopía y deuteranopía)",
    bg_color="#ffffff",
    fg_color="#18181b",
    # Grid lines
    grid_isobars="#a1a1aa",
    grid_isotherms="#d4d4d8",
    grid_zero_isotherm="#0284c7",
    grid_dry_adiabats="#d97706",
    grid_moist_adiabats="#0d9488",
    grid_mixing_ratio="#64748b",
    # Profiles: Vermilion vs Cobalt (No red-green confusion)
    temp_profile="#d55e00",             # Vermilion / Reddish-orange
    dewpoint_profile="#0072b2",         # Cobalt / Deep Sky Blue
    parcel_profile="#27272a",
    # Convective
    cape_fill="#cc79a7",                # Reddish purple
    cin_fill="#56b4e9",                 # Soft sky blue
    lcl_marker="#0072b2",
    # Winds
    wind_barbs="#18181b",
    wind_flags="#d55e00",
    # HUD
    hud_bg="#ffffff",
    hud_border="#a1a1aa",
    hud_text="#09090b",
)


THEMES: Dict[str, Theme] = {
    "modern": THEME_MODERN_LIGHT,
    "modern_light": THEME_MODERN_LIGHT,
    "default": THEME_MODERN_LIGHT,
    "dark": THEME_DARK_MODE,
    "dark_mode": THEME_DARK_MODE,
    "night": THEME_DARK_MODE,
    "classic": THEME_CLASSIC_METEO,
    "classic_meteo": THEME_CLASSIC_METEO,
    "retro": THEME_CLASSIC_METEO,
    "accessible": THEME_ACCESSIBLE,
    "colorblind": THEME_ACCESSIBLE,
    "cud": THEME_ACCESSIBLE,
}


def get_theme(name_or_theme: Union[str, Theme, None] = None) -> Theme:
    """Retrieve a Theme instance by name or alias, defaulting to modern_light."""
    if isinstance(name_or_theme, Theme):
        return name_or_theme

    if name_or_theme is None:
        return THEME_MODERN_LIGHT

    key = str(name_or_theme).strip().lower().replace("-", "_")
    if key in THEMES:
        return THEMES[key]

    valid = ", ".join(sorted(set(THEMES.keys())))
    raise ValueError(f"Tema desconocido '{name_or_theme}'. Temas disponibles: {valid}")


def list_themes() -> List[str]:
    """Return the canonical list of available theme names."""
    return ["modern_light", "dark_mode", "classic_meteo", "accessible"]
