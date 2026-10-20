# -*- coding: utf-8 -*-
"""Atmospheric thermodynamics engine for PySkewT v2.

High-precision, vectorized formulations based on standard meteorological literature
(Bolton 1980, WMO standards, AMS Glossary).
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple, Dict, Any, Optional

from .constants import (
    RD, RV, EPSILON, CP_D, KAPPA, LV0, T_ZERO, P_STD, GRAVITY
)

try:
    _trapezoid = np.trapezoid
except AttributeError:
    _trapezoid = getattr(np, "trapz", None)
    if _trapezoid is None:
        from scipy.integrate import trapezoid as _trapezoid


def celsius_to_kelvin(t_c: np.ndarray | float) -> np.ndarray | float:
    """Convert temperature from Celsius to Kelvin."""
    return np.asarray(t_c) + T_ZERO


def kelvin_to_celsius(t_k: np.ndarray | float) -> np.ndarray | float:
    """Convert temperature from Kelvin to Celsius."""
    return np.asarray(t_k) - T_ZERO


def latent_heat_vaporization(t_c: np.ndarray | float) -> np.ndarray | float:
    """Latent heat of vaporization Lv (J / kg) as a function of temperature in Celsius.
    
    Kirchhoff's equation: Lv(T) = Lv0 - 2311 * T_c
    """
    t_c = np.asarray(t_c)
    return LV0 - 2311.0 * t_c


def saturation_vapor_pressure(t_c: np.ndarray | float) -> np.ndarray | float:
    """Saturation vapor pressure over liquid water es(T) in hPa.
    
    Uses Bolton (1980) formula, accurate within 0.1% over [-35 deg C, +35 deg C].
    es(T) = 6.112 * exp(17.67 * T / (T + 243.5))
    """
    t_c = np.asarray(t_c, dtype=float)
    return 6.112 * np.exp((17.67 * t_c) / (t_c + 243.5))


def dewpoint_from_vapor_pressure(e: np.ndarray | float) -> np.ndarray | float:
    """Dew point temperature in Celsius from water vapor pressure e in hPa.
    
    Inverse of the Bolton (1980) formula.
    """
    e = np.asarray(e, dtype=float)
    # Avoid zero or negative values
    e_safe = np.maximum(e, 1e-4)
    ln_val = np.log(e_safe / 6.112)
    return (243.5 * ln_val) / (17.67 - ln_val)


def saturation_mixing_ratio(p: np.ndarray | float, t_c: np.ndarray | float) -> np.ndarray | float:
    """Saturation mixing ratio ws in kg/kg given pressure p in hPa and temperature in deg C.
    
    ws = epsilon * es(T) / (p - es(T))
    """
    p = np.asarray(p, dtype=float)
    es = saturation_vapor_pressure(t_c)
    # Ensure denominator does not blow up near p ~ es
    denom = np.maximum(p - es, 1e-3)
    return (EPSILON * es) / denom


def mixing_ratio_from_dewpoint(p: np.ndarray | float, td_c: np.ndarray | float) -> np.ndarray | float:
    """Actual mixing ratio w in kg/kg given pressure in hPa and dew point in deg C."""
    return saturation_mixing_ratio(p, td_c)


def vapor_pressure_from_mixing_ratio(p: np.ndarray | float, w: np.ndarray | float) -> np.ndarray | float:
    """Water vapor pressure e in hPa from pressure in hPa and mixing ratio w in kg/kg.
    
    e = p * w / (epsilon + w)
    """
    p = np.asarray(p, dtype=float)
    w = np.asarray(w, dtype=float)
    return (p * w) / (EPSILON + w)


def potential_temperature(p: np.ndarray | float, t_c: np.ndarray | float) -> np.ndarray | float:
    """Potential temperature theta in Kelvin from pressure in hPa and temperature in deg C.
    
    theta = (T + 273.15) * (1000 / p) ** kappa
    """
    p = np.asarray(p, dtype=float)
    tk = celsius_to_kelvin(t_c)
    return tk * (P_STD / p) ** KAPPA


def temperature_at_potential_temperature(theta_k: np.ndarray | float, p: np.ndarray | float) -> np.ndarray | float:
    """Temperature in Celsius along a dry adiabat of potential temperature theta_k (Kelvin) at pressure p (hPa).
    
    T = theta * (p / 1000) ** kappa - 273.15
    """
    p = np.asarray(p, dtype=float)
    theta_k = np.asarray(theta_k, dtype=float)
    return theta_k * (p / P_STD) ** KAPPA - T_ZERO


def virtual_temperature(t_c: np.ndarray | float, w: np.ndarray | float) -> np.ndarray | float:
    """Virtual temperature Tv in Celsius from temperature in deg C and mixing ratio in kg/kg.
    
    Tv = (T + 273.15) * (1 + 0.608 * w) - 273.15
    """
    tk = celsius_to_kelvin(t_c)
    w = np.asarray(w, dtype=float)
    return tk * (1.0 + 0.608 * w) - T_ZERO


def equivalent_potential_temperature(
    p: np.ndarray | float, t_c: np.ndarray | float, td_c: np.ndarray | float
) -> np.ndarray | float:
    """Equivalent potential temperature theta_e in Kelvin (Bolton 1980 formulation)."""
    p = np.asarray(p, dtype=float)
    t_c = np.asarray(t_c, dtype=float)
    td_c = np.asarray(td_c, dtype=float)

    tk = celsius_to_kelvin(t_c)
    w = mixing_ratio_from_dewpoint(p, td_c)  # kg/kg
    e = vapor_pressure_from_mixing_ratio(p, w)

    # Temperature at LCL (Bolton 1980 eq 22)
    t_lcl = 1.0 / (1.0 / (tk - 55.0) - np.log(np.maximum(e / saturation_vapor_pressure(t_c), 1e-4)) / 2840.0) + 55.0

    # Dry potential temperature at LCL (Bolton eq 24)
    theta_dl = tk * (P_STD / (p - e)) ** KAPPA * (tk / t_lcl) ** (0.28 * w)

    # Equivalent potential temperature (Bolton eq 38)
    # w in g/kg for the empirical factor
    r = w * 1000.0
    return theta_dl * np.exp((3036.0 / t_lcl - 1.78) * w * (1.0 + 0.448 * w))


def moist_lapse_rate(p: float, t_c: float) -> float:
    """Moist adiabatic lapse rate dT/dp (deg C / hPa).
    
    dT/dp = (Rd * T / (p * cp)) * [1 + (Lv * ws) / (Rd * T)] / [1 + (epsilon * Lv^2 * ws) / (cp * Rd * T^2)]
    """
    tk = t_c + T_ZERO
    if tk <= 50.0:
        return 0.0

    lv = latent_heat_vaporization(t_c)
    es = float(saturation_vapor_pressure(t_c))
    denom_p = max(p - es, 0.1)
    ws = (EPSILON * es) / denom_p

    numerator = (RD * tk / (p * CP_D)) * (1.0 + (lv * ws) / (RD * tk))
    denominator = 1.0 + (EPSILON * (lv ** 2) * ws) / (CP_D * RD * (tk ** 2))
    return float(numerator / denominator)


def moist_adiabat_profile(
    theta_w_c: float, p_levels: np.ndarray, p_ref: float = P_STD
) -> np.ndarray:
    """Calculate temperature in deg C along a moist pseudo-adiabat defined by wet-bulb temperature theta_w_c at p_ref.
    
    Solves dT/dp using 4th-order Runge-Kutta across p_levels.
    """
    p_levels = np.asarray(p_levels, dtype=float)
    sort_idx = np.argsort(p_levels)[::-1]
    p_arr = p_levels[sort_idx]
    
    idx_ref = int(np.argmin(np.abs(p_arr - p_ref)))
    t_arr = np.zeros_like(p_arr)
    t_arr[idx_ref] = theta_w_c

    # Integrate downward (higher p: p_ref -> max(p))
    for i in range(idx_ref, 0, -1):
        dp = p_arr[i - 1] - p_arr[i]
        p_c = p_arr[i]
        t_c = t_arr[i]
        k1 = moist_lapse_rate(p_c, t_c)
        k2 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k1)
        k3 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k2)
        k4 = moist_lapse_rate(p_c + dp, t_c + dp * k3)
        t_arr[i - 1] = t_c + (dp / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    # Integrate upward (lower p: p_ref -> min(p))
    for i in range(idx_ref, len(p_arr) - 1):
        dp = p_arr[i + 1] - p_arr[i]
        p_c = p_arr[i]
        t_c = t_arr[i]
        k1 = moist_lapse_rate(p_c, t_c)
        k2 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k1)
        k3 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k2)
        k4 = moist_lapse_rate(p_c + dp, t_c + dp * k3)
        t_arr[i + 1] = t_c + (dp / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    t_out = np.zeros_like(t_arr)
    t_out[sort_idx] = t_arr
    return t_out


def lifting_condensation_level(
    p_sfc: float, t_sfc: float, td_sfc: float
) -> Tuple[float, float]:
    """Calculate Lifting Condensation Level (LCL) pressure (hPa) and temperature (deg C).
    
    Uses exact Bolton (1980) formulas.
    """
    tk = t_sfc + T_ZERO
    tdk = td_sfc + T_ZERO
    
    # Relative humidity (fraction)
    e = float(saturation_vapor_pressure(td_sfc))
    es = float(saturation_vapor_pressure(t_sfc))
    rh = max(min(e / es, 1.0), 1e-4)

    # Temperature at LCL (Kelvin)
    t_lcl_k = 1.0 / (1.0 / (tk - 55.0) - np.log(rh) / 2840.0) + 55.0
    t_lcl_c = t_lcl_k - T_ZERO

    # Pressure at LCL (hPa) following Poisson's equation
    p_lcl = p_sfc * (t_lcl_k / tk) ** (1.0 / KAPPA)

    return float(p_lcl), float(t_lcl_c)


def parcel_ascent_profile(
    p_levels: np.ndarray, p_sfc: float, t_sfc: float, td_sfc: float
) -> np.ndarray:
    """Calculate the temperature profile of a surface air parcel lifted through p_levels.
    
    Follows dry adiabat from surface to LCL, and moist adiabat above LCL.
    """
    p_levels = np.asarray(p_levels, dtype=float)
    p_lcl, t_lcl = lifting_condensation_level(p_sfc, t_sfc, td_sfc)
    theta_sfc = (t_sfc + T_ZERO) * (P_STD / p_sfc) ** KAPPA

    t_parcel = np.zeros_like(p_levels, dtype=float)

    # Levels below LCL: dry adiabatic ascent
    mask_dry = p_levels >= p_lcl
    if np.any(mask_dry):
        t_parcel[mask_dry] = temperature_at_potential_temperature(theta_sfc, p_levels[mask_dry])

    # Levels above LCL: moist adiabatic ascent
    mask_moist = p_levels < p_lcl
    if np.any(mask_moist):
        p_moist = np.sort(p_levels[mask_moist])[::-1]
        p_dense = np.concatenate(([p_lcl], p_moist))
        t_dense = np.zeros_like(p_dense)
        t_dense[0] = t_lcl
        for i in range(len(p_dense) - 1):
            dp = p_dense[i + 1] - p_dense[i]
            p_c = p_dense[i]
            t_c = t_dense[i]
            k1 = moist_lapse_rate(p_c, t_c)
            k2 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k1)
            k3 = moist_lapse_rate(p_c + 0.5 * dp, t_c + 0.5 * dp * k2)
            k4 = moist_lapse_rate(p_c + dp, t_c + dp * k3)
            t_dense[i + 1] = t_c + (dp / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        t_parcel[mask_moist] = np.interp(p_levels[mask_moist], p_dense[::-1], t_dense[::-1])

    return t_parcel


def calculate_convective_indices(
    p: np.ndarray, temp: np.ndarray, dwpt: np.ndarray
) -> Dict[str, Any]:
    """Calculate thermodynamic and convective stability indices:
    LCL, LFC, EL, CAPE, CIN, Lifted Index (LI), K-Index, Total Totals, Precipitable Water.
    """
    # Sort from surface to top (p decreasing)
    idx_sort = np.argsort(p)[::-1]
    p_sorted = p[idx_sort]
    t_sorted = temp[idx_sort]
    td_sorted = dwpt[idx_sort]

    p_sfc = p_sorted[0]
    t_sfc = t_sorted[0]
    td_sfc = td_sorted[0]

    # LCL
    p_lcl, t_lcl = lifting_condensation_level(p_sfc, t_sfc, td_sfc)

    # Dense pressure grid for high-precision integration of CAPE / CIN
    p_grid = np.linspace(p_sfc, np.min(p_sorted), 400)
    t_env = np.interp(p_grid, p_sorted[::-1], t_sorted[::-1])
    td_env = np.interp(p_grid, p_sorted[::-1], td_sorted[::-1])
    t_parcel = parcel_ascent_profile(p_grid, p_sfc, t_sfc, td_sfc)

    # Virtual temperatures for accurate buoyancy
    w_env = mixing_ratio_from_dewpoint(p_grid, td_env)
    # Parcel mixing ratio is constant below LCL, then saturated above
    w_sfc = mixing_ratio_from_dewpoint(p_sfc, td_sfc)
    w_parcel = np.where(p_grid >= p_lcl, w_sfc, saturation_mixing_ratio(p_grid, t_parcel))

    tv_env = virtual_temperature(t_env, w_env)
    tv_parcel = virtual_temperature(t_parcel, w_parcel)

    # Buoyancy: B = g * (Tv_parcel - Tv_env) / Tv_env (or Rd * (Tv_parcel - Tv_env) dln(p))
    delta_tv = tv_parcel - tv_env

    # Find LFC: first crossing from negative to positive buoyancy above LCL
    mask_above_lcl = p_grid <= p_lcl
    p_above = p_grid[mask_above_lcl]
    delta_tv_above = delta_tv[mask_above_lcl]

    p_lfc, t_lfc = None, None
    p_el, t_el = None, None
    cape = 0.0
    cin = 0.0

    pos_buoyant = delta_tv_above > 0
    if np.any(pos_buoyant):
        first_pos_idx = np.where(pos_buoyant)[0][0]
        p_lfc = float(p_above[first_pos_idx])
        t_lfc = float(t_parcel[mask_above_lcl][first_pos_idx])

        # Equilibrium Level (EL): last crossing where buoyancy turns negative
        last_pos_idx = np.where(pos_buoyant)[0][-1]
        p_el = float(p_above[last_pos_idx])
        t_el = float(t_parcel[mask_above_lcl][last_pos_idx])

        # Integrate CAPE between LFC and EL
        mask_cape = (p_grid <= p_lfc) & (p_grid >= p_el) & (delta_tv > 0)
        if np.any(mask_cape):
            p_c = p_grid[mask_cape]
            dtv_c = delta_tv[mask_cape]
            # Standard Rd * delta_tv * -dln(p)
            cape = float(-_trapezoid(RD * dtv_c, np.log(p_c)))
            cape = max(0.0, cape)

        # Integrate CIN between surface and LFC
        mask_cin = (p_grid >= p_lfc) & (delta_tv < 0)
        if np.any(mask_cin):
            p_cn = p_grid[mask_cin]
            dtv_cn = delta_tv[mask_cin]
            cin = float(_trapezoid(RD * dtv_cn, np.log(p_cn)))
            cin = min(0.0, cin)

    # Lifted Index (LI) at 500 hPa
    t_env_500 = float(np.interp(500.0, p_sorted[::-1], t_sorted[::-1]))
    t_parcel_500 = float(np.interp(500.0, p_grid[::-1], t_parcel[::-1]))
    lifted_index = round(t_env_500 - t_parcel_500, 2)

    # K-Index = (T850 - T500) + Td850 - (T700 - Td700)
    t_850 = float(np.interp(850.0, p_sorted[::-1], t_sorted[::-1]))
    td_850 = float(np.interp(850.0, p_sorted[::-1], td_sorted[::-1]))
    t_700 = float(np.interp(700.0, p_sorted[::-1], t_sorted[::-1]))
    td_700 = float(np.interp(700.0, p_sorted[::-1], td_sorted[::-1]))
    k_index = round((t_850 - t_env_500) + td_850 - (t_700 - td_700), 1)

    # Total Totals (TT) = (T850 - T500) + (Td850 - T500)
    total_totals = round((t_850 - t_env_500) + (td_850 - t_env_500), 1)

    # Precipitable Water (PW) in mm: (1 / (g * rho_w)) * int(w dp)
    # dp in Pa = 100 * dp_hPa
    w_profile = mixing_ratio_from_dewpoint(p_sorted, td_sorted)  # kg/kg
    pw = float((100.0 / GRAVITY) * _trapezoid(w_profile, -p_sorted))  # mm (kg/m^2)

    return {
        "p_sfc": round(p_sfc, 1),
        "t_sfc": round(t_sfc, 1),
        "td_sfc": round(td_sfc, 1),
        "p_lcl": round(p_lcl, 1),
        "t_lcl": round(t_lcl, 1),
        "p_lfc": round(p_lfc, 1) if p_lfc is not None else None,
        "t_lfc": round(t_lfc, 1) if t_lfc is not None else None,
        "p_el": round(p_el, 1) if p_el is not None else None,
        "t_el": round(t_el, 1) if t_el is not None else None,
        "cape": round(cape, 1),
        "cin": round(cin, 1),
        "lifted_index": lifted_index,
        "k_index": k_index,
        "total_totals": total_totals,
        "precipitable_water": round(pw, 1),
        "parcel_p": p_grid,
        "parcel_t": t_parcel,
        "env_t": t_env,
        "env_td": td_env,
    }
