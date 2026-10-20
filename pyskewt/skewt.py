# -*- coding: utf-8 -*-
"""High-performance, publication-grade Skew-T ln-P thermodynamic diagram for PySkewT v2."""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from typing import Optional, Union, Tuple, Dict, Any, List

from .constants import P_STD, DEFAULT_P_BOT, DEFAULT_P_TOP
from .theme import Theme, get_theme
from .thermo import (
    temperature_at_potential_temperature,
    dewpoint_from_vapor_pressure,
    moist_adiabat_profile,
    calculate_convective_indices,
)
from .sounding import Sounding


class SkewT:
    """A Skew-T ln-P thermodynamic diagram for atmospheric sounding analysis."""

    def __init__(
        self,
        fig: Optional[plt.Figure] = None,
        ax: Optional[plt.Axes] = None,
        ax_wind: Optional[plt.Axes] = None,
        p_bot: float = DEFAULT_P_BOT,
        p_top: float = DEFAULT_P_TOP,
        t_min: float = -45.0,
        t_max: float = 45.0,
        beta: float = 40.0,
        with_wind_barbs: bool = True,
        theme: Union[str, Theme] = "modern",
        figsize: Tuple[float, float] = (10, 11),
    ):
        self.p_bot = float(p_bot)
        self.p_top = float(p_top)
        self.t_min = float(t_min)
        self.t_max = float(t_max)
        self.beta = float(beta)
        self.with_wind_barbs = with_wind_barbs
        self.theme = get_theme(theme)

        if fig is None or ax is None:
            if with_wind_barbs and ax_wind is None:
                self.fig, (self.ax, self.ax_wind) = plt.subplots(
                    1, 2, figsize=figsize, gridspec_kw={"width_ratios": [15, 1]}, sharey=True
                )
                self.fig.subplots_adjust(wspace=0.03, left=0.08, right=0.96, top=0.94, bottom=0.07)
            else:
                self.fig, self.ax = plt.subplots(figsize=figsize)
                self.ax_wind = ax_wind
                self.fig.subplots_adjust(left=0.08, right=0.94, top=0.94, bottom=0.07)
        else:
            self.fig = fig
            self.ax = ax
            self.ax_wind = ax_wind

        self._setup_axes()

    def _setup_axes(self) -> None:
        """Configure coordinate limits, scales, colors, and labels according to theme."""
        # Canvas background colors
        self.fig.patch.set_facecolor(self.theme.bg_color)
        self.ax.set_facecolor(self.theme.bg_color)

        ax = self.ax
        ax.set_yscale("log")
        ax.set_ylim(self.p_bot, self.p_top)
        ax.set_xlim(self.t_min, self.t_max)

        # Spines and axis border colors
        for spine in ax.spines.values():
            spine.set_color(self.theme.fg_color)
            spine.set_linewidth(0.8)

        # Standard pressure levels in hPa
        standard_p = np.array([1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100])
        mask = (standard_p <= self.p_bot) & (standard_p >= self.p_top)
        ticks_p = standard_p[mask]

        ax.set_yticks(ticks_p)
        ax.set_yticklabels(
            [f"{int(p)}" for p in ticks_p],
            fontsize=9, fontweight="bold", color=self.theme.fg_color
        )
        ax.tick_params(
            which="both", left=True, right=True, labelleft=True, labelright=True,
            direction="in", colors=self.theme.fg_color
        )

        # Bottom Temperature axis
        t_ticks = np.arange(int(self.t_min), int(self.t_max) + 1, 10)
        ax.set_xticks(t_ticks)
        ax.set_xticklabels(
            [f"{t}°C" for t in t_ticks],
            fontsize=9, fontweight="bold", color=self.theme.fg_color
        )
        ax.tick_params(bottom=True, top=False, direction="in", colors=self.theme.fg_color)

        ax.set_xlabel("Temperatura (°C) en 1000 hPa", fontsize=10, fontweight="bold", labelpad=8, color=self.theme.fg_color)
        ax.set_ylabel("Presión (hPa)", fontsize=10, fontweight="bold", labelpad=8, color=self.theme.fg_color)
        ax.grid(False)

        # Setup wind axis if enabled
        if self.ax_wind is not None:
            self.ax_wind.set_facecolor(self.theme.bg_color)
            self.ax_wind.set_ylim(self.p_bot, self.p_top)
            self.ax_wind.set_xlim(-1.2, 1.2)
            self.ax_wind.axis("off")

    def t_p_to_x(self, temp_c: np.ndarray | float, p_hpa: np.ndarray | float) -> np.ndarray:
        """Skew-T horizontal transformation: x = T - beta * ln(p / 1000)."""
        temp_c = np.asarray(temp_c, dtype=float)
        p_hpa = np.asarray(p_hpa, dtype=float)
        return temp_c - self.beta * np.log(p_hpa / P_STD)

    def draw_isobars(
        self,
        color: Optional[str] = None,
        linewidth: float = 0.7,
        levels: Optional[List[float]] = None,
    ) -> None:
        """Draw horizontal isobar lines at standard atmospheric pressure levels."""
        c = color or self.theme.grid_isobars
        if levels is None:
            levels = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100]

        lines = []
        for p in levels:
            if self.p_top <= p <= self.p_bot:
                lines.append([(self.t_min - 20, p), (self.t_max + 60, p)])

        lc = LineCollection(lines, colors=c, linewidths=linewidth, zorder=1)
        self.ax.add_collection(lc)

    def draw_isotherms(
        self,
        step: float = 10.0,
        color: Optional[str] = None,
        zero_color: Optional[str] = None,
        linewidth: float = 0.5,
        highlight_zero: bool = True,
    ) -> None:
        """Draw slanted isotherm lines across the diagram using LineCollection."""
        c = color or self.theme.grid_isotherms
        z_c = zero_color or self.theme.grid_zero_isotherm

        p_pts = np.array([self.p_bot, self.p_top])
        t_range = np.arange(-100, 80, step)

        lines = []
        zero_line = None

        for t in t_range:
            x_pts = self.t_p_to_x(t, p_pts)
            segment = np.column_stack([x_pts, p_pts])
            if highlight_zero and abs(t) < 1e-3:
                zero_line = segment
            else:
                lines.append(segment)

        if lines:
            lc = LineCollection(lines, colors=c, linewidths=linewidth, zorder=1)
            self.ax.add_collection(lc)

        # Highlight 0°C isotherm (freezing level)
        if zero_line is not None:
            self.ax.plot(
                zero_line[:, 0], zero_line[:, 1],
                color=z_c, linewidth=1.2, linestyle="-",
                label="Isoterma 0°C", zorder=2
            )

    def draw_dry_adiabats(
        self,
        theta_step: float = 10.0,
        color: Optional[str] = None,
        linewidth: float = 0.55,
        alpha: float = 0.65,
    ) -> None:
        """Draw dry adiabats (isentropes) lines of constant potential temperature theta."""
        c = color or self.theme.grid_dry_adiabats
        p_pts = np.linspace(self.p_bot, self.p_top, 80)
        theta_range = np.arange(220, 480, theta_step)

        lines = []
        for th in theta_range:
            t_c = temperature_at_potential_temperature(th, p_pts)
            x_pts = self.t_p_to_x(t_c, p_pts)
            lines.append(np.column_stack([x_pts, p_pts]))

        lc = LineCollection(lines, colors=c, linewidths=linewidth, alpha=alpha, zorder=1)
        self.ax.add_collection(lc)

    def draw_moist_adiabats(
        self,
        tw_range: Optional[List[float]] = None,
        color: Optional[str] = None,
        linewidth: float = 0.65,
        alpha: float = 0.7,
    ) -> None:
        """Draw moist pseudo-adiabats using high-order Runge-Kutta integration."""
        c = color or self.theme.grid_moist_adiabats
        if tw_range is None:
            tw_range = [-8, -4, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36]

        p_pts = np.linspace(self.p_bot, self.p_top, 90)
        lines = []
        for tw in tw_range:
            t_prof = moist_adiabat_profile(tw, p_pts)
            x_pts = self.t_p_to_x(t_prof, p_pts)
            lines.append(np.column_stack([x_pts, p_pts]))

        lc = LineCollection(lines, colors=c, linewidths=linewidth, alpha=alpha, zorder=1)
        self.ax.add_collection(lc)

    def draw_mixing_ratio_lines(
        self,
        ws_values: Optional[List[float]] = None,
        color: Optional[str] = None,
        linewidth: float = 0.55,
        alpha: float = 0.65,
    ) -> None:
        """Draw saturation mixing ratio isohumes (ws in g/kg)."""
        c = color or self.theme.grid_mixing_ratio
        if ws_values is None:
            ws_values = [0.4, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 28.0]

        # Only draw down to 400 hPa where moisture is significant
        p_pts = np.linspace(self.p_bot, 400.0, 60)
        lines = []
        for w_g in ws_values:
            w_kg = w_g * 1e-3
            e = (p_pts * w_kg) / (0.62198 + w_kg)
            t_c = dewpoint_from_vapor_pressure(e)
            x_pts = self.t_p_to_x(t_c, p_pts)
            lines.append(np.column_stack([x_pts, p_pts]))

            # Add text label near bottom
            idx_lbl = 5
            self.ax.text(
                x_pts[idx_lbl], p_pts[idx_lbl],
                f" {w_g:g}",
                color=c, fontsize=6.5, rotation=45,
                verticalalignment="bottom", alpha=0.85
            )

        lc = LineCollection(
            lines, colors=c, linewidths=linewidth, linestyles="dashed", alpha=alpha, zorder=1
        )
        self.ax.add_collection(lc)

    def draw_background(self) -> None:
        """Draw complete thermodynamic background grid (isobars, isotherms, adiabats, isohumes)."""
        self.draw_isobars()
        self.draw_isotherms()
        self.draw_dry_adiabats()
        self.draw_moist_adiabats()
        self.draw_mixing_ratio_lines()

    def plot_sounding(
        self,
        sounding: Union[Sounding, str],
        plot_parcel: bool = True,
        shade_cape_cin: bool = True,
        plot_winds: bool = True,
        title: Optional[str] = "Diagrama Skew-T ln-P",
        color_t: Optional[str] = None,
        color_td: Optional[str] = None,
        color_parcel: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Plot complete sounding profile, convective parcel ascent, CAPE/CIN shading and winds."""
        if isinstance(sounding, str):
            sounding = Sounding.from_file(sounding)

        c_t = color_t or self.theme.temp_profile
        c_td = color_td or self.theme.dewpoint_profile
        c_parcel = color_parcel or self.theme.parcel_profile

        p = sounding.pressure
        t = sounding.temperature
        td = sounding.dewpoint

        # Transform to skewed coordinate x
        x_t = self.t_p_to_x(t, p)
        x_td = self.t_p_to_x(td, p)

        # Draw temperature profile
        self.ax.plot(
            x_t, p, color=c_t, linewidth=2.0, label="Temperatura (T)", zorder=4
        )

        # Draw dew point profile
        self.ax.plot(
            x_td, p, color=c_td, linewidth=2.0, label="Punto de Rocío (Td)", zorder=4
        )

        # Compute thermodynamic indices
        indices = calculate_convective_indices(p, t, td)

        if plot_parcel and indices["p_lcl"] is not None:
            p_parcel = indices["parcel_p"]
            t_parcel = indices["parcel_t"]
            x_parcel = self.t_p_to_x(t_parcel, p_parcel)

            self.ax.plot(
                x_parcel, p_parcel, color=c_parcel, linewidth=1.5, linestyle="--",
                label="Ascenso de Parcela", zorder=5
            )

            # Mark LCL level
            p_lcl = indices["p_lcl"]
            t_lcl = indices["t_lcl"]
            x_lcl = self.t_p_to_x(t_lcl, p_lcl)
            self.ax.plot(
                x_lcl, p_lcl, marker="o", markersize=5, color=self.theme.lcl_marker, zorder=6
            )
            self.ax.text(
                x_lcl + 1.5, p_lcl, f"LCL: {p_lcl:.0f} hPa",
                fontsize=7.5, fontweight="bold", color=self.theme.lcl_marker, verticalalignment="center"
            )

            # Shade CAPE and CIN
            if shade_cape_cin and indices["p_lfc"] is not None:
                p_lfc = indices["p_lfc"]
                p_el = indices["p_el"]
                t_env = indices["env_t"]
                x_env = self.t_p_to_x(t_env, p_parcel)

                # Shaded CAPE
                mask_cape = (p_parcel <= p_lfc) & (p_parcel >= p_el) & (x_parcel > x_env)
                if np.any(mask_cape):
                    self.ax.fill_betweenx(
                        p_parcel, x_env, x_parcel,
                        where=mask_cape, color=self.theme.cape_fill, alpha=0.35,
                        label=f"CAPE ({indices['cape']:.0f} J/kg)", zorder=3
                    )

                # Shaded CIN
                mask_cin = (p_parcel >= p_lfc) & (x_parcel < x_env)
                if np.any(mask_cin):
                    self.ax.fill_betweenx(
                        p_parcel, x_env, x_parcel,
                        where=mask_cin, color=self.theme.cin_fill, alpha=0.35,
                        label=f"CIN ({indices['cin']:.0f} J/kg)", zorder=3
                    )

        # Plot wind barbs if present and enabled
        if plot_winds and sounding.u is not None and self.ax_wind is not None:
            self._plot_wind_barbs(sounding.pressure, sounding.u, sounding.v)

        # Draw indices HUD
        self.draw_indices_hud(indices)

        if title:
            self.ax.set_title(title, fontsize=12, fontweight="bold", pad=12, color=self.theme.fg_color)

        leg = self.ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
        leg.get_frame().set_facecolor(self.theme.hud_bg)
        leg.get_frame().set_edgecolor(self.theme.hud_border)
        for text in leg.get_texts():
            text.set_color(self.theme.hud_text)

        return indices

    def _plot_wind_barbs(self, p: np.ndarray, u: np.ndarray, v: np.ndarray) -> None:
        """Plot wind barbs in the dedicated vertical axis with intelligent subsampling."""
        valid = np.isfinite(p) & np.isfinite(u) & np.isfinite(v) & (p >= self.p_top) & (p <= self.p_bot)
        p_v = p[valid]
        u_v = u[valid]
        v_v = v[valid]

        if len(p_v) == 0:
            return

        target_p = np.arange(self.p_bot, self.p_top - 1, -40)
        idx_chosen = []
        for tp in target_p:
            dist = np.abs(p_v - tp)
            if np.min(dist) < 25:
                idx_chosen.append(np.argmin(dist))

        idx_chosen = np.unique(idx_chosen)
        p_sub = p_v[idx_chosen]
        u_sub = u_v[idx_chosen]
        v_sub = v_v[idx_chosen]

        self.ax_wind.barbs(
            np.zeros_like(p_sub), p_sub, u_sub, v_sub,
            length=6, barbcolor=self.theme.wind_barbs, flagcolor=self.theme.wind_flags, linewidth=0.8,
            sizes=dict(emptybarb=0.08, spacing=0.2, height=0.5)
        )

    def draw_indices_hud(self, indices: Dict[str, Any]) -> None:
        """Render a clean information HUD box with key convective stability indices."""
        cape_str = f"{indices['cape']:.0f} J/kg"
        cin_str = f"{indices['cin']:.0f} J/kg"
        lcl_str = f"{indices['p_lcl']:.0f} hPa ({indices['t_lcl']:.1f}°C)"
        lfc_str = f"{indices['p_lfc']:.0f} hPa" if indices['p_lfc'] else "N/A"
        el_str = f"{indices['p_el']:.0f} hPa" if indices['p_el'] else "N/A"
        li_str = f"{indices['lifted_index']:.1f}°C"
        k_str = f"{indices['k_index']:.1f}°C"
        pw_str = f"{indices['precipitable_water']:.1f} mm"

        text_content = (
            f"Parámetros Convectivos:\n"
            f"• CAPE: {cape_str}\n"
            f"• CIN:  {cin_str}\n"
            f"• LCL:  {lcl_str}\n"
            f"• LFC:  {lfc_str}\n"
            f"• EL:   {el_str}\n"
            f"• Lifted Index (LI): {li_str}\n"
            f"• K-Index: {k_str}\n"
            f"• Agua Precipitable: {pw_str}"
        )

        props = dict(
            boxstyle="round,pad=0.5",
            facecolor=self.theme.hud_bg,
            edgecolor=self.theme.hud_border,
            alpha=0.92
        )
        self.ax.text(
            0.02, 0.98, text_content,
            transform=self.ax.transAxes,
            fontsize=8, verticalalignment="top", fontfamily="monospace",
            color=self.theme.hud_text,
            bbox=props, zorder=10
        )

    def save(self, filepath: str, dpi: int = 300) -> None:
        """Save high-resolution diagram to disk preserving theme background."""
        self.fig.savefig(
            filepath,
            dpi=dpi,
            bbox_inches="tight",
            facecolor=self.fig.get_facecolor(),
            edgecolor="none"
        )

    def show(self) -> None:
        """Display the diagram."""
        plt.show()

    # =========================================================================
    # Legacy v1 Backwards-Compatibility Layer
    # =========================================================================
    def T(self) -> None:
        """Legacy v1 method: draw isotherms."""
        self.draw_isotherms()

    def ws(self, w_s: Optional[np.ndarray] = None) -> None:
        """Legacy v1 method: draw mixing ratios."""
        self.draw_mixing_ratio_lines()

    def TP(self, Tpot: Optional[np.ndarray] = None, p0: float = 1e5) -> None:
        """Legacy v1 method: draw dry adiabats."""
        self.draw_dry_adiabats()

    def tpe(self, unT: Optional[np.ndarray] = None, unp: float = 1e5) -> None:
        """Legacy v1 method: draw moist adiabats."""
        self.draw_moist_adiabats()

    def graficar(self, mat: np.ndarray) -> None:
        """Legacy v1 method: plot raw data matrix."""
        p = mat[:, 0]
        if np.max(p) > 2000:
            p = p / 100.0
        temp = mat[:, 2]
        dwpt = mat[:, 3]
        sounding = Sounding.from_arrays(pressure=p, temperature=temp, dewpoint=dwpt)
        self.plot_sounding(sounding, plot_parcel=True, shade_cape_cin=True, plot_winds=False)

    def dibujar(self) -> None:
        """Legacy v1 method: display plot."""
        self.draw_isobars()
