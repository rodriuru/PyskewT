# -*- coding: utf-8 -*-
"""Atmospheric sounding data container and parser for PySkewT v2."""

from __future__ import annotations
import os
import re
from typing import Optional, Union, Dict, Any
import numpy as np


class Sounding:
    """Represents an atmospheric vertical sounding (radiosonde observation)."""

    def __init__(
        self,
        pressure: np.ndarray,
        temperature: np.ndarray,
        dewpoint: np.ndarray,
        height: Optional[np.ndarray] = None,
        wind_speed: Optional[np.ndarray] = None,
        wind_direction: Optional[np.ndarray] = None,
        station_id: Optional[str] = None,
        datetime_str: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        p = np.asarray(pressure, dtype=float)
        t = np.asarray(temperature, dtype=float)
        td = np.asarray(dewpoint, dtype=float)

        # Filter out NaN or clearly corrupt values (-999, etc.)
        valid_mask = (
            np.isfinite(p) & np.isfinite(t) & np.isfinite(td) &
            (p > 0) & (p <= 1200) & (t > -120) & (t < 70) &
            (td > -150) & (td <= t + 5.0)
        )

        p = p[valid_mask]
        t = t[valid_mask]
        td = td[valid_mask]

        # Sort by pressure descending (surface -> tropopause)
        sort_idx = np.argsort(p)[::-1]
        self.pressure = p[sort_idx]
        self.temperature = t[sort_idx]
        self.dewpoint = td[sort_idx]

        if height is not None:
            h = np.asarray(height, dtype=float)[valid_mask][sort_idx]
            self.height = h
        else:
            self.height = None

        if wind_speed is not None and wind_direction is not None:
            ws = np.asarray(wind_speed, dtype=float)[valid_mask][sort_idx]
            wd = np.asarray(wind_direction, dtype=float)[valid_mask][sort_idx]
            # Replace invalid winds
            ws_valid = np.where((ws >= 0) & (ws < 300), ws, np.nan)
            wd_valid = np.where((wd >= 0) & (wd <= 360), wd, np.nan)
            self.wind_speed = ws_valid
            self.wind_direction = wd_valid

            # Meteorological wind components u (zonal, east) and v (meridional, north)
            # Wind direction is direction FROM which wind blows
            rad = np.radians(wd_valid)
            self.u = -ws_valid * np.sin(rad)
            self.v = -ws_valid * np.cos(rad)
        else:
            self.wind_speed = None
            self.wind_direction = None
            self.u = None
            self.v = None

        self.station_id = station_id
        self.datetime_str = datetime_str
        self.metadata = metadata or {}

    def __len__(self) -> int:
        return len(self.pressure)

    @classmethod
    def from_file(cls, filepath: Union[str, os.PathLike]) -> "Sounding":
        """Load a sounding from a file (University of Wyoming format or tabular CSV/whitespace)."""
        filepath = str(filepath)
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        data_rows = []
        station_id = None
        datetime_str = None

        # Look for headers or metadata in Wyoming format
        for line in lines:
            line_str = line.strip()
            if not line_str:
                continue

            # Check for header metadata like "72202 MFL Miami Observations at 00Z 15 Dec 2021"
            if "Observations at" in line or "Station number" in line:
                datetime_str = line_str

            parts = line_str.split()
            # If line has at least 3 numeric elements
            if len(parts) >= 3:
                try:
                    p = float(parts[0])
                    # Handle typical Wyoming layout:
                    # PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT
                    if len(parts) >= 8:
                        h = float(parts[1])
                        t = float(parts[2])
                        td = float(parts[3])
                        drct = float(parts[6])
                        sknt = float(parts[7])
                        data_rows.append((p, t, td, h, sknt, drct))
                    elif len(parts) >= 4:
                        # Alternative 4-column: p, h, t, td or p, t, td
                        t = float(parts[2])
                        td = float(parts[3])
                        data_rows.append((p, t, td, None, None, None))
                    elif len(parts) == 3:
                        t = float(parts[1])
                        td = float(parts[2])
                        data_rows.append((p, t, td, None, None, None))
                except ValueError:
                    continue

        if not data_rows:
            raise ValueError(f"Could not parse valid sounding data from {filepath}")

        arr = np.array(data_rows, dtype=object)
        p = np.array(arr[:, 0], dtype=float)
        t = np.array(arr[:, 1], dtype=float)
        td = np.array(arr[:, 2], dtype=float)
        h = np.array(arr[:, 3], dtype=float) if arr[0, 3] is not None else None
        ws = np.array(arr[:, 4], dtype=float) if arr[0, 4] is not None else None
        wd = np.array(arr[:, 5], dtype=float) if arr[0, 5] is not None else None

        return cls(
            pressure=p,
            temperature=t,
            dewpoint=td,
            height=h,
            wind_speed=ws,
            wind_direction=wd,
            station_id=station_id,
            datetime_str=datetime_str,
        )

    @classmethod
    def from_arrays(
        cls,
        pressure: np.ndarray,
        temperature: np.ndarray,
        dewpoint: np.ndarray,
        height: Optional[np.ndarray] = None,
        wind_speed: Optional[np.ndarray] = None,
        wind_direction: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "Sounding":
        """Instantiate a Sounding directly from NumPy arrays."""
        return cls(
            pressure=pressure,
            temperature=temperature,
            dewpoint=dewpoint,
            height=height,
            wind_speed=wind_speed,
            wind_direction=wind_direction,
            **kwargs,
        )
