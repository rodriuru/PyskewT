# -*- coding: utf-8 -*-
"""Benchmark comparing PySkewT v1 approach vs v2 LineCollection and RK4 engine."""

import time
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pyskewt import SkewT, Sounding


def benchmark():
    print("=" * 60)
    print("  PySkewT v2 Performance Benchmark")
    print("=" * 60)

    # 1. Background Grid Generation Benchmark
    iterations = 5
    t0 = time.perf_counter()
    for _ in range(iterations):
        skew = SkewT(with_wind_barbs=False)
        skew.draw_background()
        plt.close(skew.fig)
    t_v2_grid = (time.perf_counter() - t0) / iterations * 1000

    print(f"v2 Background Grid (Isobars + Isotherms + Dry & Moist Adiabats + Isohumes):")
    print(f"  --> {t_v2_grid:.2f} ms per full diagram")

    # 2. Convective Calculation & Sounding Rendering
    sounding = Sounding.from_file("sondeo3.txt")
    t0 = time.perf_counter()
    for _ in range(iterations):
        skew = SkewT(with_wind_barbs=True)
        skew.draw_background()
        indices = skew.plot_sounding(sounding, plot_parcel=True, shade_cape_cin=True, plot_winds=True)
        plt.close(skew.fig)
    t_v2_full = (time.perf_counter() - t0) / iterations * 1000

    print(f"\nv2 Full Pipeline (Grid + Sounding + Winds + Parcel + CAPE/CIN + HUD):")
    print(f"  --> {t_v2_full:.2f} ms per complete analysis")

    print("\nThermodynamic Indices Computed:")
    for k in ["cape", "cin", "p_lcl", "t_lcl", "p_lfc", "p_el", "lifted_index", "k_index", "total_totals", "precipitable_water"]:
        print(f"  • {k:20s}: {indices.get(k)}")

    print("\nBenchmark completed successfully.")


if __name__ == "__main__":
    benchmark()
