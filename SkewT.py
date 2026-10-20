# -*- coding: utf-8 -*-
"""SkewT v2 Backwards-Compatibility Module.

Preserves the original API for legacy scripts while delegating to the modern,
high-performance pyskewt engine.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pyskewt import SkewT as ModernSkewT, Sounding
import cte as C


class SkewT(ModernSkewT):
    """Backwards-compatible SkewT class."""

    def __init__(self, id: int = 1, **kwargs):
        super().__init__(**kwargs)

    def Cp(self, w):
        return 1004.0 * (1.0 + 0.87 * w)

    def Lv(self, T):
        pol = np.poly1d(np.array([8.17422613e-04, -6.00650658e-02, 9.57297089e-01, -2.36500484e+03, 2.50090732e+06]))
        return pol(T)

    def e_sat(self, T):
        return np.exp((self.Lv(T) / C.Rv) * (1.0 / C.T_lt - 1.0 / (T + C.T_lt))) * C.e_lt

    def w(self, T, p):
        return C.epsi * self.e_sat(T) / p

    def pdeTheta(self, T, Tp):
        return np.log(1e5) + C.cp * np.log(T / Tp) / C.Rd


if __name__ == "__main__":
    print("Ejecutando PySkewT (modo compatibilidad v1)...")
    m = SkewT(1)
    m.draw_background()

    archivo = "sondeo3.txt"
    try:
        sounding = Sounding.from_file(archivo)
        indices = m.plot_sounding(sounding, title="Sondeo Atmosférico - PySkewT v2")
        print(f"Diagrama generado con éxito. CAPE: {indices['cape']:.1f} J/kg, LCL: {indices['p_lcl']:.1f} hPa")
    except Exception as e:
        print(f"Error procesando sondeo: {e}")

    output_file = "example_v2.png"
    m.save(output_file, dpi=300)
    print(f"Imagen guardada en {output_file}")
