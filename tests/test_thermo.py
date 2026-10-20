# -*- coding: utf-8 -*-
"""Unit tests for atmospheric thermodynamics in pyskewt.thermo."""

import unittest
import numpy as np
from pyskewt.thermo import (
    celsius_to_kelvin,
    kelvin_to_celsius,
    saturation_vapor_pressure,
    dewpoint_from_vapor_pressure,
    saturation_mixing_ratio,
    potential_temperature,
    temperature_at_potential_temperature,
    lifting_condensation_level,
    moist_adiabat_profile,
    calculate_convective_indices,
)


class TestThermo(unittest.TestCase):

    def test_temperature_conversions(self):
        self.assertAlmostEqual(celsius_to_kelvin(0.0), 273.15, places=2)
        self.assertAlmostEqual(celsius_to_kelvin(25.0), 298.15, places=2)
        self.assertAlmostEqual(kelvin_to_celsius(273.15), 0.0, places=2)

    def test_saturation_vapor_pressure(self):
        # At 0 deg C, Bolton es is ~ 6.112 hPa
        es_0 = saturation_vapor_pressure(0.0)
        self.assertAlmostEqual(float(es_0), 6.112, places=2)

        # At 20 deg C, standard es is ~ 23.37 hPa
        es_20 = saturation_vapor_pressure(20.0)
        self.assertAlmostEqual(float(es_20), 23.37, delta=0.1)

        # At 35 deg C, es is ~ 56.2 hPa
        es_35 = saturation_vapor_pressure(35.0)
        self.assertAlmostEqual(float(es_35), 56.24, delta=0.2)

    def test_dewpoint_inversion(self):
        temps = [-20.0, -5.0, 0.0, 15.0, 30.0]
        for t in temps:
            e = saturation_vapor_pressure(t)
            td = dewpoint_from_vapor_pressure(e)
            self.assertAlmostEqual(float(td), t, places=2)

    def test_potential_temperature(self):
        # At 1000 hPa, theta = T (in Kelvin)
        theta = potential_temperature(1000.0, 20.0)
        self.assertAlmostEqual(float(theta), 293.15, places=2)

        # Inversion check
        t_back = temperature_at_potential_temperature(theta, 1000.0)
        self.assertAlmostEqual(float(t_back), 20.0, places=2)

    def test_lifting_condensation_level(self):
        # Saturated parcel (T == Td) should condense immediately at surface
        p_lcl, t_lcl = lifting_condensation_level(1000.0, 20.0, 20.0)
        self.assertAlmostEqual(p_lcl, 1000.0, delta=1.0)
        self.assertAlmostEqual(t_lcl, 20.0, delta=0.5)

        # Unsaturated parcel
        p_lcl2, t_lcl2 = lifting_condensation_level(1014.0, 28.2, 25.2)
        self.assertTrue(960.0 < p_lcl2 < 980.0)
        self.assertTrue(23.0 < t_lcl2 < 26.0)

    def test_moist_adiabat_cooling(self):
        p_levels = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
        t_prof = moist_adiabat_profile(20.0, p_levels)
        # Saturated parcel should get colder as pressure drops
        self.assertTrue(np.all(np.diff(t_prof) < 0))
        self.assertAlmostEqual(t_prof[0], 20.0, places=1)


if __name__ == "__main__":
    unittest.main()
