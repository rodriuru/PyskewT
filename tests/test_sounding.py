# -*- coding: utf-8 -*-
"""Unit tests for Sounding container and parser."""

import unittest
import numpy as np
from pyskewt.sounding import Sounding


class TestSounding(unittest.TestCase):

    def test_from_file_wyoming(self):
        sounding = Sounding.from_file("sondeo3.txt")
        self.assertGreater(len(sounding), 50)
        self.assertAlmostEqual(sounding.pressure[0], 1014.0, places=1)
        self.assertAlmostEqual(sounding.temperature[0], 28.2, places=1)
        self.assertAlmostEqual(sounding.dewpoint[0], 25.2, places=1)

        # Check wind components
        self.assertIsNotNone(sounding.u)
        self.assertIsNotNone(sounding.v)
        self.assertEqual(len(sounding.u), len(sounding.pressure))

    def test_from_arrays(self):
        p = np.array([1000.0, 850.0, 700.0, 500.0])
        t = np.array([25.0, 15.0, 5.0, -10.0])
        td = np.array([20.0, 10.0, -2.0, -25.0])
        sounding = Sounding.from_arrays(p, t, td)
        self.assertEqual(len(sounding), 4)
        self.assertAlmostEqual(sounding.pressure[0], 1000.0)


if __name__ == "__main__":
    unittest.main()
