# -*- coding: utf-8 -*-
"""Unit tests for SkewT diagram rendering and export."""

import os
import unittest
import matplotlib.pyplot as plt
from pyskewt.skewt import SkewT
from pyskewt.sounding import Sounding


class TestSkewT(unittest.TestCase):

    def setUp(self):
        self.tmp_output = "temp_test_skewt.png"

    def tearDown(self):
        if os.path.exists(self.tmp_output):
            os.remove(self.tmp_output)
        plt.close("all")

    def test_diagram_creation_and_save(self):
        skew = SkewT(with_wind_barbs=True)
        skew.draw_background()

        sounding = Sounding.from_file("sondeo3.txt")
        indices = skew.plot_sounding(
            sounding, plot_parcel=True, shade_cape_cin=True, plot_winds=True
        )

        self.assertIn("cape", indices)
        self.assertIn("cin", indices)
        self.assertGreater(indices["cape"], 0)

        skew.save(self.tmp_output, dpi=100)
        self.assertTrue(os.path.exists(self.tmp_output))
        self.assertGreater(os.path.getsize(self.tmp_output), 1000)

    def test_legacy_interface(self):
        from SkewT import SkewT as LegacySkewT
        leg = LegacySkewT(1)
        leg.T()
        leg.ws()
        leg.TP()
        leg.tpe()
        self.assertIsNotNone(leg.fig)


if __name__ == "__main__":
    unittest.main()
