# -*- coding: utf-8 -*-
"""Unit tests for the theme and color palette system in pyskewt.theme."""

import os
import unittest
import matplotlib.pyplot as plt
from pyskewt.theme import Theme, get_theme, list_themes, THEMES
from pyskewt.skewt import SkewT
from pyskewt.sounding import Sounding


class TestTheme(unittest.TestCase):

    def setUp(self):
        self.tmp_output = "temp_theme_test.png"

    def tearDown(self):
        if os.path.exists(self.tmp_output):
            os.remove(self.tmp_output)
        plt.close("all")

    def test_list_themes(self):
        canonical = list_themes()
        self.assertEqual(len(canonical), 4)
        self.assertIn("modern_light", canonical)
        self.assertIn("dark_mode", canonical)
        self.assertIn("classic_meteo", canonical)
        self.assertIn("accessible", canonical)

    def test_get_theme_valid_and_aliases(self):
        t_modern = get_theme("modern")
        self.assertEqual(t_modern.name, "modern_light")

        t_dark = get_theme("dark")
        self.assertEqual(t_dark.name, "dark_mode")

        t_night = get_theme("night")
        self.assertEqual(t_night.name, "dark_mode")

        t_classic = get_theme("classic")
        self.assertEqual(t_classic.name, "classic_meteo")

        t_acc = get_theme("accessible")
        self.assertEqual(t_acc.name, "accessible")

        t_cb = get_theme("colorblind")
        self.assertEqual(t_cb.name, "accessible")

        # Passing Theme instance directly
        self.assertEqual(get_theme(t_modern), t_modern)

    def test_get_theme_invalid(self):
        with self.assertRaises(ValueError):
            get_theme("tema_inexistente")

    def test_skewt_themes_render_and_save(self):
        sounding = Sounding.from_file("sondeo3.txt")

        for theme_name in ["modern", "dark", "classic", "accessible"]:
            skew = SkewT(theme=theme_name, with_wind_barbs=True)
            skew.draw_background()
            indices = skew.plot_sounding(sounding, plot_parcel=True, shade_cape_cin=True, plot_winds=True)
            self.assertIn("cape", indices)

            skew.save(self.tmp_output, dpi=80)
            self.assertTrue(os.path.exists(self.tmp_output))
            self.assertGreater(os.path.getsize(self.tmp_output), 1000)
            os.remove(self.tmp_output)
            plt.close(skew.fig)


if __name__ == "__main__":
    unittest.main()
