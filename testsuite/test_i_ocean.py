"""Tests for i.ocean."""

import grass.script as gs
from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestIOceanBasic(TestCase):
    """Test basic i.ocean functionality."""

    input_map = "test_iocean_input"
    output_map = "test_iocean_output"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        # Small 20×20 region (geographic, lat ~35°N → temperate style)
        cls.runModule(
            "g.region",
            n=35.1, s=34.9, e=10.1, w=9.9,
            rows=20, cols=20,
        )
        # Create a simple mask: left half = ocean (1), right half = land (null)
        gs.mapcalc(
            f"{cls.input_map} = if(col() <= 10, 1, null())",
            overwrite=True, quiet=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        cls.runModule(
            "g.remove", flags="f", type="raster",
            name=[cls.input_map], quiet=True,
        )

    def tearDown(self):
        self.runModule(
            "g.remove", flags="f", type="raster",
            name=[self.output_map], quiet=True,
        )

    # ── basic run ────────────────────────────────────────────────────────────

    def test_output_created(self):
        """Module runs without error and creates the output map."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
        )
        self.assertRasterExists(self.output_map)

    def test_output_null_outside_ocean(self):
        """Non-ocean pixels must be null in the output."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        # The right half is null → total non-null pixels = 20*10 = 200
        self.assertEqual(int(float(stats["n"])), 200)

    def test_flat_depth_default_value(self):
        """Without depth input the depth index should be 500 everywhere."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        self.assertAlmostEqual(float(stats["mean"]), 500.0, places=1)
        self.assertAlmostEqual(float(stats["stddev"]), 0.0, places=1)

    # ── style / palette ──────────────────────────────────────────────────────

    def test_explicit_style_accepted(self):
        """Explicit style option is accepted without error."""
        for style in ("tropical", "subtropical", "temperate", "subpolar", "polar"):
            with self.subTest(style=style):
                self.assertModule(
                    "i.ocean",
                    input=self.input_map,
                    output=self.output_map,
                    style=style,
                    overwrite=True,
                )
                self.assertRasterExists(self.output_map)

    def test_explicit_latitude(self):
        """Supplying an explicit latitude does not crash the module."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            latitude="-60",  # subpolar southern hemisphere
        )
        self.assertRasterExists(self.output_map)

    # ── flags ────────────────────────────────────────────────────────────────

    def test_wave_flag_adds_variance(self):
        """Wave texture (-w) should introduce non-zero standard deviation."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            flags="w",
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        self.assertGreater(float(stats["stddev"]), 0.0)

    def test_shore_depth_flag(self):
        """Shore-distance depth (-s) should produce a non-uniform depth map."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            flags="s",
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        self.assertGreater(float(stats["stddev"]), 0.0)

    def test_lat_gradient_flag(self):
        """Latitude gradient (-l) should introduce north–south variation."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            flags="l",
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        self.assertGreater(float(stats["stddev"]), 0.0)

    def test_all_flags_combined(self):
        """All three flags together should run without error."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            flags="wsl",
        )
        self.assertRasterExists(self.output_map)

    # ── bad input ────────────────────────────────────────────────────────────

    def test_wrong_ocean_value_fails(self):
        """A wrong ocean_value that matches no pixel must raise a fatal error."""
        self.assertModuleFail(
            "i.ocean",
            input=self.input_map,
            output=self.output_map,
            ocean_value="99",
        )

    def test_no_input_no_depth_fails(self):
        """Omitting both input and depth must be rejected."""
        self.assertModuleFail(
            "i.ocean",
            output=self.output_map,
        )


class TestIOceanWithDepth(TestCase):
    """Test i.ocean with a user-supplied depth map alongside an input mask."""

    input_map = "test_iocean_input_d"
    depth_map = "test_iocean_depth"
    output_map = "test_iocean_output_d"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule(
            "g.region",
            n=35.1, s=34.9, e=10.1, w=9.9,
            rows=20, cols=20,
        )
        gs.mapcalc(
            f"{cls.input_map} = if(col() <= 10, 1, null())",
            overwrite=True, quiet=True,
        )
        # Synthetic depth: increases with column number (shallow near shore)
        gs.mapcalc(
            f"{cls.depth_map} = if(col() <= 10, col() * 200.0, null())",
            overwrite=True, quiet=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        cls.runModule(
            "g.remove", flags="f", type="raster",
            name=[cls.input_map, cls.depth_map], quiet=True,
        )

    def tearDown(self):
        self.runModule(
            "g.remove", flags="f", type="raster",
            name=[self.output_map], quiet=True,
        )

    def test_depth_map_produces_gradient(self):
        """A varying depth map should produce a range of depth-index values."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            depth=self.depth_map,
            output=self.output_map,
        )
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        self.assertGreater(float(stats["stddev"]), 10.0)
        # Index must stay in [0, 1000]
        self.assertGreaterEqual(float(stats["min"]), 0.0)
        self.assertLessEqual(float(stats["max"]), 1000.0)

    def test_depth_map_with_waves(self):
        """Depth map + waves should still produce valid output."""
        self.assertModule(
            "i.ocean",
            input=self.input_map,
            depth=self.depth_map,
            output=self.output_map,
            flags="w",
        )
        self.assertRasterExists(self.output_map)


class TestIOceanDepthOnly(TestCase):
    """Test i.ocean using the depth map as the sole input (no mask raster)."""

    depth_map = "test_iocean_depth_only"
    output_map = "test_iocean_output_only"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule(
            "g.region",
            n=35.1, s=34.9, e=10.1, w=9.9,
            rows=20, cols=20,
        )
        # Full-coverage depth map: col 1-10 = ocean (depth > 0), col 11-20 = 0
        gs.mapcalc(
            f"{cls.depth_map} = if(col() <= 10, col() * 200.0, 0.0)",
            overwrite=True, quiet=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        cls.runModule(
            "g.remove", flags="f", type="raster",
            name=[cls.depth_map], quiet=True,
        )

    def tearDown(self):
        self.runModule(
            "g.remove", flags="f", type="raster",
            name=[self.output_map], quiet=True,
        )

    def test_depth_as_sole_input(self):
        """Depth map alone (no input mask) should produce valid ocean output."""
        self.assertModule(
            "i.ocean",
            depth=self.depth_map,
            output=self.output_map,
            depth_min="1",  # only cells with depth >= 1 are ocean
        )
        self.assertRasterExists(self.output_map)
        stats = gs.parse_command(
            "r.univar", map=self.output_map, flags="g", quiet=True
        )
        # Only the 10 ocean columns should be non-null (10 * 20 = 200 pixels)
        self.assertEqual(int(float(stats["n"])), 200)

    def test_depth_only_with_waves_and_shore(self):
        """Depth-only input with flags -w -s -l should complete without error."""
        self.assertModule(
            "i.ocean",
            depth=self.depth_map,
            output=self.output_map,
            depth_min="1",
            flags="wsl",
        )
        self.assertRasterExists(self.output_map)

    def test_depth_only_bad_threshold_fails(self):
        """A depth_min higher than any value in the depth map must be rejected."""
        self.assertModuleFail(
            "i.ocean",
            depth=self.depth_map,
            output=self.output_map,
            depth_min="99999",
        )


if __name__ == "__main__":
    test()
