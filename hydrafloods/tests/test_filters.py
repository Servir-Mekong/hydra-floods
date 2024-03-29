import pytest

import ee
import hydrafloods as hf

ee.Initialize()

# set constants for testing
IMG = ee.Image("projects/servir-mekong/unitTests/dummyv1")
GEOM = IMG.geometry()
SCALE = 1


class TestFilters:
    def test_lee_sigma(self):
        filtered = hf.lee_sigma(IMG, keep_bands=None)
        result = filtered.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 3) for k, v in result.items()}

        expected = {"b1": 2.751}

        assert result == expected

    def test_gamma_map(self):
        filtered = hf.gamma_map(IMG, keep_bands=None)
        result = filtered.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 3) for k, v in result.items()}

        expected = {"b1": 0.518}

        assert result == expected

    def test_refined_lee(self):
        filtered = hf.refined_lee(IMG, keep_bands=None)
        result = filtered.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 3) for k, v in result.items()}

        expected = {"b1": 0.514}

        assert result == expected

    def test_p_median(self):
        filtered = hf.p_median(IMG)
        result = filtered.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 3) for k, v in result.items()}

        expected = {"b1": 0.503}

        assert result == expected

    def test_perona_malik(self):
        filtered1 = hf.perona_malik(IMG, method=1)
        result1 = filtered1.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result1 = {k: round(v, 3) for k, v in result1.items()}

        expected1 = {"b1": -0.253}

        filtered2 = hf.perona_malik(IMG, method=2)
        result2 = filtered2.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result2 = {k: round(v, 3) for k, v in result2.items()}

        expected2 = {"b1": -0.328}

        assert (result1 == expected1) and (result2 == expected2)
