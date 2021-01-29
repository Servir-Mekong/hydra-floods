import pytest

import ee
import hydrafloods as hf

ee.Initialize()

# set constants for testing
IMG = ee.Image("projects/servir-mekong/unitTests/dummyv1")
GEOM = IMG.geometry()
SCALE = 1

class TestFuzzy:
    def test_fuzzy_gaussian(self):
        fuzzed = hf.fuzzy_gaussian(IMG,0.25,0.1)
        result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {"b1": 0.985699}

        assert result == expected

    def test_fuzzy_near(self):
        fuzzed = hf.fuzzy_near(IMG,0.25,0.1)
        result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {'b1': 0.985921}

        assert result == expected

    def test_fuzzy_large(self):
        fuzzed = hf.fuzzy_large(IMG,0.25,5)
        result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {'b1': 0.746005}

        assert result == expected

    # def test_fuzzy_mslarge(self):
    #     fuzzed = hf.fuzzy_mslarge(IMG,1,1)
    #     result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
    #     result = {k: round(v, 6) for k, v in result.items()}

    #     expected = {'b1': 0.968378}

    #     assert result == expected
    
    def test_fuzzy_small(self):
        fuzzed = hf.fuzzy_small(IMG,0.25,5)
        result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {'b1': 0.253995}

        assert result == expected

    # def test_fuzzy_mssmall(self):
    #     fuzzed = hf.fuzzy_mssmall(IMG,1,1)
    #     result = fuzzed.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
    #     result = {k: round(v, 6) for k, v in result.items()}

    #     expected = {'b1': 0.968378}

    #     assert result == expected

    def test_fuzzy_linear(self):
        fuzzed1 = hf.fuzzy_linear(IMG,-1,1)
        result1 = fuzzed1.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result1 = {k: round(v, 6) for k, v in result1.items()}

        fuzzed2 = hf.fuzzy_linear(IMG,1.5,-1)
        result2 = fuzzed2.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result2 = {k: round(v, 6) for k, v in result2.items()}

        expected1 = {'b1': 0.752007}
        expected2 = {'b1': 0.398394}

        assert (result1 == expected1) and (result2 == expected2)

    def test_fuzzy_or(self):
        f1 = hf.fuzzy_gaussian(IMG,0.25,0.1)
        f2 = hf.fuzzy_small(IMG,0.25,5)
        overlay = hf.fuzzy_or([f1,f2])
        result = overlay.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {"fuzzy_or": 0.985997}

        assert result == expected

    def test_fuzzy_and(self):
        f1 = hf.fuzzy_gaussian(IMG,0.25,0.1)
        f2 = hf.fuzzy_small(IMG,0.25,5)
        overlay = hf.fuzzy_and([f1,f2])
        result = overlay.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {"fuzzy_and": 0.253697}

        assert result == expected

    def test_fuzzy_gamma(self):
        f1 = hf.fuzzy_gaussian(IMG,0.25,0.1)
        f2 = hf.fuzzy_small(IMG,0.25,5)
        overlay = hf.fuzzy_gamma([f1,f2],gamma=0.4)
        result = overlay.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {"fuzzy_gamma": 0.322515}

        assert result == expected

    def test_fuzzy_weighted(self):
        f1 = hf.fuzzy_gaussian(IMG,0.25,0.1)
        f2 = hf.fuzzy_small(IMG,0.25,5)
        overlay = hf.fuzzy_weighted([f1,f2],[0.3,0.7])
        result = overlay.reduceRegion(ee.Reducer.mean(), GEOM, SCALE).getInfo()
        result = {k: round(v, 6) for k, v in result.items()}

        expected = {"fuzzy_weighted": 0.473506}

        assert result == expected

    