import pytest

import ee
import hydrafloods as hf

ee.Initialize()

TEST_REGION = ee.Geometry.Rectangle([102.3335, 10.4045, 107.6277, 14.6900])
TEST_START_TIME = "2019-02-07"
TEST_END_TIME = "2019-02-08"

S1 = ee.Image(
    hf.Sentinel1(TEST_REGION, TEST_START_TIME, TEST_END_TIME).collection.first()
)

GEOM = S1.geometry().centroid()


def test_lee_sigma():
    filtered = hf.lee_sigma(S1)
    value = filtered.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected = {
        "VH": -13.66781582837678,
        "VV": -7.876199146226998,
        "angle": 38.63777542114258,
    }

    assert value.getInfo() == expected


def test_gamma_map():
    filtered = hf.gamma_map(S1)
    value = filtered.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected = {
        "VH": -13.977827072143555,
        "VV": -7.706578731536865,
        "angle": 38.63688659667969,
    }

    assert value.getInfo() == expected


def test_refined_lee():
    filtered = hf.refined_lee(S1)
    value = filtered.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected = {
        "VH": -14.386202713274955,
        "VV": -7.370427354980791,
        "angle": 38.63797261920891,
    }

    assert value.getInfo() == expected


def test_p_median():
    filtered = hf.p_median(S1)
    value = filtered.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected = {
        "VH": -14.252587399378688,
        "VV": -7.465488774911076,
        "angle": 38.63777542114258,
    }

    assert value.getInfo() == expected


def test_perona_malik():
    filtered1 = hf.perona_malik(S1, method=1)
    value1 = filtered1.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected1 = {
        "VH": -13.829632515470857,
        "VV": -7.3728422483001115,
        "angle": 38.63777451159317,
    }

    filtered2 = hf.perona_malik(S1, method=2)
    value2 = filtered2.reduceRegion(ee.Reducer.mean(), GEOM, 10)

    expected2 = {
        "VH": -13.829609585803633,
        "VV": -7.372836610955799,
        "angle": 38.637774542076784,
    }

    assert (value1.getInfo() == expected1) and (value2.getInfo() == expected2)
