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


def test_edge_otsu():
    kwargs = dict(
        initial_threshold=-16, edge_buffer=300, scale=150, return_threshold=True
    )
    threshold = (
        hf.edge_otsu(S1, **kwargs).reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()
    )

    thresh_expected = {"constant": -14.203254920711124}

    kwargs["return_threshold"] = False
    water = (
        hf.edge_otsu(S1, **kwargs).reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()
    )

    water_expected = {"water": 0}

    assert (threshold == thresh_expected) and (water == water_expected)

def test_edge_otsu():
    kwargs = dict(
        initial_threshold=-16, scale=150, return_threshold=True
    )
    threshold = (
        hf.bmax_otsu(S1, **kwargs).reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()
    )

    thresh_expected = {"constant": -14.459723114655215}

    kwargs["return_threshold"] = False
    water = (
        hf.bmax_otsu(S1, **kwargs).reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()
    )

    water_expected = {"water": 0}

    assert (threshold == thresh_expected) and (water == water_expected)

def test_kmeans_extent():
    hand = ee.Image("MERIT/Hydro/v1_0_1").select("hnd")

    water = hf.kmeans_extent(S1,hand,band="VV",initial_threshold=-16).reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()

    water_expected = {"water": 0}

    assert water == water_expected

def test_multidim_semisupervised():
    index_img = hf.add_indices(S1,indices=["ndpi","vv_vh_ratio","nvvi","nvhi"])

    multidim = hf.multidim_semisupervised(
        index_img,
        bands=["VV","VH","ndpi"],
        rank_band="VH",
        ranking="min",
        n_samples=500,
    )

    water_proba = multidim.reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()

    proba_expected = {'water_proba': 0.02067818018845265}

    multidim = hf.multidim_semisupervised(
        index_img,
        bands=["VV","VH","ndpi"],
        rank_band="VH",
        ranking="min",
        n_samples=500,
        proba_threshold=0.5
    )

    water = multidim.reduceRegion(ee.Reducer.mean(), GEOM, 10).getInfo()

    water_expected = {"water": 0}

    assert (water_proba == proba_expected) and (water == water_expected)
