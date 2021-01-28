import pytest

import ee
import hydrafloods as hf

ee.Initialize()

TEST_REGION = ee.Geometry.Rectangle([102.3335, 10.4045, 107.6277, 14.6900])
TEST_START_TIME = "2019-02-07"
TEST_END_TIME = "2019-02-14"


def test_sentinel1():
    s1 = hf.Sentinel1(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert s1.n_images == 24


def test_sentinel2():
    s2 = hf.Sentinel2(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert s2.n_images == 101


def test_landsat8():
    lc8 = hf.Landsat8(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert lc8.n_images == 11


def test_landsat7():
    lc7 = hf.Landsat7(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert lc7.n_images == 7


def test_modis():
    modis = hf.Modis(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert modis.n_images == 6


def test_viirs():
    viirs = hf.Viirs(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    assert viirs.n_images == 7


def test_from_imgcollection():
    us = hf.country_bbox("United States")

    # get the public SkySat ImageCollection for the USA
    # and compute a NDVI band
    ic = ee.ImageCollection("SKYSAT/GEN-A/PUBLIC/ORTHO/MULTISPECTRAL").filterBounds(us)

    planet = hf.Dataset.from_imgcollection(ic)

    assert planet.n_images == 33


def test_merge():
    lc8 = hf.Landsat8(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    s2 = hf.Sentinel2(TEST_REGION, TEST_START_TIME, TEST_END_TIME)

    merged = lc8.merge(s2)
    assert merged.n_images == 112


def test_join():
    lc8 = hf.Landsat8(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    s1 = hf.Sentinel1(TEST_REGION, TEST_START_TIME, TEST_END_TIME)

    joined = lc8.join(s1)
    first = joined.collection.first()

    assert (joined.n_images == 8) and (
        first.bandNames().getInfo()
        == ["blue", "green", "red", "nir", "swir1", "swir2", "VV", "VH", "angle"]
    )


def test_temporalagg():
    terra = hf.Modis(TEST_REGION, TEST_START_TIME, TEST_END_TIME)
    # get the aqua MODIS dataset
    # note calling the asset_id explicitly
    aqua = hf.Modis(
        TEST_REGION, TEST_START_TIME, TEST_END_TIME, asset_id="MODIS/006/MYD09GA"
    )

    merged = terra.merge(aqua)

    year_starts = ["2015-01-01", "2016-01-01", "2017-01-01", "2018-01-01", "2019-01-01"]
    agg = merged.aggregate_time(
        dates=year_starts, period=365, reducer=ee.Reducer.median()
    )

    assert agg.dates == [
        "2015-01-01 00:00:00.000",
        "2016-01-01 00:00:00.000",
        "2017-01-01 00:00:00.000",
        "2018-01-01 00:00:00.000",
        "2019-01-01 00:00:00.000",
    ]

