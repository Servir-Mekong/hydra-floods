import pytest

import ee
import hydrafloods as hf

ee.Initialize()

TEST_REGION = ee.Geometry.Rectangle([102.3335, 10.4045, 107.6277, 14.6900])
# define two time periods to test for newer/older data
TEST_START_TIME_E2 = "2022-02-07"
TEST_END_TIME_E2 = "2022-02-14"

TEST_START_TIME_E1 = "2009-02-07"
TEST_END_TIME_E1 = "2009-02-14"



class TestDatasets:
    def test_sentinel1(self):
        s1 = hf.Sentinel1(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        assert s1.n_images == 20

    def test_sentinel2(self):
        s2 = hf.Sentinel2(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        assert s2.n_images == 100

    def test_sentinel2_dedupe(self):
        s2 = hf.Sentinel2(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        s2_dedupe = s2.deduplicate()
        assert s2_dedupe.n_images == 5

    def test_landsat9(self):
        lc9 = hf.Landsat9(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        assert lc9.n_images == 11

    def test_landsat8(self):
        lc8 = hf.Landsat8(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        assert lc8.n_images == 8

    def test_landsat7(self):
        le7 = hf.Landsat7(TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1)
        assert le7.n_images == 5

    def test_landsat5(self):
        lt5 = hf.Landsat5(TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1)
        assert lt5.n_images == 7

    def test_modis_terra(self):
        modis_t = hf.Modis(TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1)
        assert modis_t.n_images == 7

    def test_modis_aqua(self):
        modis_a = hf.Modis(TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1, asset_id="MODIS/006/MYD09GA")
        assert modis_a.n_images == 7

    def test_viirs(self):
        viirs = hf.Viirs(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        assert viirs.n_images == 7

    def test_from_imgcollection(self):
        us = hf.country_bbox("United States")

        # get the public SkySat ImageCollection for the USA
        # and compute a NDVI band
        ic = ee.ImageCollection("SKYSAT/GEN-A/PUBLIC/ORTHO/MULTISPECTRAL").filterBounds(
            us
        )

        planet = hf.Dataset.from_imgcollection(ic)

        assert planet.n_images == 33

    def test_merge(self):
        lc8 = hf.Landsat8(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)
        viirs = hf.Viirs(TEST_REGION, TEST_START_TIME_E2, TEST_END_TIME_E2)

        merged = lc8.merge(viirs)
        assert merged.n_images == 15

    def test_join(self):
        TEST_REGION_JOIN = ee.Geometry.Rectangle([0,0,120,30],"EPSG:4326", False)
        lc8 = hf.Landsat8(TEST_REGION_JOIN, TEST_START_TIME_E2, TEST_END_TIME_E2)
        s1 = hf.Sentinel1(TEST_REGION_JOIN, TEST_START_TIME_E2, TEST_END_TIME_E2)

        joined = lc8.join(s1)
        first = joined.collection.first()

        assert (joined.n_images == 105) and (
            first.bandNames().getInfo()
            == ["blue", "green", "red", "nir", "swir1", "swir2", "VV", "VH", "angle"]
        )

    def test_temporalagg(self):
        terra = hf.Modis(TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1)
        # get the aqua MODIS dataset
        # note calling the asset_id explicitly
        aqua = hf.Modis(
            TEST_REGION, TEST_START_TIME_E1, TEST_END_TIME_E1, asset_id="MODIS/006/MYD09GA"
        )

        merged = terra.merge(aqua)

        year_starts = [
            "2015-01-01",
            "2016-01-01",
            "2017-01-01",
            "2018-01-01",
            "2019-01-01",
        ]
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
