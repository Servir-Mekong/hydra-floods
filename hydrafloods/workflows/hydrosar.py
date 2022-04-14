import os
import ee
import time
from ee.ee_exception import EEException
import gcsfs
import logging
import datetime
from functools import partial
import warnings
import hydrafloods as hf
from hydrafloods import (
    datasets,
    timeseries,
    filtering,
    ml,
    utils,
    geeutils,
    indices,
    thresholding,
    decorators,
    corrections,
    fuzzy,
)


# set global constant variables in module that are used throughout
MERIT = ee.Image("MERIT/Hydro/v1_0_1")
DEM = MERIT.select("elv")
HAND = MERIT.select("hnd")
HAND_THRESH = 15  # hand threshold
HAND_MASK = HAND.lt(HAND_THRESH)
SLOPE = ee.Terrain.slope(HAND.unmask(0))
SAND = ee.Image("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02")


def hydro30(
    region,
    target_date,
    speckle_filter="gamma_map",
    grid_size=0.1,
    scale=90,
    combine_scenes=True,
):

    if not isinstance(target_date, ee.Date):
        target_date = ee.Date(target_date)

    s1 = datasets.Sentinel1(region, target_date, target_date.advance(1, "day"))

    filter_opts = dict(
        lee_sigma=filtering.lee_sigma,
        gamma_map=filtering.gamma_map,
        refined_lee=filtering.refined_lee,
    )
    filter_func = filter_opts[speckle_filter]

    proc = (
        (corrections.slope_correction, dict(elevation=DEM, buffer=60)),
        filter_func,
        (_water_classification, dict(grid_size=grid_size, scale=scale)),
    )

    processed = s1.pipe(proc)

    if combine_scenes:
        output = processed.collection.mosaic()
    else:
        output = processed.collection

    return output


@decorators.keep_attrs
def _water_classification(img, grid_size=0.1, scale=90):
    def _per_band_water(img, tiles, bandname):
        # !!! this function assumes one band image !!!
        def _calc_t(tile):
            histogram = band.updateMask(band.lt(0)).reduceRegion(
                ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
                .combine("mean", None, True)
                .combine("variance", None, True),
                tile.geometry(),
                scale=scale,
                bestEffort=True,
                tileScale=16,
            )

            threshold = hf.otsu(histogram.get(bandname.cat("_histogram")))

            return tile.set("t", threshold)

        max_buckets = 255
        min_bucket_width = 0.001
        max_raw = 1e6
        bandname = ee.String(bandname)
        band = img.select(bandname)

        thresh_tiles = tiles.map(_calc_t)
        db_t = ee.Number(thresh_tiles.aggregate_array("t").reduce(ee.Reducer.mean()))

        water_init = band.lt(db_t)

        # db_stats = band.updateMask(band.lt(0)).reduceRegion(
        #     reducer=ee.Reducer.median(),
        #     geometry = region,
        #     scale=scale,
        #     bestEffort=True,
        #     tileScale=8
        # )
        #
        # db_lower = db_stats.get(bandname)
        # db_mid = db_t.add(db_lower).divide(2)

        fuzzy_db = fuzzy.fuzzy_large(band.select(bandname), db_t, 5)

        connected_objects = (
            water_init.select(bandname)
            .selfMask()
            .connectedComponents(connectedness=ee.Kernel.square(1.5), maxSize=1024)
        )

        connected_n = connected_objects.select("labels").connectedPixelCount(
            maxSize=1024
        )

        fuzzy_connected = fuzzy.fuzzy_zmf(connected_n.unmask(0), 10, 3)

        fuzziness = [
            fuzzy_slope,
            fuzzy_hand,
            fuzzy_db.select(bandname),
            fuzzy_connected,
            fuzzy_sand,
        ]

        weights = [1 / len(fuzziness) for _ in range(len(fuzziness))]

        fuzzy_combo = hf.fuzzy_weighted(fuzziness, weights).updateMask(
            band.select([0]).mask()
        )

        fuzzy_mask = fuzzy_combo.gt(0.60)

        water_final = water_init.And(fuzzy_mask)

        return water_final.rename(bandname.cat("_water")).uint8()

    def tile_stats(tile):
        stats = img.reduceRegion(
            reducer=stat_reducer,
            geometry=tile.geometry(),
            scale=scale,
            bestEffort=True,
            tileScale=16,
        )

        vv_denom = ee.Number(stats.get("VV_stdDev"))
        vh_denom = ee.Number(stats.get("VH_stdDev"))

        vv_denom = ee.Algorithms.If(vv_denom, vv_denom, ee.Number(0.00001))
        vh_denom = ee.Algorithms.If(vh_denom, vh_denom, ee.Number(0.00001))

        vv_ratio = ee.Number(img_stats.get("VV_p50")).divide(vv_denom)
        vh_ratio = ee.Number(img_stats.get("VH_p50")).divide(vh_denom)

        stats = stats.set("VV_ratio", vv_ratio).set("VH_ratio", vh_ratio)

        return tile.set(stats)

    #

    db_bands = img.bandNames().removeAll(["angle", "local_inc_angle"])
    hd_bands = ["hnd", "hand_mask"]
    img = img.select(db_bands).addBands(HAND_MASK).addBands(HAND)

    region = img.geometry()

    # this only creates tiles that are within the image so some tiles along edges will drop
    tiles = hf.tile_region(region, contain_geom=region, grid_size=grid_size)

    stat_reducer = ee.Reducer.stdDev().combine(ee.Reducer.mean(), None, True)

    img_stats = img.reduceRegion(
        reducer=ee.Reducer.percentile([50, 90]),
        geometry=region,
        scale=scale,
        bestEffort=True,
        tileScale=16,
    )

    stat_tiles = tiles.map(tile_stats).filter(ee.Filter.gte("hnd_mean", 0.8))

    vv_sel_tiles = stat_tiles.limit(5, "VV_ratio", False)
    vh_sel_tiles = stat_tiles.limit(5, "VH_ratio", False)

    # calculate the fuzzy membership of HAND Slope
    fuzzy_slope = fuzzy.fuzzy_zmf(SLOPE, 0, 15)

    hand_stats = HAND.updateMask(
        HAND.neq(0).And(HAND.lt(ee.Number(img_stats.get("hnd_p90"))))
    ).reduceRegion(
        reducer=ee.Reducer.median().combine(ee.Reducer.stdDev(), None, True),
        geometry=region,
        scale=scale,
        bestEffort=True,
        tileScale=16,
    )

    hand_lower = ee.Number(hand_stats.get("hnd_median"))
    hand_std = ee.Number(hand_stats.get("hnd_stdDev"))
    hand_upper = hand_lower.add(ee.Number(3.0).multiply(hand_std)).add(5.0)

    fuzzy_hand = fuzzy.fuzzy_zmf(HAND, hand_lower, hand_upper).unmask(1)

    fuzzy_sand = fuzzy.fuzzy_small(SAND.select([0]).unmask(0), 37.5, 8)

    vv_water = _per_band_water(img, vv_sel_tiles, "VV")
    vh_water = _per_band_water(img, vh_sel_tiles, "VH")

    overlap = vv_water.Or(vh_water)

    # post processing opening of the segmented water
    opened_closed_water = filtering.close_binary(filtering.open_binary(overlap))

    out_water = opened_closed_water.rename("water").updateMask(img.select([0]).mask())

    return out_water
