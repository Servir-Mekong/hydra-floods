import ee
from ee.ee_exception import EEException
import random
from hydrafloods import geeutils, decorators


@decorators.carry_metadata
def bmax_otsu(
    img,
    band=None,
    region=None,
    reduction_scale=90,
    initial_threshold=0,
    invert=False,
    grid_size=0.1,
    bmax_threshold=0.75,
    max_boxes=100,
    seed=7,
    max_buckets=255,
    min_bucket_width=0.001,
    max_raw=1e6,
    return_threshold=False
):

    def calcBmax(feature):
        segment = img  # .clip(feature)
        initial = segment.lt(initial_threshold)
        p1 = ee.Number(
            initial.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=feature.geometry(),
                bestEffort=True,
                scale=reduction_scale,
            ).get(histBand)
        )
        p1 = ee.Number(ee.Algorithms.If(p1, p1, 0.99))
        p2 = ee.Number(1).subtract(p1)

        m = (
            segment.updateMask(initial)
            .rename("m1")
            .addBands(segment.updateMask(initial.Not()).rename("m2"))
        )

        mReduced = m.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=feature.geometry(),
            bestEffort=True,
            scale=reduction_scale,
        )

        m1 = ee.Number(mReduced.get("m1"))
        m2 = ee.Number(mReduced.get("m2"))

        m1 = ee.Number(ee.Algorithms.If(m1, m1, -25))
        m2 = ee.Number(ee.Algorithms.If(m2, m2, 0))

        sigmab = p1.multiply(p2.multiply(m1.subtract(m2).pow(2)))
        sigmat = ee.Number(
            segment.reduceRegion(
                reducer=ee.Reducer.variance(),
                geometry=feature.geometry(),
                bestEffort=True,
                scale=reduction_scale,
            ).get(histBand)
        )
        sigmat = ee.Number(ee.Algorithms.If(sigmat, sigmat, 2))
        bmax = sigmab.divide(sigmat)
        return feature.set({"bmax": bmax})

    if band is None:
        img = img.select([0])
        histBand = ee.String(img.bandNames().get(0))

    else:
        histBand = ee.String(band)
        img = img.select(histBand)

    if region is None:
        region = img.geometry()

    grid = geeutils.tile_region(region,intersect_geom=region,grid_size=0.1)

    bmaxes = (
        grid.map(calcBmax)
        .filter(ee.Filter.gt("bmax", bmax_threshold))
        .randomColumn("random", seed)
    )

    nBoxes = ee.Number(bmaxes.size())
    randomThresh = ee.Number(max_boxes).divide(nBoxes)
    selection = bmaxes.filter(ee.Filter.lt("random", randomThresh))

    histogram = img.reduceRegion(
        ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
        .combine("mean", None, True)
        .combine("variance", None, True),
        selection,
        reduction_scale,
        bestEffort=True,
        tileScale=16,
    )

    threshold = otsu(histogram.get(histBand.cat("_histogram")))

    if return_threshold is True:
        return ee.Image(threshold)
    else:
        water = ee.Image(ee.Algorithms.If(invert, img.gt(threshold), img.lt(threshold)))
        return water.rename("water").uint8()


@decorators.carry_metadata
def edge_otsu(
    img,
    initial_threshold=0,
    canny_threshold=0.05,  # threshold for canny edge detection
    canny_sigma=0,  # sigma value for gaussian filter
    canny_lt=0.05,  # lower threshold for canny detection
    smoothing=100,  # amount of smoothing in meters
    connected_pixels=200,  # maximum size of the neighborhood in pixels
    edge_length=50,  # minimum length of edges from canny detection
    edge_buffer=100,
    band=None,
    region=None,
    reduction_scale=90,
    invert=False,
    seed=7,
    max_buckets=255,
    min_bucket_width=0.001,
    max_raw=1e6,
    return_threshold=False
):

    if band is None:
        img = img.select([0])
        histBand = ee.String(img.bandNames().get(0))

    else:
        histBand = ee.String(band)
        img = img.select(histBand)

    if region is None:
        region = img.geometry()

    binary = img.lt(initial_threshold).rename("binary")

    # get canny edges
    canny = ee.Algorithms.CannyEdgeDetector(binary, canny_threshold, canny_sigma)
    # process canny edges
    connected = (
        canny.mask(canny).lt(canny_lt).connectedPixelCount(connected_pixels, True)
    )
    edges = connected.gte(edge_length)
    edgeBuffer = edges.focal_max(edge_buffer, "square", "meters")

    # mask out areas to get histogram for Otsu
    histogram_image = img.updateMask(edgeBuffer)

    histogram = histogram_image.reduceRegion(
        ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
        .combine("mean", None, True)
        .combine("variance", None, True),
        region,
        reduction_scale,
        bestEffort=True,
        tileScale=16,
    )

    threshold = otsu(histogram.get(histBand.cat("_histogram")))

    if return_threshold is True:
        return threshold
    else:
        water = ee.Image(ee.Algorithms.If(invert, img.gt(threshold), img.lt(threshold)))
        return water.rename("water").uint8()


def otsu(histogram):
    counts = ee.Array(ee.Dictionary(histogram).get("histogram"))
    means = ee.Array(ee.Dictionary(histogram).get("bucketMeans"))
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sums = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean = sums.divide(total)
    indices = ee.List.sequence(1, size)
    # Compute between sum of squares, where each mean partitions the data.

    def bss_function(i):
        aCounts = counts.slice(0, 0, i)
        aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
        aMeans = means.slice(0, 0, i)
        aMean = (
            aMeans.multiply(aCounts)
            .reduce(ee.Reducer.sum(), [0])
            .get([0])
            .divide(aCount)
        )
        bCount = total.subtract(aCount)
        bMean = sums.subtract(aCount.multiply(aMean)).divide(bCount)
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(
            bCount.multiply(bMean.subtract(mean).pow(2))
        )

    bss = indices.map(bss_function)
    output = means.sort(bss).get([-1])
    return output
