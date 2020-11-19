import ee
from ee.ee_exception import EEException
import random
from hydrafloods import geeutils, decorators


@decorators.carry_metadata
def bmax_otsu(
    img,
    band=None,
    region=None,
    scale=90,
    initial_threshold=0,
    invert=False,
    grid_size=0.1,
    bmax_threshold=0.75,
    max_boxes=100,
    seed=7,
    max_buckets=255,
    min_bucket_width=0.001,
    max_raw=1e6,
    return_threshold=False,
):
    """Implementation of the B-Max Otsu thresholding algorithm.
    Detailed explanation of algorithm can be found at https://doi.org/10.3390/rs12152469

    args:
        img (ee.Image): input image to thresholding algorithm
        band (str | None,optional): band name to use for thresholding, if set to `None` will use first band in image. default = None
        region (ee.Geometry | None, optional): region to determine threshold, if set to `None` will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90
        initial_threshold (float, optional): initial estimate of water/no-water for estimating the probabilities of classes in segment. default = 0
        invert (bool, optional): boolean switch to determine if to threshold greater than (True) or less than (False). default = False
        grid_size (float, optional): size in decimal degrees to tile image/region to check for bimodality. default = 0.1
        bmax_threshold (float, optional): value 0-1 to determine if a value of bmax is bimodal or not. default = 0.75
        max_boxes (int, optional): maximum number of tiles/boxes to use when determining threshold. default = 100
        seed (int, optional): random number generator seed for randomly selected max_boxes. default = 7
        max_buckets (int, optional): The maximum number of buckets to use when building a histogram; will be rounded up to a power of 2. default = 255
        min_bucket_width (float, optional): The minimum histogram bucket width to allow any power of 2. default = 0.001
        max_raw (int, optional): The number of values to accumulate before building the initial histogram. default = 1e6
        return_threshold (bool, optional): boolean switch, if set to true then function will return threshold number, else thresholded image. default = False

    returns:
        ee.Image: thresholded image based (if return_threshold==False) or threshold value (if return_threshold==True) based on the threshold determined by the algorithm
    """

    def calcBmax(feature):
        """Closure function to calculate Bmax for each feature covering image
        """
        segment = img 
        initial = segment.lt(initial_threshold)
        p1 = ee.Number(
            initial.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=feature.geometry(),
                bestEffort=True,
                scale=scale,
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
            scale=scale,
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
                scale=scale,
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

    grid = geeutils.tile_region(region, intersect_geom=region, grid_size=0.1)

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
        scale,
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
    band=None,
    region=None,
    scale=90,
    initial_threshold=0,
    invert=False,
    canny_threshold=0.05, 
    canny_sigma=0, 
    canny_lt=0.05,  
    connected_pixels=200,  
    edge_length=50, 
    edge_buffer=100,
    max_buckets=255,
    min_bucket_width=0.001,
    max_raw=1e6,
    return_threshold=False,
):
    """Implementation of the Edge Otsu thresholding algorithm.
    Detailed explanation of algorithm can be found at https://doi.org/10.3390/rs12152469

    args:
        img (ee.Image): input image to thresholding algorithm
        band (str | None,optional): band name to use for thresholding, if set to `None` will use first band in image. default = None
        region (ee.Geometry | None, optional): region to determine threshold, if set to `None` will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90
        initial_threshold (float, optional): initial estimate of water/no-water for estimating the edges. default = 0
        invert (bool, optional): boolean switch to determine if to threshold greater than (True) or less than (False). default = False
        canny_threshold (float, optional): threshold for canny edge detection. default = 0.05
        canny_sigma (float, optional): sigma value for gaussian filter in canny edge detection. default = 0
        canny_lt (float, optional): lower threshold for canny detection. default = 0.05
        connected_pixels (int, optional): maximum size of the neighborhood in pixels to determine if connected. default = 200
        edge_length (int, optional): minimum length of edges from canny detection to be considered edge. default = 50
        edge_buffer (int, optional): distance in meters to buffer edges on a side for histogram sampling. default = 100
        max_buckets (int, optional): The maximum number of buckets to use when building a histogram; will be rounded up to a power of 2. default = 255
        min_bucket_width (float, optional): The minimum histogram bucket width to allow any power of 2. default = 0.001
        max_raw (int, optional): The number of values to accumulate before building the initial histogram. default = 1e6
        return_threshold (bool, optional): boolean switch, if set to true then function will return threshold number, else thresholded image. default = False

    returns:
        ee.Image: thresholded image based (if return_threshold==False) or threshold value (if return_threshold==True) based on the threshold determined by the algorithm
    """

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
        scale,
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
    """Otsu's method threhsolding algorithm.
    Computes single intensity threshold that separate histogram into two classes, foreground and background

    args:
        histogram (ee.Dictionary): computed object from ee.Reducer.histogram with keys "histogram" and "bucketMeans"

    returns:
        ee.Number: value of maximum inter-class intensity variance based on histogram
    """
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
    return ee.Number(output)


def kmeans_extent(img, hand, initial_threshold=0, region=None, band=None, samples=500, seed=7, invert=False, scale=90):
    """Water thresholding methodology using image values and HAND.
    Method from https://doi.org/10.1016/j.rse.2020.111732

    args:
        img (ee.Image): input image to thresholding algorithm
        hand (ee.Image): Height Above Nearest Drainage image used as axis in clustering
        initial_threshold (float, optional): initial estimate of water/no-water for stratified sampling. default = 0
        region (ee.Geometry | None, optional): region to sample values for KMeans clustering, if set to `None` will use img.geometry(). default = None
        band (str | None,optional): band name to use for thresholding, if set to `None` will use first band in image. default = None
        samples (int, optional): number of stratified samples to gather for clustering from the initial water/no-water map. default=500
        seed (int, optional): random number generator seed for sampling. default = 7
        invert (bool, optional): boolean switch to determine if class 1 is greater than initial_threshold then water (True),
             or less than initial_water then water (False). default = False
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90

    returns:
        ee.Image: clustered image from KMeans clusterer. Classes assumed to be water/no-water

    """
    if region is None:
        region = img.geometry()

    if band is None:
        img = img.select([0])
        band = ee.String(img.bandNames().get(0))

    else:
        img = img.select(band)

    hand_band = ee.String(hand.bandNames().get(0))

    # convert invert to a number 0 or 1 for ee server-side 
    invert = 1 if invert else 0

    strata = img.gt(initial_threshold).rename("strata")

    samples = ee.Image.cat([img, hand, strata]).stratifiedSample(
        numPoints=samples,
        classBand="strata",
        region=region,
        scale=scale,
        classValues=[0, 1],
        classPoints=[samples, samples],
        tileScale=16,
        seed=seed,
    )

    clusterer = ee.Clusterer.wekaKMeans(2, 2).train(samples, [band, hand_band])

    test = samples.cluster(clusterer,"classes")
    classes = ee.Array(test.aggregate_array("classes"))
    vals = ee.Array(test.aggregate_array(band)).mask(classes)

    class_mean = ee.Number(vals.reduce(ee.Reducer.mean(),[0]).get([0]))

    water = ee.Image.cat([img, hand.unmask(0)]).cluster(clusterer)

    do_inversion = class_mean.gt(initial_threshold).And(ee.Number(invert))

    water = ee.Image(ee.Algorithms.If(do_inversion, water.Not(), water))

    return water.rename("water").uint8()

