import ee
from ee.ee_exception import EEException
import random
from hydrafloods import geeutils, decorators, ml, fuzzy


@decorators.carry_metadata
def bmax_otsu(
    img,
    band=None,
    region=None,
    scale=90,
    initial_threshold=0,
    thresh_no_data=None,
    invert=False,
    grid_size=0.1,
    bmax_threshold=0.75,
    iters=1,
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
        thresh_no_data (float, optional): threshold to be used when histogram is empty (likely little to no water present in image). default = -0.2
        invert (bool, optional): boolean switch to determine if to threshold greater than (True) or less than (False). default = False
        grid_size (float, optional): size in decimal degrees to tile image/region to check for bimodality. default = 0.1
        bmax_threshold (float, optional): value 0-1 to determine if a value of bmax is bimodal or not. default = 0.75
        iters (int, optional): number of iterations to successively shrink the bmax tiles to refine threshold. 1 uses defined grid_size. default = 1
        max_boxes (int, optional): maximum number of tiles/boxes to use when determining threshold. default = None
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
        p1 = ee.Number(ee.Algorithms.If(p1, p1, 0.9999))
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

    grid = geeutils.tile_region(region, centroid_within=region, grid_size=grid_size)

    bmaxes = grid.map(calcBmax).filter(ee.Filter.gt("bmax", bmax_threshold))

    if iters > 1:
        for i in range(iters - 1):
            grid_size /= 2.0
            grid = geeutils.tile_region(
                bmaxes.geometry(), centroid_within=bmaxes, grid_size=grid_size
            )

            bmaxes = grid.map(calcBmax).filter(ee.Filter.gt("bmax", bmax_threshold))

    if max_boxes is not None:
        selection = bmaxes.randomColumn("random", seed).limit(max_boxes, "random")
    else:
        selection = bmaxes

    histogram = img.reduceRegion(
        ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
        .combine("mean", None, True)
        .combine("variance", None, True),
        selection,
        scale,
        bestEffort=True,
        tileScale=16,
    )

    if thresh_no_data is not None:
        threshold = ee.Number(
            ee.Algorithms.If(
                ee.Dictionary(histogram.get(histBand.cat("_histogram"))).contains(
                    "bucketMeans"
                ),
                otsu(histogram.get(histBand.cat("_histogram"))),
                thresh_no_data,
            )
        )
    else:
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
    thresh_no_data=None,
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
        thresh_no_data (float, optional): threshold to be used when histogram is empty (likely little to no water present in image). default = -0.2
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

    if thresh_no_data is not None:
        threshold = ee.Number(
            ee.Algorithms.If(
                ee.Dictionary(histogram.get(histBand.cat("_histogram"))).contains(
                    "bucketMeans"
                ),
                otsu(histogram.get(histBand.cat("_histogram"))),
                thresh_no_data,
            )
        )
    else:
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


def kmeans_extent(
    img,
    hand,
    initial_threshold=0,
    region=None,
    band=None,
    n_samples=500,
    seed=7,
    invert=False,
    scale=90,
):
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
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90

    returns:
        ee.Image: clustered image from KMeans clusterer. Classes assumed to be water/no-water

    """

    def _cluster_center(x):
        return ee.Number(
            (
                band_arr.mask(classes.eq(ee.Number(x))).reduce(ee.Reducer.mean(), [0])
            ).get([0])
        )

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

    img = img.addBands(hand.unmask(0))

    strata = img.select(band).gt(initial_threshold).rename("strata")

    samples = img.addBands(strata).stratifiedSample(
        numPoints=n_samples,
        classBand="strata",
        region=region,
        scale=scale,
        classValues=[0, 1],
        classPoints=[n_samples, n_samples],
        tileScale=16,
        seed=seed,
    )

    clusterer = ee.Clusterer.wekaKMeans(2, 1).train(samples, [band, hand_band])

    samples = samples.cluster(clusterer, "classes")
    classes = samples.aggregate_array("classes")
    unique = classes.distinct().sort()
    classes = ee.Array(classes)
    band_arr = ee.Array(samples.aggregate_array(band))

    class_means = unique.map(_cluster_center)

    min_mean = class_means.reduce(ee.Reducer.min())

    water_class = class_means.indexOf(min_mean)

    water = img.cluster(clusterer)

    water = ee.Image(ee.Algorithms.If(water_class.eq(0), water.Not(), water))

    return water.rename("water").uint8()


def multidim_semisupervised(
    img,
    bands,
    rank_band=None,
    ranking="min",
    region=None,
    n_samples=500,
    seed=7,
    scale=90,
    proba_threshold=None,
):
    """Implementation of the Automatic water detection from 
    multidimensional hierarchical clustering algorithm.
    Method details: https://doi.org/10.1016/j.rse.2020.112209
    Note: this is a similar method, not exact from the paper

    args:
        img (ee.Image): input image to thresholding algorithm
        bands (list | ee.List): band names to use for the semi supervised classification
        rank_band (str, optional): band name used to rank which unserpervised class is water. If None then first band name in `bands` is used. default = None
        ranking (str, optional): method to rank the classes by `rank_band`. Options are 'min' or 'max'. If 'min', then the lowest class mean is considered water. default = 'min'
        region (ee.Geometry | None, optional): region to sample values for KMeans clustering, if set to `None` will use img.geometry(). default = None
        samples (int, optional): number of stratified samples to gather for clustering from the initial water/no-water map. default=500
        seed (int, optional): random number generator seed for sampling. default = 7
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90
        proba_threshold (float, optional): probability threshold to create a binary water mask, should be in range of 0-1. If None, then the probability values are returned. default = None
    """

    if region is None:
        region = img.geometry()

    if bands is None:
        bands = img.bandNames()
    else:
        bands = ee.List(bands)

    if rank_band is None:
        rank_band = ee.String(bands.get(0))

    samples = img.select(bands).sample(
        region=region,
        scale=scale,
        numPixels=n_samples,
        seed=seed,
        dropNulls=True,
        tileScale=16,
    )

    classifier = ml.unsupervised_rf(
        100, samples, features=bands, rank_feature=rank_band, ranking=ranking
    )

    probas = img.select(bands).classify(classifier)

    if proba_threshold is None:
        output = probas.rename("water_proba")
    else:
        output = probas.gt(proba_threshold).rename("water").uint8()

    return output


@decorators.carry_metadata
def fuzzy_otsu(
    img,
    band=None,
    region=None,
    scale=90,
    initial_threshold=-15.5,
    grid_size=0.1,
    max_boxes=20,
    seed=0,
    max_buckets=255,
    min_bucket_width=0.001,
    max_raw=1e6,
):
    """ Implementation of Otsu thresholding algorithm.
    Segment the grids into water and land representation using initial threshold and randomly select features to calculate 
    minimum and maximum threshold of the image. Calculate Midpoint of min and max threshold using Fuzzy_Gaussian to map water.


    args:
        img (ee.Image): input image to thresholding algorithm
        band (str | None,optional): band name to use for thresholding, if set to `None` will use first band in image. default = None
        region (ee.Geometry | None, optional): region to determine threshold, if set to `None` will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90
        initial_threshold (float, optional): initial estimate of water/no-water for estimating the probabilities of classes in segment. default = -15.5
        grid_size (float, optional): size in decimal degrees to tile image/region to check for bimodality. default = 0.1
        max_boxes (int, optional): maximum number of tiles/boxes to use when determining threshold. default = 20
        seed (int, optional): random number generator seed for randomly selected max_boxes. default = 7
        max_buckets (int, optional): The maximum number of buckets to use when building a histogram; will be rounded up to a power of 2. default = 255
        min_bucket_width (float, optional): The minimum histogram bucket width to allow any power of 2. default = 0.001
        max_raw (int, optional): The number of values to accumulate before building the initial histogram. default = 1e6
        
    returns:
        ee.Image: thresholded water image.    
    """

    # Function to Segment Land and Water Grids
    def segment_grid(feature):
        counts = initImg.reduceRegion(
            reducer=ee.Reducer.sum(),
            geometry=feature.geometry(),
            bestEffort=True,
            scale=scale,
        )

        return feature.set(counts)

    if region is None:
        region = img.geometry()

    if band is None:
        img = img.select([0])
        histBand = ee.String(img.bandNames().get(0))

    else:
        histBand = ee.String(band)
        img = img.select(histBand)

    grid = geeutils.tile_region(region, intersect_geom=region, grid_size=grid_size)

    # Prepare good representation of water and land
    initWater = img.lt(initial_threshold).rename("water")
    initImg = initWater.addBands(initWater.Not().rename("land"))

    grid = grid.map(segment_grid)

    # Select water and land grid
    selection = grid.filter(
        ee.Filter.And(ee.Filter.gt("water", 1000), ee.Filter.gt("land", 1000))
    )

    # Randomly select features
    selection = selection.randomColumn("random")
    nBoxes = ee.Number(selection.size())
    randomSelection = ee.Number(max_boxes).divide(nBoxes)
    selection = selection.filter(ee.Filter.lt("random", randomSelection))

    # Create histogram for selected features to calculate min threshold
    histogram = img.reduceRegion(
        ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
        .combine("mean", None, True)
        .combine("variance", None, True),
        selection.geometry(),
        scale,
        bestEffort=True,
        tileScale=16,
    )

    min_threshold = otsu(histogram.get(histBand.cat("_histogram")))

    # Create histogram to calculate max threshold
    histogram2 = img.updateMask(img.lt(min_threshold)).reduceRegion(
        ee.Reducer.histogram(max_buckets, min_bucket_width, max_raw)
        .combine("mean", None, True)
        .combine("variance", None, True),
        selection.geometry(),
        scale,
        bestEffort=True,
        tileScale=16,
    )

    max_threshold = otsu(histogram2.get(histBand.cat("_histogram")))

    # Calculate Midpoint of min and max threshold using Fuzzy_Gaussian
    midpoint = ee.Number(ee.Number(min_threshold).add(ee.Number(max_threshold))).divide(2)
    spread = 0.2

    gauss = fuzzy.fuzzy_gaussian(img, midpoint, spread).clip(region)

    waterImg = img.lt(midpoint)

    waterImg = waterImg.where(waterImg.eq(0), gauss)

    return ee.Image(waterImg)

