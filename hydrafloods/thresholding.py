import ee
from ee.ee_exception import EEException
import random
from hydrafloods import geeutils


def bmaxOtsu(collection,target_date,region,
             smoothing=100,
             qualityBand=None,
             reductionScale=90,
             initThresh=0,
             reverse=False,
             gridSize=0.1,
             bmaxThresh=0.75,
             maxBoxes=100,
             seed=7):

    def constuctGrid(i):
        def contructXGrid(j):
            j = ee.Number(j)
            box = ee.Feature(ee.Geometry.Rectangle(j,i,j.add(gridSize),i.add(gridSize)))
            out = ee.Algorithms.If(geom.contains(box.geometry()),box,None)
            return ee.Feature(out)
        i = ee.Number(i)
        out = ee.List.sequence(west,east.subtract(gridSize),gridSize).map(contructXGrid)
        return out


    def calcBmax(feature):
        segment = target.clip(feature)
        initial = segment.lt(initThresh)
        p1 = ee.Number(initial.reduceRegion(
            reducer= ee.Reducer.mean(),
            geometry= feature.geometry(),
            bestEffort= True,
            scale= reductionScale,
        ).get(histBand))
        p1 = ee.Number(ee.Algorithms.If(p1,p1,0.99))
        p2 = ee.Number(1).subtract(p1)

        m = segment.updateMask(initial).rename('m1').addBands(
            segment.updateMask(initial.Not()).rename('m2')
        )

        mReduced = m.reduceRegion(
            reducer= ee.Reducer.mean(),
            geometry= feature.geometry(),
            bestEffort= True,
            scale= reductionScale,
        )

        m1 = ee.Number(mReduced.get('m1'))
        m2 = ee.Number(mReduced.get('m2'))

        m1 = ee.Number(ee.Algorithms.If(m1,m1,globalLow))
        m2 = ee.Number(ee.Algorithms.If(m2,m2,globalHigh))

        sigmab = p1.multiply(p2.multiply(m1.subtract(m2).pow(2)))
        sigmat = ee.Number(segment.reduceRegion(
            reducer= ee.Reducer.variance(),
            geometry= feature.geometry(),
            bestEffort= True,
            scale= reductionScale,
        ).get(histBand))
        sigmat = ee.Number(ee.Algorithms.If(sigmat,sigmat,2))
        bmax = sigmab.divide(sigmat)
        return feature.set({'bmax':bmax})


    tDate = ee.Date(target_date)
    targetColl = collection.filterDate(tDate,tDate.advance(1,'day'))

    if qualityBand == None:
        target = targetColl.mosaic().select([0])
        histBand = ee.String(target.bandNames().get(0))

    else:
        histBand = ee.String(qualityBand)
        target = targetColl.qualityMosaic(qualityBand)\
            .select(histBand)

    geom = targetColl.map(geeutils.getGeoms).union().geometry()\
        .intersection(region,1)

    theoretical = target.reduceRegion(
        reducer= ee.Reducer.percentile([10,90]),
        geometry= geom,
        bestEffort= True,
        scale= 5000
    )
    globalLow = theoretical.get(histBand.cat('_p10'))
    globalHigh = theoretical.get(histBand.cat('_p90'))

    bounds = geom.bounds()
    coords = ee.List(bounds.coordinates().get(0))
    gridSize = ee.Number(gridSize)

    west = ee.Number(ee.List(coords.get(0)).get(0))
    south = ee.Number(ee.List(coords.get(0)).get(1))
    east = ee.Number(ee.List(coords.get(2)).get(0))
    north = ee.Number(ee.List(coords.get(2)).get(1))

    west = west.subtract(west.mod(gridSize))
    south = south.subtract(south.mod(gridSize))
    east = east.add(gridSize.subtract(east.mod(gridSize)))
    north = north.add(gridSize.subtract(north.mod(gridSize)))

    grid = ee.FeatureCollection(
      ee.List.sequence(south,north.subtract(gridSize),gridSize).map(constuctGrid).flatten()
    )

    bmaxes = grid.map(calcBmax).filter(ee.Filter.gt('bmax',bmaxThresh)).randomColumn('random',seed)

    nBoxes = ee.Number(bmaxes.size())
    randomThresh = ee.Number(maxBoxes).divide(nBoxes)
    selection = bmaxes.filter(ee.Filter.lt('random',randomThresh))

    histogram =  target.reduceRegion(ee.Reducer.histogram(255, 1)\
                                .combine('mean', None, True)\
                                .combine('variance', None,True),selection,reductionScale,bestEffort=True,
                                tileScale=16)

    threshold = otsu(histogram.get(histBand.cat('_histogram')))

    water = target.gt(threshold).clip(geeutils.LAND.geometry())

    return water.rename('water')

def globalOtsu(collection,target_date,region,
               canny_threshold=0.05, # threshold for canny edge detection
               canny_sigma=0,        # sigma value for gaussian filter
               canny_lt=0.05,        # lower threshold for canny detection
               smoothing=100,        # amount of smoothing in meters
               connected_pixels=200, # maximum size of the neighborhood in pixels
               edge_length=50,       # minimum length of edges from canny detection
               smooth_edges=100,
               qualityBand=None,
               reductionScale=90,
               initThresh=0,
               reverse=False,
               seed=7):


    tDate = ee.Date(target_date)
    targetColl = collection.filterDate(tDate,tDate.advance(1,'day'))

    if qualityBand == None:
        histBand = ee.String(target.bandNames().get(0))
        target = targetColl.mosaic()\
            .select(histBand)
    else:
        histBand = ee.String(qualityBand)
        target = targetColl.qualityMosaic(qualityBand)\
            .select(histBand)

    canny = ee.Algorithms.CannyEdgeDetector(target,canny_threshold,canny_sigma)

    connected = canny.mask(canny).lt(canny_lt).connectedPixelCount(connected_pixels, True)
    edges = connected.gte(edge_length)

    edgeBuffer = edges.focal_max(smooth_edges, 'square', 'meters')

    binary = edgeBuffer.rename('binary')

    samps = binary.stratifiedSample(numPoints=20,
        classBand='binary',
        region=region,
        scale=reductionScale,
        geometries=True,
        seed=seed,
        tileScale=16
    )
    sampleRegion = samps.geometry().buffer(2500)

    histogram_image = target.mask(edgeBuffer)

    histogram =  histogram_image.reduceRegion(ee.Reducer.histogram(255, 2)\
                                .combine('mean', None, True)\
                                .combine('variance', None,True),sampleRegion,reductionScale,bestEffort=True,
                                tileScale=16)

    threshold = otsu(histogram.get(histBand.cat('_histogram')))

    water = target.gt(threshold).clip(geeutils.LAND.geometry())

    return water.rename('water')

def bootstrapOtsu(collection,target_date, reductionPolygons,
                  neg_buffer=-1500,     # negative buffer for masking potential bad data
                  upper_threshold=-14,  # upper limit for water threshold
                  canny_threshold=7,    # threshold for canny edge detection
                  canny_sigma=1,        # sigma value for gaussian filter
                  canny_lt=7,           # lower threshold for canny detection
                  smoothing=100,        # amount of smoothing in meters
                  connected_pixels=200, # maximum size of the neighborhood in pixels
                  edge_length=50,       # minimum length of edges from canny detection
                  smooth_edges=100,
                  qualityBand=None,
                  reverse=False,
                  reductionScale=90):

    tDate = ee.Date(target_date)
    targetColl = collection.filterDate(tDate,tDate.advance(1,'day'))

    nImgs = targetColl.size().getInfo()
    if nImgs <= 0:
        raise EEException('Selected date has no imagery, please try processing another date')

    collGeom = targetColl.geometry()
    polygons = reductionPolygons.filterBounds(collGeom)

    nPolys = polygons.size().getInfo()
    if nPolys > 0:
        ids = ee.List(polygons.aggregate_array('id'))
        random_ids = []
        for i in range(3):
            random_ids.append(random.randint(0, ids.size().subtract(1).getInfo()))
        random_ids = ee.List(random_ids)

        def getRandomIds(i):
            return ids.get(i)

        ids = random_ids.map(getRandomIds)
        polygons = polygons.filter(ee.Filter.inList('id', ids))

        if qualityBand == None:
            target   = targetColl.mosaic().set('system:footprint', collGeom.dissolve())
            target   = target.clip(target.geometry().buffer(neg_buffer))
            smoothed = target.focal_median(smoothing, 'circle', 'meters')
            histBand = ee.String(target.bandNames().get(0))
        else:
            target   = targetColl.qualityMosaic(qualityBand).set('system:footprint', collGeom.dissolve())
            target   = target.clip(target.geometry().buffer(neg_buffer))
            smoothed = target.focal_median(smoothing, 'circle', 'meters')
            histBand = ee.String(qualityBand)

        canny = ee.Algorithms.CannyEdgeDetector(smoothed,canny_threshold,canny_sigma)

        connected = canny.mask(canny).lt(canny_lt).connectedPixelCount(connected_pixels, True)
        edges = connected.gte(edge_length)

        edgeBuffer = edges.focal_max(smooth_edges, 'square', 'meters')

        histogram_image = smoothed.mask(edgeBuffer)
        histogram = histogram_image.reduceRegion(ee.Reducer.histogram(255, 2),polygons.geometry(),reductionScale,bestEffort=True)

        threshold = ee.Number(otsu_function(histogram.get(histBand))).min(upper_threshold)
    else:
        threshold = upper_threshold

    water = smoothed.lt(threshold).clip(geeutils.LAND.geometry())

    return water

def otsu(histogram):
    counts = ee.Array(ee.Dictionary(histogram).get('histogram'))
    means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sums = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean = sums.divide(total)
    indices = ee.List.sequence(1, size)
    #Compute between sum of squares, where each mean partitions the data.

    def bss_function(i):
        aCounts = counts.slice(0, 0, i)
        aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
        aMeans = means.slice(0, 0, i)
        aMean = aMeans.multiply(aCounts).reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount)
        bCount = total.subtract(aCount)
        bMean = sums.subtract(aCount.multiply(aMean)).divide(bCount)
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(
               bCount.multiply(bMean.subtract(mean).pow(2)))

    bss = indices.map(bss_function)
    output = means.sort(bss).get([-1])
    return output
