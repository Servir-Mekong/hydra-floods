import ee
from ee.ee_exception import EEException
import random
from . import geeutils


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

    threshold = geeutils.otsu_function(histogram.get(histBand.cat('_histogram')))

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

    print(canny_threshold,canny_sigma,canny_lt,qualityBand)

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

    threshold = geeutils.otsu_function(histogram.get(histBand.cat('_histogram')))

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


def SurfaceWaterAlgorithm(aoi,images, pcnt_perm, pcnt_temp, water_thresh, ndvi_thresh, hand_mask):

    STD_NAMES   = ['blue2', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2']

    # calculate percentile images
    prcnt_img_perm = images.reduce(ee.Reducer.percentile([float(pcnt_perm)])).rename(STD_NAMES)
    prcnt_img_temp = images.reduce(ee.Reducer.percentile([float(pcnt_temp)])).rename(STD_NAMES)

    # MNDWI
    MNDWI_perm = prcnt_img_perm.normalizedDifference(['green', 'swir1'])
    MNDWI_temp = prcnt_img_temp.normalizedDifference(['green', 'swir1'])

    # water
    water_perm = MNDWI_perm.gt(float(water_thresh))
    water_temp = MNDWI_temp.gt(float(water_thresh))

    # get NDVI masks
    NDVI_perm_pcnt = prcnt_img_perm.normalizedDifference(['nir', 'red'])
    NDVI_temp_pcnt = prcnt_img_temp.normalizedDifference(['nir', 'red'])
    NDVI_mask_perm = NDVI_perm_pcnt.gt(float(ndvi_thresh))
    NDVI_mask_temp = NDVI_temp_pcnt.gt(float(ndvi_thresh))

    # combined NDVI and HAND masks
    full_mask_perm = NDVI_mask_perm.add(hand_mask)
    full_mask_temp = NDVI_mask_temp.add(hand_mask)

    # apply NDVI and HAND masks
    water_perm_masked = water_perm.updateMask(full_mask_perm.Not())
    water_temp_masked = water_temp.updateMask(full_mask_perm.Not())

    # single image with permanent and temporary water
    water_complete = water_perm_masked.add(water_temp_masked).clip(aoi)

    #return water_complete.updateMask(water_complete)
    return water_complete

def getLandsatCollection(aoi,time_start, time_end, month_index=None, climatology=True, defringe=True, cloud_thresh=None):

    # Landsat band names
    LC457_BANDS = ['B1',    'B1',   'B2',    'B3',  'B4',  'B5',    'B7']
    LC8_BANDS   = ['B1',    'B2',   'B3',    'B4',  'B5',  'B6',    'B7']
    STD_NAMES   = ['blue2', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2']


    # filter Landsat collections on bounds and dates
    L4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA').filterBounds(aoi).filterDate(time_start, ee.Date(time_end).advance(1, 'day'))
    L5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').filterBounds(aoi).filterDate(time_start, ee.Date(time_end).advance(1, 'day'))
    L7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterBounds(aoi).filterDate(time_start, ee.Date(time_end).advance(1, 'day'))
    L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterBounds(aoi).filterDate(time_start, ee.Date(time_end).advance(1, 'day'))

    # apply cloud masking
    if int(cloud_thresh) >= 0:
        # helper function: cloud busting
        # (https://code.earthengine.google.com/63f075a9e212f6ed4770af44be18a4fe, Ian Housman and Carson Stam)
        def bustClouds(img):
            t = img
            cs = ee.Algorithms.Landsat.simpleCloudScore(img).select('cloud')
            out = img.mask(img.mask().And(cs.lt(ee.Number(int(cloud_thresh)))))
            return out.copyProperties(t)
        # apply cloud busting function
        L4 = L4.map(bustClouds)
        L5 = L5.map(bustClouds)
        L7 = L7.map(bustClouds)
        L8 = L8.map(bustClouds)

    # select bands and rename
    L4 = L4.select(LC457_BANDS, STD_NAMES)
    L5 = L5.select(LC457_BANDS, STD_NAMES)
    L7 = L7.select(LC457_BANDS, STD_NAMES)
    L8 = L8.select(LC8_BANDS, STD_NAMES)

    # apply defringing
    if defringe == 'true':
        # helper function: defringe Landsat 5 and/or 7
        # (https://code.earthengine.google.com/63f075a9e212f6ed4770af44be18a4fe, Bonnie Ruefenacht)
        k = ee.Kernel.fixed(41, 41, \
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        fringeCountThreshold = 279  # define number of non null observations for pixel to not be classified as a fringe
        def defringeLandsat(img):
            m   = img.mask().reduce(ee.Reducer.min())
            sum = m.reduceNeighborhood(ee.Reducer.sum(), k, 'kernel')
            sum = sum.gte(fringeCountThreshold)
            img = img.mask(img.mask().And(sum))
            return img
        L5 = L5.map(defringeLandsat)
        L7 = L7.map(defringeLandsat)

    # merge collections
    images = ee.ImageCollection(L4.merge(L5).merge(L7).merge(L8))

    # filter on selected month
    if climatology:
        if month_index != None:
            images = images.filter(ee.Filter.calendarRange(int(month_index), int(month_index), 'month'))
        else:
            raise ValueError('Month needs to be defined to calculate climatology')

    return images

def JRCAlgorithm(geom,startDate, endDate, month=None):

    IMAGE_COLLECTION = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory')

    myjrc = IMAGE_COLLECTION.filterBounds(geom).filterDate(startDate, endDate)

    if month != None:
        myjrc = myjrc.filter(ee.Filter.calendarRange(int(month), int(month), 'month'))

    # calculate total number of observations
    def calcObs(img):
        # observation is img > 0
        obs = img.gt(0)
        return ee.Image(obs).set('system:time_start', img.get('system:time_start'))

    # calculate the number of times water
    def calcWater(img):
        water = img.select('water').eq(2)
        return ee.Image(water).set('system:time_start', img.get('system:time_start'))

    observations = myjrc.map(calcObs)

    water = myjrc.map(calcWater)

    # sum the totals
    totalObs = ee.Image(ee.ImageCollection(observations).sum().toFloat())
    totalWater = ee.Image(ee.ImageCollection(water).sum().toFloat())

    # calculate the percentage of total water
    returnTime = totalWater.divide(totalObs).multiply(100)

    # make a mask
    water = returnTime.gt(75).rename(['water'])
    water = water.updateMask(water)

    return water


def getHistoricalMap(geom,iniTime,endTime,
                  climatology=True,
                  month=None,
                  defringe=True,
                  pcnt_perm=40,
                  pcnt_temp=8,
                  water_thresh=0.35,
                  ndvi_thresh=0.5,
                  hand_thresh=30,
                  cloud_thresh=10,
                  algorithm='SWT'):

    def spatialSelect(feature):
        test = ee.Algorithms.If(geom.contains(feature.geometry()),feature,None)
        return ee.Feature(test)

    countries = geeutils.LAND.filterBounds(geom).map(spatialSelect,True)

    if climatology:
        if month == None:
            raise ValueError('Month needs to be defined to calculate climatology')

    if algorithm == 'SWT':
        # get images
        images = getLandsatCollection(geom,iniTime, endTime, climatology, month, defringe, cloud_thresh)

        # Height Above Nearest Drainage (HAND)
        HAND = ee.Image('users/arjenhaag/SERVIR-Mekong/HAND_MERIT')

        # get HAND mask
        HAND_mask = HAND.gt(float(hand_thresh))


        water = SurfaceWaterAlgorithm(geom,images, pcnt_perm, pcnt_temp, water_thresh, ndvi_thresh, HAND_mask).clip(countries)
        waterMap = getTileLayerUrl(water.updateMask(water.eq(2)).visualize(min=0,max=2,palette='#ffffff,#9999ff,#00008b'))

    elif algorithm == 'JRC':
        water = JRCAlgorithm(geom,iniTime,endTime).clip(countries)
        waterMap = getTileLayerUrl(water.visualize(min=0,max=1,bands='water',palette='#ffffff,#00008b'))

    else:
        raise NotImplementedError('Selected algorithm string not available. Options are: "SWT" or "JRC"')

    return waterMap
