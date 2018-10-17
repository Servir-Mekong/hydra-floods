import ee
import datetime

try:
    ee.Initialize()
except EEException as e:
    from oauth2client.service_account import ServiceAccountCredentials
    credentials = ee.ServiceAccountCredentials(
    service_account_email='',
    filename='',
    )
    ee.Initialize(credentials)

landShp = ee.FeatureCollection('USDOS/LSIB/2013')


def getTileLayerUrl(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    tile_url_template = "https://earthengine.googleapis.com/map/{mapid}/{{z}}/{{x}}/{{y}}?token={token}"
    return tile_url_template.format(**map_id)


def otsu_function(histogram):
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


def s1WaterMap(geom,iniTime,endTime,
               canny_threshold=7,    # threshold for canny edge detection
               canny_sigma=1,        # sigma value for gaussian filter
               canny_lt=7,           # lower threshold for canny detection
               smoothing=100,        # amount of smoothing in meters
               connected_pixels=200, # maximum size of the neighborhood in pixels
               edge_length=50,       # minimum length of edges from canny detection
               smooth_edges=100):

    def spatialSelect(feature):
        test = ee.Algorithms.If(geom.contains(feature.geometry()),feature,None)
        return ee.Feature(test)

    collection = ee.ImageCollection('COPERNICUS/S1_GRD')\
                 .filterBounds(geom)\
                 .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
                 .select('VV').filterDate(iniTime, endTime)

    after = collection.mosaic()

    smoothed = after.focal_median(smoothing, 'circle', 'meters')

    canny = ee.Algorithms.CannyEdgeDetector(smoothed,canny_threshold,canny_sigma)

    connected = canny.mask(canny).lt(canny_lt).connectedPixelCount(connected_pixels, True)
    edges = connected.gte(edge_length)

    edgeBuffer = edges.focal_max(smooth_edges, 'square', 'meters');

    imageEdge = smoothed.mask(edges)

    polygon = ee.Geometry.MultiPolygon(
        [[[[94.32174682617188, 22.916025590238387],
           [94.47761535644531, 22.915393135699805],
           [94.4879150390625, 23.02223600607765],
           [94.32723999023438, 23.026027690875438]]],
         [[[94.14939880371094, 23.187706723229915],
           [94.21119689941406, 23.18581317530292],
           [94.20913696289062, 23.252702265171013],
           [94.14665222167969, 23.250178757002324]]]])

    histogram_image = smoothed.mask(edgeBuffer)
    histogram = histogram_image.reduceRegion(ee.Reducer.histogram(255, 2)\
                                .combine('mean', None, True)\
                                .combine('variance', None,True),polygon,10,bestEffort=True)

    threshold = otsu_function(histogram.get('VV_histogram'));

    countries = landShp.filterBounds(geom).map(spatialSelect,True)

    water = smoothed.mask(smoothed.lt(threshold)).clip(countries)

    mapUrl = getTileLayerUrl(water.visualize(palette='#9999ff'))

    return mapUrl



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
        obs = img.gt(0);
        return ee.Image(obs).set('system:time_start', img.get('system:time_start'));

    # calculate the number of times water
    def calcWater(img):
        water = img.select('water').eq(2);
        return ee.Image(water).set('system:time_start', img.get('system:time_start'));

    observations = myjrc.map(calcObs)

    water = myjrc.map(calcWater)

    # sum the totals
    totalObs = ee.Image(ee.ImageCollection(observations).sum().toFloat());
    totalWater = ee.Image(ee.ImageCollection(water).sum().toFloat());

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

    countries = landShp.filterBounds(geom).map(spatialSelect,True)

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


def _accumulate(ic,today,ndays=1):

    eeDate = ee.Date(today)

    ic_filtered = ic.filterDate(eeDate.advance(-ndays,'day'),eeDate)

    accum_img = ee.Image(ic_filtered.sum())

    return accum_img.updateMask(accum_img.gt(1))

def getPrecipMap(accumulation=1):
    if accumulation not in [1,3,7]:
        raise NotImplementedError('Selected accumulation value is not yet impleted, options are: 1, 3, 7')

    dt = datetime.datetime.utcnow() - datetime.timedelta(1)
    today = dt.strftime('%Y-%m-%d')

    ic = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/operational').select(['hourlyPrecipRateGC'])

    ranges = {1:[1,100],3:[1,250],7:[1,500]}
    crange = ranges[accumulation]

    accum = _accumulate(ic,today,accumulation)

    precipMap = getTileLayerUrl(accum.visualize(min=crange[0],max=crange[1],
                                                palette='#000080,#0045ff,#00fbb2,#67d300,#d8ff22,#ffbe0c,#ff0039,#c95df5,#fef8fe'
                                               )
                               )
    return precipMap

def getAdminMap(geom):
    def spatialSelect(feature):
        test = ee.Algorithms.If(geom.contains(feature.geometry()),feature,None)
        return ee.Feature(test)

    countries = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').map(spatialSelect,True)

    # Create an empty image into which to paint the features, cast to byte.
    empty = ee.Image().byte()

    # Paint all the polygon edges with the same number and width, display.
    outline = empty.paint(
      featureCollection=countries,
      width=2
    )

    return getTileLayerUrl(outline.visualize())
