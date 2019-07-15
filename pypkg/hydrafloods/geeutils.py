import ee
from ee.ee_exception import EEException
import math
import string
import random
import datetime

try:
    ee.Initialize()
except EEException as e:
    from oauth2client.service_account import ServiceAccountCredentials
    credentials = ServiceAccountCredentials.from_p12_keyfile(
    service_account_email='',
    filename='',
    )
    ee.Initialize(credentials)

landShp = ee.FeatureCollection('USDOS/LSIB/2013')
S1_polygons = ee.FeatureCollection('projects/servir-mekong/hydrafloods/S1_polygons')

# helper function to convert qa bit image to flag
def extractBits(image, start, end, newName):
    # Compute the bits we need to extract.
    pattern = 0
    for i in range(start,end):
       pattern += int(math.pow(2, i))

    # Return a single band image of the extracted QA bits, giving the band
    # a new name.
    return image.select([0], [newName])\
                  .bitwiseAnd(pattern)\
                  .rightShift(start)


def getTileLayerUrl(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    tile_url_template = "https://earthengine.googleapis.com/map/{mapid}/{{z}}/{{x}}/{{y}}?token={token}"
    return tile_url_template.format(**map_id)


def exportImage(image,region,assetId,description=None,scale=90,crs='EPSG:4326'):
    if (description == None) or (type(description) != str):
        description = ''.join(random.SystemRandom().choice(
        string.ascii_letters) for _ in range(8)).lower()
    # get serializable geometry for export
    exportRegion = region.bounds().getInfo()['coordinates']

    # set export process
    export = ee.batch.Export.image.toAsset(image,
      description = description,
      assetId = assetId,
      scale = scale,
      region = exportRegion,
      maxPixels = 1e13,
      crs = crs
    )
    # start export process
    export.start()

    return

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


def rescaleBands(img):
    def individualBand(b):
        b = ee.String(b)
        minKey = b.cat('_min')
        maxKey = b.cat('_max')
        return img.select(b).unitScale(
            ee.Number(minMax.get(minKey)),
            ee.Number(minMax.get(maxKey)))

    bandNames = img.bandNames()
    geom = img.geometry()

    minMax = img.reduceRegion(
        reducer = ee.Reducer.minMax(),
        geometry = geom,
        scale = 90,
        bestEffort = True
    )

    rescaled = bandNames.map(individualBand)
    bandSequence = ee.List.sequence(0,bandNames.length().subtract(1))

    return ee.ImageCollection(rescaled).toBands().select(bandSequence,bandNames)


def logitTransform(img):
    return img.divide(ee.Image(1).subtract(img)).log()


def toNatural(img):
    return ee.Image(10.0).pow(img.select(0).divide(10.0))


def toDB(img):
    return ee.Image(img).log10().multiply(10.0)


def despeckle(img):
  """ Refined Lee Speckle Filter """
  t = ee.Date(img.get('system:time_start'))
  # angles = img.select('angle')
  geom = img.geometry()
  # The RL speckle filter
  img = toNatural(img)
  # img must be in natural units, i.e. not in dB!
  # Set up 3x3 kernels
  weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
  kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

  mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
  variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

  # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0],\
                            [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

  sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False)

  # Calculate mean and variance for the sampled windows and store as 9 bands
  sample_mean = mean3.neighborhoodToBands(sample_kernel)
  sample_var = variance3.neighborhoodToBands(sample_kernel)

  # Determine the 4 gradients for the sampled windows
  gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

  # And find the maximum gradient amongst gradient bands
  max_gradient = gradients.reduce(ee.Reducer.max())

  # Create a mask for band pixels that are the maximum gradient
  gradmask = gradients.eq(max_gradient)

  # duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask)

  # Determine the 8 directions
  directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).\
                                  subtract(sample_mean.select(7))).multiply(1)

  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).\
                                   gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))

  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).\
                                   gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))

  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).\
                                   gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
  # The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).Not().multiply(5))
  directions = directions.addBands(directions.select(1).Not().multiply(6))
  directions = directions.addBands(directions.select(2).Not().multiply(7))
  directions = directions.addBands(directions.select(3).Not().multiply(8))

  # Mask all values that are not 1-8
  directions = directions.updateMask(gradmask)

  # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum())

  sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

  # Calculate localNoiseVariance
  sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

  # Set up the 7*7 kernels for directional statistics
  rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

  diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],\
                          [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]])

  rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
  diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

  # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
  dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

  # and add the bands for rotated kernels
  i = 1
  while i < 4:
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
    i+=1
  # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum())
  dir_var = dir_var.reduce(ee.Reducer.sum())

  # A finally generate the filtered value
  varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))

  b = varX.divide(dir_var)

  # Get a multi-band image bands.
  result = dir_mean.add(b.multiply(img.subtract(dir_mean))).arrayProject([0])\
    .arrayFlatten([['sum']])\
    .float()


  return toDB(result)

def addIndices(img):
    ndvi = img.normalizedDifference(['nir','red']).rename('ndvi')
    mndwi = img.normalizedDifference(['green','swir1']).rename('mndwi')
    nwi = img.expression('((b-(n+s+w))/(b+(n+s+w))*100)',{
        'b':img.select('blue'),
        'n':img.select('nir'),
        's':img.select('swir1'),
        'w':img.select('swir2')
    }).rename('nwi')
    aewinsh = img.expression('4.0 * (g-s) - ((0.25*n) + (2.75*w))',{
        'g':img.select('green'),
        's':img.select('swir1'),
        'n':img.select('nir'),
        'w':img.select('swir2')
    }).rename('aewinsh')
    aewish = img.expression('b+2.5*g-1.5*(n+s)-0.25*w',{
        'b':img.select('blue'),
        'g':img.select('green'),
        'n':img.select('nir'),
        's':img.select('swir1'),
        'w':img.select('swir2')
    }).rename('aewish')
    tcwet = img.expression('0.1509*b + 0.1973*g + 0.3279*r + 0.3406*n - 0.7112*s - 0.4572*w',{
        'b':img.select('blue'),
        'g':img.select('green'),
        'r':img.select('red'),
        'n':img.select('nir'),
        's':img.select('swir1'),
        'w':img.select('swir2')
    }).rename('tcwet')

    return ee.Image.cat([img,ndvi,mndwi,nwi,aewinsh,aewish,tcwet])


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
