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

LAND = ee.FeatureCollection('USDOS/LSIB/2013')
#S1_polygons = ee.FeatureCollection('projects/servir-mekong/hydrafloods/S1_polygons')

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


def getGeoms(img):
    return img.geometry()


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
