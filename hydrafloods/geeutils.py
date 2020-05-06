import ee
from ee.ee_exception import EEException
import math
import string
import random
import datetime
from hydrafloods import decorators

try:
    ee.Initialize()
except EEException as e:
    from oauth2client.service_account import ServiceAccountCredentials
    credentials = ServiceAccountCredentials.from_p12_keyfile(
    service_account_email='',
    filename='',
    )
    ee.Initialize(credentials)

INITIME = ee.Date('1971-01-01T00:00:00')
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


def addTimeBand(img):
    t = ee.Date(img.get('system:time_start'))
    nDays = t.difference(INITIME, 'day')
    time = ee.Image(nDays).int16().rename('time')
    return img.addBands(time)


def getTileLayerUrl(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    tile_url_template = "https://earthengine.googleapis.com/map/{mapid}/{{z}}/{{x}}/{{y}}?token={token}"
    return tile_url_template.format(**map_id)


def exportImage(image, region, assetId, description=None, scale=1000, crs='EPSG:4326', pyramiding=None):
    if (description == None) or (type(description) != str):
        description = ''.join(random.SystemRandom().choice(
            string.ascii_letters) for _ in range(8)).lower()
    # get serializable geometry for export
    exportRegion = region.bounds(maxError=1).getInfo()['coordinates']

    if pyramiding is None:
        pyramiding = {'.default': 'mean'}

    # set export process
    export = ee.batch.Export.image.toAsset(image,
                                           description=description,
                                           assetId=assetId,
                                           scale=scale,
                                           region=exportRegion,
                                           maxPixels=1e13,
                                           crs=crs,
                                           pyramidingPolicy=pyramiding
                                           )
    # start export process
    export.start()

    return


def batchExport(collection, collectionAsset,  region=None, prefix=None, suffix=None, scale=1000, crs='EPSG:4326', metadata=None, pyramiding=None,verbose=False):
    if type(collection) is not ee.imagecollection.ImageCollection:
        try:
            collection = getattr(collection,'collection')
        except Exception as e:
            raise TypeError('argument collection needs to be either of type ee.ImageCollection '
                            'or hydrafloods.hfCollection')

    n = collection.size()
    exportImages = collection.sort('system:time_start', False).toList(n)
    nIter = n.getInfo()

    for i in range(nIter):
        img = ee.Image(exportImages.get(i))
        if metadata is not None:
            img = img.set(metadata)

        t = img.get('system:time_start').getInfo()
        date = datetime.datetime.utcfromtimestamp(
            t / 1e3).strftime("%Y%m%d")

        if region is None:
            region = img.geometry()

        exportName = date
        if prefix is not None:
            exportName = f"{prefix}_" + exportName
        if suffix is not None:
            exportName = exportName + f"_{suffix}"

        description = exportName
        if verbose:
            print(f"running export for {description}")

        if not collectionAsset.endswith('/'):
            collectionAsset += '/'

        exportName=collectionAsset + description

        exportImage(img, region, exportName, description=description,
                    scale=scale, crs=crs, pyramiding=pyramiding)

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


@decorators.carryMetadata
def powerToDb(img):
    return ee.Image(10).multiply(img.log10())

@decorators.carryMetadata
def dbToPower(img):
    return ee.Image(10).pow(img.divide(10))

@decorators.carryMetadata
def addIndices(img):
    ndvi = img.normalizedDifference(['nir','red']).rename('ndvi')
    mndwi = img.normalizedDifference(['green','swir2']).rename('mndwi')
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

    indices = ee.Image.cat([ndvi,mndwi,nwi,aewinsh,aewish,tcwet])

    return img.addBands(indices)
