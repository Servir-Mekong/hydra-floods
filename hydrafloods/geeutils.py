import ee
from ee.ee_exception import EEException
import math
import string
import random
import datetime
from hydrafloods import decorators


# helper function to convert qa bit image to flag
def extract_bits(image, start, end=None, new_name=None):
    newame = new_name if new_name is not None else f'{start}Bits'

    if (start == end) or (end is None):
        #perform a bit shift with bitwiseAnd
        return image.select([0],[new_name])\
            .bitwiseAnd(1<<start)
    else:
        # Compute the bits we need to extract.
        pattern = 0
        for i in range(start, end):
            pattern += int(math.pow(2, i))

        # Return a single band image of the extracted QA bits, giving the band
        # a new name.
        return image.select([0], [new_name])\
            .bitwiseAnd(pattern)\
            .rightShift(start)


def get_geoms(img):
    return img.geometry()


def get_tile_layer_url(ee_image_object):
    map_id = ee.Image(ee_image_object).getMapId()
    tile_url_template = (
        "https://earthengine.googleapis.com/map/{mapid}/{{z}}/{{x}}/{{y}}?token={token}"
    )
    return tile_url_template.format(**map_id)


def export_image(
    image,
    region,
    asset_id,
    description=None,
    scale=1000,
    crs="EPSG:4326",
    pyramiding=None,
):
    if (description == None) or (type(description) != str):
        description = "".join(
            random.SystemRandom().choice(string.ascii_letters) for _ in range(8)
        ).lower()
    # get serializable geometry for export
    export_region = region.bounds(maxError=1).getInfo()["coordinates"]

    if pyramiding is None:
        pyramiding = {".default": "mean"}

    # set export process
    export = ee.batch.Export.image.toAsset(
        image,
        description=description,
        assetId=asset_id,
        scale=scale,
        region=export_region,
        maxPixels=1e13,
        crs=crs,
        pyramidingPolicy=pyramiding,
    )
    # start export process
    export.start()

    return


def export_table(
    table,
    region,
    assetId,
    description=None,
    scale=1000,
    crs="EPSG:4326",
    pyramiding=None,
):
    if (description == None) or (type(description) != str):
        description = "".join(
            random.SystemRandom().choice(string.ascii_letters) for _ in range(8)
        ).lower()
    # get serializable geometry for export
    exportRegion = region.bounds(maxError=1).getInfo()["coordinates"]

    if pyramiding is None:
        pyramiding = {".default": "mean"}

    # set export process
    export = ee.batch.Export.image.toAsset(
        image,
        description=description,
        assetId=assetId,
        scale=scale,
        region=exportRegion,
        maxPixels=1e13,
        crs=crs,
        pyramidingPolicy=pyramiding,
    )
    # start export process
    export.start()

    return


def batch_export(
    collection,
    collection_asset,
    region=None,
    prefix=None,
    suffix=None,
    scale=1000,
    crs="EPSG:4326",
    metadata=None,
    pyramiding=None,
    verbose=False,
):
    if type(collection) is not ee.imagecollection.ImageCollection:
        try:
            collection = getattr(collection, "collection")
        except Exception as e:
            raise TypeError(
                "argument collection needs to be either of type ee.ImageCollection "
                "or hydrafloods.hfCollection"
            )

    n = collection.size()
    exportImages = collection.sort("system:time_start", False).toList(n)
    nIter = n.getInfo()

    for i in range(nIter):
        img = ee.Image(exportImages.get(i))
        if metadata is not None:
            img = img.set(metadata)

        t = img.get("system:time_start").getInfo()
        date = datetime.datetime.utcfromtimestamp(t / 1e3).strftime("%Y%m%d")

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

        if not collectionAsset.endswith("/"):
            collectionAsset += "/"

        exportName = collection_asset + description

        exportImage(
            img,
            region,
            exportName,
            description=description,
            scale=scale,
            crs=crs,
            pyramiding=pyramiding,
        )

    return


@decorators.carry_metadata
def rescale_bands(img):
    def _individual_band(b):
        b = ee.String(b)
        minKey = b.cat("_min")
        maxKey = b.cat("_max")
        return img.select(b).unitScale(
            ee.Number(minMax.get(minKey)), ee.Number(minMax.get(maxKey))
        )

    bandNames = img.bandNames()
    geom = img.geometry()

    minMax = img.reduceRegion(
        reducer=ee.Reducer.minMax(), geometry=geom, scale=90, bestEffort=True
    )

    rescaled = bandNames.map(_individual_band)
    bandSequence = ee.List.sequence(0, bandNames.length().subtract(1))

    return ee.ImageCollection(rescaled).toBands().select(bandSequence, bandNames)


@decorators.carry_metadata
def power_to_db(img):
    return ee.Image(10).multiply(img.log10())


@decorators.carry_metadata
def db_to_power(img):
    return ee.Image(10).pow(img.divide(10))


@decorators.carry_metadata
def add_indices(img):
    ndvi = img.normalizedDifference(["nir", "red"]).rename("ndvi")
    mndwi = img.normalizedDifference(["green", "swir1"]).rename("mndwi")
    # nwi = img.expression(
    #     "((b-(n+s+w))/(b+(n+s+w))*100)",
    #     {
    #         "b": img.select("blue"),
    #         "n": img.select("nir"),
    #         "s": img.select("swir1"),
    #         "w": img.select("swir2"),
    #     },
    # ).rename("nwi")
    # aewinsh = img.expression(
    #     "4.0 * (g-s) - ((0.25*n) + (2.75*w))",
    #     {
    #         "g": img.select("green"),
    #         "s": img.select("swir1"),
    #         "n": img.select("nir"),
    #         "w": img.select("swir2"),
    #     },
    # ).rename("aewinsh")
    # aewish = img.expression(
    #     "b+2.5*g-1.5*(n+s)-0.25*w",
    #     {
    #         "b": img.select("blue"),
    #         "g": img.select("green"),
    #         "n": img.select("nir"),
    #         "s": img.select("swir1"),
    #         "w": img.select("swir2"),
    #     },
    # ).rename("aewish")
    # tcwet = img.expression(
    #     "0.1509*b + 0.1973*g + 0.3279*r + 0.3406*n - 0.7112*s - 0.4572*w",
    #     {
    #         "b": img.select("blue"),
    #         "g": img.select("green"),
    #         "r": img.select("red"),
    #         "n": img.select("nir"),
    #         "s": img.select("swir1"),
    #         "w": img.select("swir2"),
    #     },
    # ).rename("tcwet")

    out_img = ee.Image.cat([img,ndvi, mndwi])#, nwi, aewinsh, aewish, tcwet])

    return out_img
