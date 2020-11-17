import ee
from ee.ee_exception import EEException
import math
import string
import random
import datetime
from hydrafloods import decorators


# helper function to convert qa bit image to flag
def extract_bits(image, start, end=None, new_name=None):
    """Function to conver qa bits to binary flag image

    args:
        image (ee.Image): qa image to extract bit from
        start (int): starting bit for flag
        end (int | None, optional): ending bit for flag, if None then will only use start bit. default = None
        new_name (str | None, optional): output name of resulting image, if None name will be {start}Bits. default = None

    returns:
        ee.Image: image with extract bits
    """

    newname = new_name if new_name is not None else f"{start}Bits"

    if (start == end) or (end is None):
        # perform a bit shift with bitwiseAnd
        return image.select([0], [newname]).bitwiseAnd(1 << start)
    else:
        # Compute the bits we need to extract.
        pattern = 0
        for i in range(start, end):
            pattern += int(math.pow(2, i))

        # Return a single band image of the extracted QA bits, giving the band
        # a new name.
        return image.select([0], [newname]).bitwiseAnd(pattern).rightShift(start)


def get_geoms(img):
    """Helper function to get geometry from image

    args:
        img (ee.Image): image to get geometry from
    
    returns:
        ee.Geometry: geometry of image
    """
    return img.geometry()


def export_image(
    image,
    region,
    asset_id,
    description=None,
    scale=1000,
    crs="EPSG:4326",
    pyramiding=None,
):
    """Function to wrap image export with EE Python API

    args:
        image (ee.Image): image to export
        region (ee.Geometry): region to export image
        asset_id (str): asset ID to export image to
        description (str | None, optional): description to identify image export/
            if None then description will be random string. default = None
        scale (int, optional): resolution in meters to export image to. default = 1000
        crs (str, optional): epsg code to export image to. default = "EPSG:4326"
        pyramiding (dict | None, optional): dictionary defining band pyramiding scheme.
            if None then "mean" will be used as default for all bands. default = None

    """
    if (description == None) or (type(description) != str):
        description = "".join(
            random.SystemRandom().choice(string.ascii_letters) for _ in range(8)
        ).lower()
    # get serializable geometry for export
    export_region = region.bounds(maxError=10).getInfo()["coordinates"]

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


def batch_export(
    collection,
    collection_asset,
    region=None,
    prefix=None,
    suffix=None,
    scale=1000,
    crs="EPSG:4326",
    pyramiding=None,
    metadata=None,
    verbose=False,
):
    """Function to export each image in a collection
    Wraps `export_image` will set YYYYMMdd formatted time in file name

    args:
        collection (ee.ImageCollection): image collection to export
        collection_asset (str): image collection asset ID to export to
        region (ee.Geometry): region to export image
        prefix (str): prefix string to add before time info in name
        suffix (str): suffix string to add after time info in name
        scale (int, optional): resolution in meters to export image to. default = 1000
        crs (str, optional): epsg code to export image to. default = "EPSG:4326"
        pyramiding (dict | None, optional): dictionary defining band pyramiding scheme.
            if None then "mean" will be used as default for all bands. default = None
        metadata (dict | None, optional):
        verbose (bool, optional):

    """
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

        if not collection_asset.endswith("/"):
            collection_asset += "/"

        exportName = collection_asset + description

        export_image(
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
def rescale(img, scale=0.0001, offset=0):
    """Function to linearly rescale units using user defined scale and offset

    args:
        img (ee.Image): image to rescale
        scale (float,optional): scale value (i.e. slope of linear equation). default = 0.0001
        offset (float, optional): offset value (i.e. y-intercept). default = 0

    returns:
        ee.Image: rescaled image
    """
    return img.multiply(scale).add(offset)


@decorators.carry_metadata
def power_to_db(img):
    """Function to convert SAR units from power to dB

    args:
        img (ee.Image): SAR power image to convert to dB

    returns:
        ee.Image: dB SAR image
    """
    return ee.Image(10).multiply(img.log10())


@decorators.carry_metadata
def db_to_power(img):
    """Function to convert SAR units from dB to power

    args:
        img (ee.Image): SAR dB image to convert to power

    returns:
        ee.Image: power SAR image
    """
    return ee.Image(10).pow(img.divide(10))


@decorators.carry_metadata
def ndvi(img):
    """Function to calculate Normalized Difference Vegetation Index (NDVI).
    Expects image has "nir" and "red" bands.

    args:
        img (ee.Image): image to calculate NDVI

    returns:
        ee.Image: NDVI image
    """
    return img.normalizedDifference(["nir", "red"]).rename("ndvi")


@decorators.carry_metadata
def evi(img):
    """Function to calculate Enhanced Vegetation Index (EVI).
    Expects image has "blue", "red", and "nir" bands.

    args:
        img (ee.Image): image to calculate EVI

    returns:
        ee.Image: EVI image
    """
    return img.expression(
        "2.5*(nir-red)/(nir+6.0*red-7.5*blue+1)",
        {
            "blue": img.select("blue"),
            "red": img.select("red"),
            "nir": img.select("nir"),
        },
    ).rename("evi")


@decorators.carry_metadata
def mndwi(img):
    """Function to calculate modified Difference Water Index (MNDWI).
    Expects image has "green" and "swir1" bands.

    args:
        img (ee.Image): image to calculate MNDWI

    returns:
        ee.Image: MNDWI image
    """
    return img.normalizedDifference(["green", "swir1"]).rename("mndwi")


@decorators.carry_metadata
def nwi(img):
    """Function to calculate new water index (NWI).
    Expects image has "blue", "nir", "swir1" and "swir2" bands.

    args:
        img (ee.Image): image to calculate NWI

    returns:
        ee.Image: NWI image
    """
    return img.expression(
        "((b-(n+s+w))/(b+(n+s+w)))",
        {
            "b": img.select("blue"),
            "n": img.select("nir"),
            "s": img.select("swir1"),
            "w": img.select("swir2"),
        },
    ).rename("nwi")


@decorators.carry_metadata
def gwi(img):
    """Function to calculate general water index (GWI)
    Expects image has "green", "red", "nir", and "swir1" bands.

    args:
        img (ee.Image): image to calculate GWI

    returns:
        ee.Image: GWI image
    """
    return img.expression(
        "(g+r)-(n+s)",
        {
            "g": img.select("green"),
            "r": img.select("red"),
            "n": img.select("nir"),
            "s": img.select("swir1"),
        },
    ).rename("gwi")


@decorators.carry_metadata
def aewinsh(img):
    """Function to calculate automated water extraction index (AEWI) no shadow
    Expects image has "green", "nir", "swir1" and "swir2" bands.

    args:
        img (ee.Image): image to calculate AEWInsh

    returns:
        ee.Image:  AEWInsh image
    """
    return img.expression(
        "4.0 * (g-s) - ((0.25*n) + (2.75*w))",
        {
            "g": img.select("green"),
            "s": img.select("swir1"),
            "n": img.select("nir"),
            "w": img.select("swir2"),
        },
    ).rename("aewinsh")


@decorators.carry_metadata
def aewish(img):
    """Function to calculate automated water extraction index (AEWI) shadow
    Expects image has "blue", "green", "nir", "swir1" and "swir2" bands.

    args:
        img (ee.Image): image to calculate AEWIsh

    returns:
        ee.Image:  AEWIsh image
    """
    return img.expression(
        "b+2.5*g-1.5*(n+s)-0.25*w",
        {
            "b": img.select("blue"),
            "g": img.select("green"),
            "n": img.select("nir"),
            "s": img.select("swir1"),
            "w": img.select("swir2"),
        },
    ).rename("aewish")


@decorators.carry_metadata
def lswi(img):
    """Function to calculate land surface water index (LSWI).
    Expects image has "nir" and "swir1" bands.

    args:
        img (ee.Image): image to calculate LSWI

    returns:
        ee.Image: LSWI image
    """
    return img.expression(
        "(nir-swir)/(nir+swir)", {"nir": img.select("nir"), "swir": img.select("swir1")}
    ).rename("lswi")


@decorators.carry_metadata
def wri(img):
    """Function to calculate water ratio index (WRI).
    Expects image has "green", "red", "nir" and "swir1" bands.

    args:
        img (ee.Image): image to calculate WRI

    returns:
        ee.Image: WRI image
    """
    return img.expression(
        "(green+red)/(nir+swir)",
        {
            "green": img.select("green"),
            "red": img.select("red"),
            "nir": img.select("nir"),
            "swir": img.select("swir1"),
        },
    ).rename("wri")


@decorators.carry_metadata
def rfi(img):
    """Function to calculate SAR RFI index.
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate RFI

    returns:
        ee.Image: RFI image
    """
    return img.expression(
        "(4*VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("rfi")


@decorators.carry_metadata
def vv_vh_ratio(img):
    """Function to calculate ratio between VV and VH bands.
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate ration

    returns:
        ee.Image: ratio image name 'ratio'
    """
    return img.expression(
        "(VV/VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("ratio")


@decorators.carry_metadata
def vv_vh_abs_sum(img):
    """Function to calculate the absolute value of the sum of VV and VH bands.
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to apply calculation

    returns:
        ee.Image: image name 'vv_vh_abs_sum'
    """
    return img.select("VV").add(img.select("VH")).abs().rename("vv_vh_abs_sum")


@decorators.carry_metadata
def ndpi(img):
    """Function to calculate nomalized difference polarization index (NDPI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NDPI

    returns:
        ee.Image: NPDI image
    """
    return img.expression(
        "(VV-VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("ndpi")


@decorators.carry_metadata
def add_indices(img, indices=["mndwi"]):
    """Function to calculate multiple band indices and add to image as bands

    args:
        img (ee.Image): image to calculate indices from
        indices (list[str], optional): list of strings of index names to calculate. 
            can use any named index function in geeutils. default = ["mndwi"]

    returns:
        ee.Image: image object with added indices
    """
    # create a dict to look up index functions
    # need to watch out for namespacing...
    # does globals grab variable names too???
    local_funcs = globals()

    # loop through each index and append to images list
    cat_bands = [img]
    for index in indices:
        cat_bands.append(local_funcs[index](img))

    # return images as concatenated bands
    return ee.Image.cat(cat_bands)


def tile_region(region, grid_size=0.1, intersect_geom=None, contain_geom=None):
    """Function to create a feature collection of tiles covering a region

    args:
        region (ee.Geometry): region to create tile grid over
        grid_size (float, optional): resolution in decimal degrees to create tiles. default = 0.1
        intersect_geom (ee.Geometry | None, optional): geometry object to filter tiles that intesect with 
            geometry useful for filtering tiles that are created over oceans with no data. default = None
        contain_geom (ee.Geometry | None, optional): geometry object to filter tiles that are contained within 
            geometry useful for filtering tiles that are only in an area. default = None

    returns:
        ee.FeatureCollection: collection of feature tiles at a given grid_size over a region
    """
    # nesting grid construction along y and then x coordinates
    def constuctGrid(i):
        """Closure function to contruct grid
        """

        def contructXGrid(j):
            j = ee.Number(j)
            box = ee.Feature(
                ee.Geometry.Rectangle(
                    [j, i, j.add(grid_res), i.add(grid_res)],
                    "epsg:4326",
                    geodesic=False,
                )
            )
            if contain_geom is not None:
                out = ee.Algorithms.If(
                    box.contains(contain_geom, maxError=10), box, None
                )
            elif intersect_geom is not None:
                out = ee.Algorithms.If(
                    box.intersects(intersect_geom, maxError=10), box, None
                )
            else:
                out = box
            return ee.Feature(out)

        i = ee.Number(i)
        out = ee.List.sequence(west, east.subtract(grid_res), grid_res).map(
            contructXGrid, True
        )
        return out

    if (contain_geom is not None) and (intersect_geom is not None):
        raise ValueError(
            "contains and intersection keywords are mutually exclusive, please define only one"
        )

    bounds = region.bounds(maxError=10)
    coords = ee.List(bounds.coordinates().get(0))
    grid_res = ee.Number(grid_size)

    west = ee.Number(ee.List(coords.get(0)).get(0))
    south = ee.Number(ee.List(coords.get(0)).get(1))
    east = ee.Number(ee.List(coords.get(2)).get(0))
    north = ee.Number(ee.List(coords.get(2)).get(1))

    west = ee.Algorithms.If(
        west.lt(0),
        west.subtract(west.mod(grid_res).add(grid_res)),
        west.subtract(west.mod(grid_res)),
    )
    south = ee.Algorithms.If(
        south.lt(0),
        south.subtract(south.mod(grid_res).add(grid_res)),
        south.subtract(south.mod(grid_res)),
    )
    east = east.add(grid_res.subtract(east.mod(grid_res)))
    north = north.add(grid_res.subtract(north.mod(grid_res)))

    grid = ee.FeatureCollection(
        ee.List.sequence(south, north.subtract(grid_res), grid_res)
        .map(constuctGrid)
        .flatten()
    )

    return grid


def country_bbox(country_name, max_error=100):
    """Function to get a bounding box geometry of a country

    args:
        country_name (str): US-recognized country name
        max_error (float,optional): The maximum amount of error tolerated when 
            performing any necessary reprojection. default = 100

    returns:
        ee.Geometry: geometry of country bounding box    
    """

    all_countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
    return (
        all_countries.filter(ee.Filter.eq("country_na", country_name))
        .geometry(max_error)
        .bounds(max_error)
    )
