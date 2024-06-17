# module to host index functions for image transformations
import ee
from hydrafloods import decorators


@decorators.keep_attrs
def ndvi(img):
    """Function to calculate Normalized Difference Vegetation Index (NDVI).
    Expects image has "nir" and "red" bands.

    args:
        img (ee.Image): image to calculate NDVI

    returns:
        ee.Image: NDVI image
    """
    return img.normalizedDifference(["nir", "red"]).rename("ndvi")


@decorators.keep_attrs
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


@decorators.keep_attrs
def mndwi(img):
    """Function to calculate Modified Difference Water Index (MNDWI).
    Expects image has "green" and "swir1" bands.

    args:
        img (ee.Image): image to calculate MNDWI

    returns:
        ee.Image: MNDWI image
    """
    return img.normalizedDifference(["green", "swir1"]).rename("mndwi")


@decorators.keep_attrs
def nwi(img):
    """Function to calculate New Water Index (NWI).
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


@decorators.keep_attrs
def gwi(img):
    """Function to calculate General Water Index (GWI)
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


@decorators.keep_attrs
def aewinsh(img):
    """Function to calculate Automated Water Extraction Index no shadow (AEWInsh)
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


@decorators.keep_attrs
def aewish(img):
    """Function to calculate Automated Water Extraction Index shadow (AEWIsh)
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


@decorators.keep_attrs
def lswi(img):
    """Function to calculate Land Surface Water Index (LSWI).
    Expects image has "nir" and "swir1" bands.

    args:
        img (ee.Image): image to calculate LSWI

    returns:
        ee.Image: LSWI image
    """
    return img.expression(
        "(nir-swir)/(nir+swir)", {"nir": img.select("nir"), "swir": img.select("swir1")}
    ).rename("lswi")


@decorators.keep_attrs
def wri(img):
    """Function to calculate Water Ratio Index (WRI).
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


@decorators.keep_attrs
def mbwi(img, factor=3):
    """Function to calculate Multi Band Water Index (MBWI).
    Expects image has "green", "red", "nir", "swir1", and "swir2" bands.
    https://doi.org/10.1016/j.jag.2018.01.018

    args:
        img (ee.Image): image to calculate MBWI
        factor (int,optional): factor to scale green band for index. default=3

    returns:
        ee.Image: MBWI image
    """
    return img.expression(
        "((factor*green)-red-nir-swir1-swir2)",
        {
            "factor": factor,
            "green": img.select("green"),
            "red": img.select("red"),
            "nir": img.select("nir"),
            "swir1": img.select("swir1"),
            "swir2": img.select("swir2"),
        },
    ).rename("mbwi")


@decorators.keep_attrs
def mwi(img):
    """Function to calculate the Modified Water Index (MWI)
    Expect image has "green", "red", "nir", and "swir1" bands
    https://doi.org/10.1007/978-3-662-45737-5_51 , https://ieeexplore.ieee.org/document/9011209
    """
    # calculate ndvi
    ndvi_img = ndvi(img)
    # calculate mndwi
    mndwi_img = mndwi(img)

    # mwi = 1 - (ndvi - mndwi)
    return ee.Image.constant(1).subtract(ndvi_img.subtract(mndwi_img)).rename("mwi")


@decorators.keep_attrs
def rvi(img):
    """Function to calculate SAR Radar Vegetation Index (RVI) index.
    See more here: https://pro.arcgis.com/en/pro-app/latest/help/analysis/raster-functions/sar-indices-function.htm
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate RVI

    returns:
        ee.Image: RVI image
    """
    return img.expression(
        "(4*VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("rvi")


@decorators.keep_attrs
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


@decorators.keep_attrs
def vv_vh_abs_sum(img):
    """Function to calculate the absolute value of the sum of VV and VH bands.
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to apply calculation

    returns:
        ee.Image: image name 'vv_vh_abs_sum'
    """
    return img.select("VV").add(img.select("VH")).abs().rename("vv_vh_abs_sum")


@decorators.keep_attrs
def ndpi(img):
    """Function to calculate Normalized Difference Polarization Index (NDPI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NDPI

    returns:
        ee.Image: NPDI image
    """
    return img.expression(
        "(VV-VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("ndpi")


@decorators.keep_attrs
def nvvi(img):
    """Function to calculate Normalized VV Index (NVVI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NVVI

    returns:
        ee.Image: NVVI image
    """
    return img.expression(
        "(VV)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("nvvi")


@decorators.keep_attrs
def nvhi(img):
    """Function to calculate Normalized VH Index (NVHI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NVHI

    returns:
        ee.Image: NVHI image
    """
    return img.expression(
        "(VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("nvhi")
