# module to host index functions for image transformations
import ee
from ee.ee_exception import EEException
import math
import string
import random
import datetime
from hydrafloods import decorators


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
def mbwi(img,factor=3):
    """Function to calculate multi band water index (MBWI).
    Expects image has "green", "red", "nir", "swir1", and "swir2" bands.

    args:
        img (ee.Image): image to calculate MBWI
        factor (int,optional): factor to scale green band for index. default=2

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
            "swir2": img.select("swir2")
        },
    ).rename("mbwi")


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
def nvvi(img):
    """Function to calculate nomalized VV index (NVVI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NVVI

    returns:
        ee.Image: NVVI image
    """
    return img.expression(
        "(VV)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("nvvi")


@decorators.carry_metadata
def nvhi(img):
    """Function to calculate nomalized VH index (NVHI).
    Expects image has "VV" and "VH" bands.

    args:
        img (ee.Image): image to calculate NVHI

    returns:
        ee.Image: NVHI image
    """
    return img.expression(
        "(VH)/(VV+VH)", {"VV": img.select("VV"), "VH": img.select("VH")}
    ).rename("nvhi")
