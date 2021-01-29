import ee
from ee.ee_exception import EEException
import math
import datetime
from hydrafloods import decorators


@decorators.carry_metadata
def fuzzy_gaussian(img, midpoint, spread):
    """Fuzzy membership function through a Gaussian or normal distribution based around a 
    user-specified midpoint (which is assigned a membership of 1) with a defined spread.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzygaussian-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        midpoint (float): user-defined value with a fuzzy membership of 1
        spread (float): the spread of the Gaussian function. The spread generally 
            ranges from 0.01 to 1, with the larger the value results in a steeper 
            distribution around the midpoint

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    bandNames = img.bandNames()

    spread = ee.Image.constant(spread)
    midpoint = ee.Image.constant(midpoint)
    e = ee.Image.constant(math.e)

    gauss = e.expression(
        "e ** (-s * (x-m)**2)", {"e": e, "x": img, "s": spread, "m": midpoint}
    )

    return gauss.rename(bandNames)


@decorators.carry_metadata
def fuzzy_near(img, midpoint, spread):
    """Fuzzy membership function around a specific value which is defined by a user-defined midpoint
     (which is assigned a membership of 1), with a defined spread decreasing to zero.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzylarge-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        midpoint (float): user-defined value with a fuzzy membership of 1
        spread (float): the spread of the Near function. The spread generally 
            ranges from 0.001 to 1, with the larger the value results in a steeper 
            distribution from the midpoint

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    bandNames = img.bandNames()

    spread = ee.Image.constant(spread)
    midpoint = ee.Image.constant(midpoint)
    one = ee.Image.constant(1)

    near = one.expression(
        "o / (o + s * (x - m)**2)", {"o": one, "x": img, "s": spread, "m": midpoint}
    )

    return near.rename(bandNames)


@decorators.carry_metadata
def fuzzy_large(img, midpoint, spread):
    """Fuzzy membership function where the larger input values have membership closer to 1. 
    The function is defined by a user-specified midpoint (which is assigned a membership of 0.5) 
    with a defined spread.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzylarge-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        midpoint (float): user-defined value with a fuzzy membership of 0.5
        spread (float): the spread of the Large function. The spread generally 
            ranges from 1 to 10, with the larger the value results in a steeper 
            distribution from the midpoint

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    bandNames = img.bandNames()

    spread = ee.Image.constant(spread)
    midpoint = ee.Image.constant(midpoint)
    one = ee.Image.constant(1)

    large = one.expression(
        "o / ( o + (x / m) ** (s * -1))",
        {"o": one, "x": img, "s": spread, "m": midpoint},
    )

    return large.rename(bandNames)


@decorators.carry_metadata
def fuzzy_mslarge(img, mean_scaling, std_scaling, region=None, scale=90):
    """Fuzzy membership through a function based on the mean and standard deviation, 
    with the larger values having a membership closer to 1.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzymslarge-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        mean_scaling (float): multiplier for the mean of the input values in the MSLarge function equation.
        std_scaling (float): multiplier for the standard deviation of the input values in the MSLarge function equation.
        region (ee.Geometry | None, optional): region to calculate mean/stdDev image statistics, 
            if set to None will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    if region is None:
        region = img.geometry()

    bandNames = img.bandNames()

    mean_scale = ee.Image.constant(mean_scaling)
    std_scale = ee.Image.constant(std_scaling)
    one = ee.Image.constant(1)

    reducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(), None, True)

    stats = ee.Dictionary(
        img.reduceRegion(reducer, region, scale, bestEffort=True, tileScale=16)
    ).toImage()

    mslarge = one.expression(
        "o - (b * s ) / (x - (a * m) + (b * s))",
        {
            "o": one,
            "x": img,
            "a": stats.select(".*(mean)$"),
            "b": stats.select(".*(stdDev)$"),
            "m": mean_scale,
            "s": std_scale,
        },
    )

    mask = img.gte(stats.select(".*(mean)$").multiply(mean_scale))

    return mslarge.multiply(mask).rename(bandNames)


@decorators.carry_metadata
def fuzzy_small(img, midpoint, spread):
    """Fuzzy membership function with the smaller input values having membership closer to 1. 
    The function is defined by a user-specified midpoint (which is assigned a membership of 0.5) 
    with a defined spread.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzysmall-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        midpoint (float): user-defined value with a fuzzy membership of 0.5
        spread (float): the spread of the Small function. The spread generally 
            ranges from 1 to 10, with the larger the value results in a steeper 
            distribution from the midpoint

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    bandNames = img.bandNames()

    spread = ee.Image.constant(spread)
    midpoint = ee.Image.constant(midpoint)
    one = ee.Image.constant(1)

    small = one.expression(
        "o / (o + (x / m) ** s)", {"o": one, "x": img, "s": spread, "m": midpoint}
    )

    return small.rename(bandNames)


@decorators.carry_metadata
def fuzzy_mssmall(img, mean_scale, std_scale, region=None, scale=90):
    """Fuzzy membership through a function based on the mean and standard deviation, 
    with the smaller values having a membership closer to 1.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzymssmall-class.htm

    args:
        img (ee.Image): The input image whose values will be scaled from 0 to 1
        mean_scaling (float): multiplier for the mean of the input values in the MSSmall function equation.
        std_scaling (float): multiplier for the standard deviation of the input values in the MSSmall function equation.
        region (ee.Geometry | None, optional): region to calculate mean/stdDev image statistics, 
            if set to None will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90

    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    if region is None:
        region = img.geometry()

    bandNames = img.bandNames()

    mean_scale = ee.Image.constant(mean_scale)
    std_scale = ee.Image.constant(std_scale)
    one = ee.Image.constant(1)

    reducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(), None, True)

    stats = ee.Dictionary(
        img.reduceRegion(reducer, region, scale, bestEffort=True, tileScale=16)
    ).toImage()

    mssmall = one.expression(
        "(b * s ) / (x - (a * m) + (b * s))",
        {
            "o": one,
            "x": img,
            "a": stats.select(".*(mean)$"),
            "b": stats.select(".*(stdDev)$"),
            "m": mean_scale,
            "s": std_scale,
        },
    )

    mask = img.gte(stats.select(".*(mean)$").multiply(mean_scale))

    return mssmall.multiply(mask).rename(bandNames)


@decorators.carry_metadata
def fuzzy_linear(img, minimum, maximum):
    """Fuzzy membership function through a linear transformation between the user-specified minimum 
    value, a membership of 0, to the user-defined maximum value, which is assigned a membership of 1.
    If the minimum value is less than the maximum, the linear function will have a positive slope. 
    If the minimum value is greater than the maximum, the linear function will have a negative slope.
    https://desktop.arcgis.com/en/arcmap/10.3/analyze/arcpy-spatial-analyst/fuzzylinear-class.htm
        
    args:
        img (ee.Image): The input image whose values will be scaled linearly from 0 to 1
        minimum (float): The value that will have a membership of 0
        maximum (float): The value that will have a membership of 1
    
    returns:
        ee.Image: output floating-point image with values ranging from 0 to 1
    """
    bandNames = img.bandNames()

    invert = minimum > maximum

    if invert:
        min_val = ee.Number(maximum)
        max_val = ee.Number(minimum)
    else:
        min_val = ee.Number(minimum)
        max_val = ee.Number(maximum)

    linear = img.unitScale(min_val, max_val).clamp(min_val, max_val)

    if invert:
        linear = ee.Image.constant(1).subtract(linear)

    return linear.rename(bandNames)


def fuzzy_or(img_list):
    """Fuzzy Or overlay returning the maximum value of the input images.
    This technique is useful when you want to identify the highest membership values 
    for any of the input criteria

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    in_img = ee.Image.cat([img_list])

    out_img = in_img.reduce(ee.Reducer.max()).rename("fuzzy_or")

    return out_img


def fuzzy_and(img_list):
    """Fuzzy And overlay returning the minimum value of the input images.
    This technique is useful when you want to identify the least common denominator
    for the membership of all the input criteria

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    in_img = ee.Image.cat([img_list])

    out_img = in_img.reduce(ee.Reducer.min()).rename("fuzzy_and")

    return out_img


def fuzzy_product(img_list):
    """Fuzzy Product overlay which multiplies each of the fuzzy values togetherfor all the input images.
    The resulting product will be less than any of the input, and when a member of many sets is input,
    the value can be very small.

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    in_img = ee.Image.cat([img_list])

    out_img = in_img.reduce(ee.Reducer.product()).rename("fuzzy_product")

    return out_img


def fuzzy_sum(img_list):
    """Fuzzy Product overlay add the fuzzy values for all the input images.
    The resulting sum is an increasing linear combination function that is based on the number of 
    criteria entered into the analysis

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    one = ee.Image.constant(1)

    in_img = one.subtract(ee.Image.cat([img_list]))

    out_img = one.subtract(in_img.reduce(ee.Reducer.product())).rename("fuzzy_sum")

    return out_img


def fuzzy_gamma(img_list, gamma=0.5):
    """Fuzzy Gamma is an algebraic product of fuzzy Product and fuzzy Sum, 
    which are both raised to the power of gamma. If the specified gamma is 1, 
    the output is the same as fuzzy Sum; if gamma is 0, the output is the same as fuzzy Product

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together
        gamma (float): the gamma value to be used when weighting Product vs Sum

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    gamma = ee.Image.constant(gamma)
    one = ee.Image.constant(1)
    a = fuzzy_sum(img_list)
    b = fuzzy_product(img_list)

    out_img = a.pow(gamma).multiply(b.pow(one.subtract(gamma))).rename("fuzzy_gamma")

    return out_img


def fuzzy_weighted(img_list, weights):
    """Overlays several rasters, multiplying each by their given weight and summing them together.

    args:
        img_list (list[ee.Image]): list of ee.Image objects to overlay together
        weights (list[float]): list of weight values to assign images, must be in order

    returns:
        ee.Image: output floating-point image with overlayed values
    """
    weights = ee.Image.constant(weights)

    in_img = ee.Image.cat([img_list])

    out_img = in_img.multiply(weights).reduce(ee.Reducer.sum()).rename("fuzzy_weighted")

    return out_img
