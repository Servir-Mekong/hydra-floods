import ee
import math
import copy
from functools import partial
from hydrafloods import decorators, datasets


def add_time_band(img, offset="year",apply_mask=False):
    """Function to add time band to image. Expects image has `system:time_start` property.
    Added band will have name "time"

    args:
        img (ee.Image): image with `system:time_start` to add time band too
        offset (str, optional): units to calculate from 1970-01-01. default = "year"
        apply_mask (bool, optional): boolean switch to apply image mask to time band. if False, then time band will have full coverage

    returns:
        ee.Image: image object with time band added
    """
    t = img.date()
    tDiff = t.difference(ee.Date("1970-01-01"), offset)
    time = ee.Image(tDiff).float().rename("time")
    if apply_mask:
        time = time.updateMask(img.select([0]).mask())
    return img.addBands(time)


def prep_inputs(collection, keep_bands=None, apply_mask=False):
    """Helper function to prepare inputs into time series algorithms.
    Will check if each image in collection has a time and constant band, if not it will add to collection
    Wraps `add_time_band()` for adding the time band

    args:
        collection (ee.ImageCollection): image collection to check and add time/constant band too
        keep_bands (list[str] | None, optional): regex name or list of band names to drop during prep and include in the result
            Will check if "time" or "constant" in list, if not then add to collection. default = None
        apply_mask (bool, optional): boolean switch to apply image mask to time band. if False, then time band will have full coverage

    returns:
        ee.ImageCollection: collection with time and constant bands included
    """
    if not isinstance(collection, ee.ImageCollection):
        collection = collection.collection
    first = ee.Image(collection.first())
    outCollection = copy.deepcopy(collection)

    if keep_bands is None:
        keep_bands = []

    tband = partial(add_time_band,apply_mask=apply_mask)

    if "time" not in keep_bands:
        outCollection = outCollection.map(tband)
    if "constant" not in keep_bands:
        outCollection = outCollection.map(lambda x: x.addBands(ee.Image(1)))

    out_band_order = ee.List(["constant", "time"]).cat(keep_bands)
    return outCollection.select(out_band_order)


def _get_names(base, n):
    """Private helper function to get band names based on base and number of iterations
    i.e. band_1, band_2, ... , band_n
    """
    return ee.List([f"{base}_{i}" for i in range(1, n + 1)])


def add_harmonic_coefs(image, n_cycles=2):
    """Function to add harmonic coefficients as bands to images
    Harmonic coefficients are calculated as sin and cos of frequency within year

    args:
        image (ee.Image): image object to add harmonic coefficiencts to. Expects that image has time band
        n_cycles (int, optional): number of interannual cycles to include. default = 2

    returns:
        ee.Image: image with harmonic coefficient bands added
    """
    cosNames = _get_names("cos", n_cycles)
    sinNames = _get_names("sin", n_cycles)
    frequencyImg = ee.Image.constant(ee.List.sequence(1, n_cycles))
    timeRadians = image.select("time").multiply(2 * math.pi)
    cosines = timeRadians.multiply(frequencyImg).cos().rename(cosNames)
    sines = timeRadians.multiply(frequencyImg).sin().rename(sinNames)

    return image.addBands(cosines).addBands(sines)


def fit_linear_trend(
    collection,
    independents=["constant", "time"],
    dependent=None,
    regression_method="ols",
    output_err=False,
):
    """Function to fit a linear reducer on image collection along the time dimension.

    args:
        collection (ee.ImageCollection): image collection to fit trend line for
        independents (list[str], optional): list of band names to use for fitting regression. default = ["constant", "time"]
        dependent (str | None, optional): band name of values to fit, if None then uses the first band. default = None
        regression_method (str, optional): name of regression reducer to use. options are "simple", "ols", "robust", "ridge", "sen". default = "ols"
        output_err (bool, optional): switch to output regression x and y errors, if true output will have "mean_x", "residual_x" and "residual_y" bands. 
            Useful for estimating confidence intervals. default = False

    returns:
        ee.Image: output image with regression coeffients as bands named "offset" and "slope". 
            Will have "mean_x", "residual_x" and "residual_y" if ouput_err == True
    """

    independents = ee.List(independents)

    if dependent is None:
        dependent = ee.String(ee.Image(collection.first()).bandNames().get(0))
    else:
        dependent = ee.String(dependent)

    methods = ["simple", "ols", "robust", "ridge", "sen"]
    reducers = [
        ee.Reducer.linearFit(),
        ee.Reducer.linearRegression(independents.length(), 1),
        ee.Reducer.robustLinearRegression(independents.length(), 1),
        ee.Reducer.ridgeRegression(independents.length(), 1),
        ee.Reducer.sensSlope(),
    ]
    method_lookup = {v: reducers[i] for i, v in enumerate(methods)}
    if regression_method not in methods:
        raise ValueError()
    else:
        lin_reducer = method_lookup[regression_method]

    if not isinstance(collection, ee.ImageCollection):
        collection = collection.collection

    reduction_coll = prep_inputs(
        collection, keep_bands=[dependent], apply_mask=True
    ).select(independents.add(dependent))

    n = reduction_coll.select(dependent).reduce("count").rename("n")

    if regression_method in ["sen", "ridge"]:
        independents = ee.List(["constant"]).cat(independents)

    if regression_method == "sen":
        lr = reduction_coll.reduce(lin_reducer, 16).select(
            ["offset", "slope"], independents
        )
    else:
        lr_arr = reduction_coll.reduce(lin_reducer, 16)
        lr = (
            lr_arr.select(["coefficients"])
            .arrayProject([0])
            .arrayFlatten([independents])
        )

    if output_err and regression_method is not "sen":
        linear_mse = (
            lr_arr.select("residuals").arrayProject([0]).arrayFlatten([["residual_y"]])
        )

        t_mean = reduction_coll.select("time").mean().rename("mean_x")
        x_err = (
            reduction_coll.select("time")
            .map(lambda x: x.subtract(t_mean).pow(2))
            .reduce(ee.Reducer.sum(),16).rename("residual_x")
        )

        lr = ee.Image.cat([lr,linear_mse,t_mean,x_err])

    return lr.addBands(n)


def fit_harmonic_trend(
    collection,
    n_cycles=2,
    independents=["constant", "time"],
    dependent=None,
    output_err=False,
):
    """Function to fit a harmonic trend on image collection along the time dimension.
    Uses ee.Reducer.linearRegression to solve for coefficients

    args:
        collection (ee.ImageCollection | hydrafloods.Dataset): image collection to fit trend line for
        n_cycles (int, optional): number of interannual cycles to model. default = 2
        independents (list[str], optional): list of band names to use for fitting regression. default = ["constant", "time"]
        dependent (str | None, optional): band name of values to fit, if None then uses the first band. default = None
        output_err (bool, optional): switch to output regression x and y errors, if true output will have "mean_x", "residual_x" and "residual_y" bands. 
            Useful for estimating confidence intervals. default = False

    returns:
        ee.Image: output image with regression coeffients as bands named "constant", "time", "cos_n", and "sin_n" where n is a sequnces of cycles. 
            Will have "mean_x", "residual_x" and "residual_y" if ouput_err == True

    raises:
        ValueError: if collection is not of type ee.ImageCollection or hydrafloods.Dataset
    """
    if not isinstance(collection, ee.ImageCollection):
        try:
            collection = collection.collection
        except AttributeError:
            raise ValueError("collection argument expected type ee.ImageCollection or hydrafloods.Dataset,"+
                             f"got {type(collection)}")

    if dependent is None:
        dependent = ee.String(ee.Image(collection.first()).bandNames().get(0))
    else:
        dependent = ee.String(dependent)

    independents = (
        ee.List(independents)
        .cat(_get_names("cos", n_cycles))
        .cat(_get_names("sin", n_cycles))
    )

    _add_coefs = partial(add_harmonic_coefs, n_cycles=n_cycles)

    harmonic_collection = prep_inputs(
        collection, keep_bands=[dependent], apply_mask=True
    ).map(_add_coefs)

    harmonic_trend = harmonic_collection.select(independents.add(dependent)).reduce(
        ee.Reducer.linearRegression(numX=independents.length(), numY=1), 16
    )

    n = harmonic_collection.select(dependent).reduce(ee.Reducer.count(),16).rename("n")

    # Turn the array image into a multi-band image of coefficients
    harmonic_coefficients = (
        harmonic_trend.select("coefficients")
        .arrayProject([0])
        .arrayFlatten([independents])
    ).addBands(n)

    if output_err:
        harmonic_mse = (
            harmonic_trend.select("residuals")
            .arrayProject([0])
            .arrayFlatten([["residual_y"]])
        )

        t_mean = harmonic_collection.select("time").reduce(ee.Reducer.mean(),16).rename("mean_x")
        x_err = (
            harmonic_collection.select("time")
            .map(lambda x: x.subtract(t_mean).pow(2))
            .reduce(ee.Reducer.sum(),16).rename("residual_x")
        )

        harmonic_coefficients = ee.Image.cat(
            [harmonic_coefficients, harmonic_mse, t_mean, x_err]
        )

    return harmonic_coefficients



def predict_harmonics(
    collection, harmonics, n_cycles=2, independents=["constant", "time"]
):
    """Helper function to apply harmonic trend prediction

    args:
        collection (ee.ImageCollection): image collection to apply prediction on
        harmonics (ee.Image): harmonic coefficient image with coeffiecient bands
        n_cycles (int, optional): number of interannual cycles to model, note n_cycles must equal 
            the number of cycle coefficients in `harmonics`. default = 2
        independents (list[str], optional): list of independent band names to use in model. These are
            other than the harmonic coefficient names. default = ["constant", "time"]

    returns
        ee.ImageCollection: image collection with predicted values
    """
    @decorators.carry_metadata
    def _apply_prediction(image):
        """Closure function to apply prediction
        """
        return (
            image.select(independents)
            .multiply(harmonics)
            .reduce("sum")
            .rename("predicted")
        )

    independents = (
        ee.List(independents)
        .cat(_get_names("cos", n_cycles))
        .cat(_get_names("sin", n_cycles))
    )

    _add_coefs = partial(add_harmonic_coefs, n_cycles=n_cycles)

    harmonic_collection = prep_inputs(collection).map(_add_coefs)

    # Compute fitted values.
    predicted = harmonic_collection.map(_apply_prediction)

    return predicted


def get_dummy_img(t):
    """Helper function to get an image readily available to use for predictions
    Resulting image will include time based on year and constant band of 1

    args:
        t (str | ee.Date): string or ee.Date string to create dummy image for

    returns:
        ee.Image: image object with time and constant bands
    """
    if not isinstance(t, ee.Date):
        t = ee.Date(t)
    img = ee.Image.constant(1).set("system:time_start", t.millis())
    time_band = (
        ee.Image(t.difference(ee.Date("1970-01-01"), "year")).float().rename("time")
    )
    return img.addBands(time_band)


def get_dummy_collection(start_time, end_time):
    """Helper function to create image collection of images to apply time series prediction on
    Creates daily imagery between start_time and end_time

    args:
        start_time (str | ee.Date): string or ee.Date string to start creating images on (inclusive)
        end_time (str | ee.Date): string or ee.Date string to end creating images on (exclusive)

    returns:
        ee.ImageCollection: image collection with image containing time and constant bands
    """
    def _gen_image(i):
        t = start_time.advance(ee.Number(i), "day")
        return ee.Image().rename("blank").set("system:time_start", t.millis())

    if not isinstance(start_time, ee.Date):
        start_time = ee.Date(start_time)
    if not isinstance(end_time, ee.Date):
        end_time = ee.Date(end_time)

    n = end_time.difference(start_time, "day")
    coll = ee.ImageCollection(ee.List.sequence(0, n).map(_gen_image))

    return prep_inputs(coll)

def temporal_smoothing(collection,reducer,days=10):
    """Function to apply moving window reducer in time on image collection
    
    args:
        collection (ee.ImageCollection): image collection to apply moving window reducer in time on
        reducer (ee.Reducer): earth engine reducer object to apply
        days (int,optional): size of moving time window in days to apply reducer. default = 10

    returns:
        ee.ImageCollection: image collection with reducer applied in time
    """
    @decorators.carry_metadata
    def _smooth(img):
        """Closure function to apply smoothing in between window
        """
        t = img.date()
        band_names = img.bandNames()
        t_start = t.advance(-days//2, "day")
        t_stop = t.advance(days//2, "day")
        return collection.filterDate(t_start,t_stop).reduce(reducer,8).rename(band_names)

    return collection.map(_smooth)

def temporal_iqr_filter(collection):
    """Function to filter values outside of IQR in time on image collection
    
    args:
        collection (ee.ImageCollection): image collection to apply filter in time on

    returns:
        ee.ImageCollection: image collection with values filtered 
    """
    @decorators.carry_metadata
    def _filter(img):
        """Closure function to apply smoothing in between window
        """
        mask = img.gt(low_bounds).And(img.lt(upper_bounds))
        return img.updateMask(mask)

    percentiles = collection.reduce(ee.Reducer.percentile([25,75]))
    iqr = percentiles.select(1).subtract(percentiles.select(0))
    low_bounds = percentiles.select(0).subtract(iqr.multiply(1.5))
    upper_bounds = percentiles.select(1).add(iqr.multiply(1.5))

    return collection.map(_filter)