import ee
import math
import copy
from functools import partial
from hydrafloods import decorators, datasets


def add_time_band(img, offset="year",apply_mask=False):
    t = img.date()
    tDiff = t.difference(ee.Date("1970-01-01"), offset)
    time = ee.Image(tDiff).float().rename("time")
    if apply_mask:
        time = time.updateMask(img.select([0]).mask())
    return img.addBands(time)


def prep_inputs(collection, keep_bands=None, apply_mask=False):
    if not isinstance(collection, ee.ImageCollection):
        collection = collection.collection
    first = ee.Image(collection.first())
    outCollection = copy.deepcopy(collection)

    tband = partial(add_time_band,apply_mask=apply_mask)

    if "time" not in keep_bands:
        outCollection = outCollection.map(tband)
    if "constant" not in keep_bands:
        outCollection = outCollection.map(lambda x: x.addBands(ee.Image(1)))

    out_band_order = ee.List(["constant", "time"]).cat(keep_bands)
    return outCollection.select(out_band_order)


def _get_names(base, n):
    return ee.List([f"{base}_{i}" for i in range(1, n + 1)])


def add_harmonic_coefs(image, n_cycles=2):
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
    if not isinstance(collection, ee.ImageCollection):
        collection = collection.collection

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


def calc_xy_errs(img, regression_img=None):
    t = img.date()
    dum = get_dummy_img(t)
    tPred = dum.multiply(regression_img).reduce("sum")
    xerr = img.select("time").subtract(tmean).pow(2).rename("x")
    yerr = img.select(label).subtract(tPred).pow(2).rename("y")
    return xerr.addBands(yerr)


def predict_harmonics(
    collection, harmonics, n_cycles=2, independents=["constant", "time"]
):
    @decorators.carry_metadata
    def _apply_prediction(image):
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
    if not isinstance(t, ee.Date):
        t = ee.Date(t)
    img = ee.Image.constant(1).set("system:time_start", t.millis())
    time_band = (
        ee.Image(t.difference(ee.Date("1970-01-01"), "year")).float().rename("time")
    )
    return img.addBands(time_band)


def get_dummy_collection(t1, t2):
    def _gen_image(i):
        t = t1.advance(ee.Number(i), "day")
        return ee.Image().rename("blank").set("system:time_start", t.millis())

    if not isinstance(t1, ee.Date):
        t1 = ee.Date(t1)
    if not isinstance(t2, ee.Date):
        t2 = ee.Date(t2)

    n = t2.difference(t1, "day")
    coll = ee.ImageCollection(ee.List.sequence(0, n).map(_gen_image))

    return prep_inputs(coll)
