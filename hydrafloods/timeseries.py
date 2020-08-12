import ee
import math
import copy
from functools import partial
from hydrafloods import decorators, datasets


def add_time_band(img, offset="year"):
    t = img.date()
    tDiff = t.difference(ee.Date("1970-01-01"), offset)
    time = ee.Image(tDiff).float().rename("time")
    return img.addBands(time)


def prep_inputs(collection):
    first = ee.Image(collection.first())
    bands = first.bandNames().getInfo()
    outCollection = copy.deepcopy(collection)
    if "time" not in bands:
        outCollection = outCollection.map(add_time_band)
    if "constant" not in bands:
        outCollection = outCollection.map(lambda x: x.addBands(ee.Image(1)))

    return outCollection


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


def fit_harmonic_trend(
    collection, n_cycles=2, independents=["constant", "time"], dependent=None
):
    if isinstance(
        collection,
        (
            datasets.Dataset,
            datasets.Viirs,
            datasets.Modis,
            datasets.Landsat8,
            datasets.Sentinel2,
            datasets.Sentinel1,
        ),
    ):
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

    harmonic_collection = prep_inputs(collection).map(_add_coefs)

    harmonic_trend = harmonic_collection.select(independents.add(dependent)).reduce(
        ee.Reducer.linearRegression(numX=independents.length(), numY=1)
    )

    # Turn the array image into a multi-band image of coefficients
    harmonic_coefficients = (
        harmonic_trend.select("coefficients")
        .arrayProject([0])
        .arrayFlatten([independents])
    )

    return harmonic_coefficients


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
    img = ee.Image(1).set("system:time_start", t.millis())
    time_band = (
        ee.Image(img.date().difference(ee.Date("1970-01-01"), "year"))
        .float()
        .rename("time")
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
