import ee

from hydrafloods import (
    timeseries,
    decorators,
    geeutils,
    thresholding
)


@decorators.keep_attrs
def extract_flood(observation, reference="seasonal", permanent_threshold=75):
    """Function used to extract flooded area based off of different JRC datasets.
    Expects that the water image has water = 1, land = 0, unknown = mask

    args:
        observation (ee.Image): image object representing the observation of water to extract floods
        reference (str): string to define how permanent water is defined. 'yearly' means use past 5 years
            of JRC yearly data and extract permanent class. 'seasonal' means use monthly water occurrence
            for the month of observation. 'occurrence' means use full JRC occurrence dataset irregardless of
            time. default = seasonal
        permanent_threshold (int): threshold value in % to define permanent water. Only used when
            reference equals 'seasonal' or 'occurrence'. default = 75

    returns:
        ee.Image: the extracted flood image where floods = 1

    raises:
        NotImplementedError: when user provides an invalid option for reference parameter
    """

    if reference == "yearly":
        end_time = observation.date()
        permanent_water = (
            ee.ImageCollection(
                "JRC/GSW1_2/YearlyHistory"
            )  # get the JRC historical dataset
            .filterDate(
                "1985-01-01", end_time
            )  # filter for historical data up to date of interest
            .limit(5, "system:time_start", False)  # grab the 5 latest images
            .map(
                lambda x: x.select("waterClass").eq(3)
            )  # extract out the permanent water class
            .sum()  # check if a pixel has been classified as permanent water in the past 5 years
            .unmask(0)
            .gt(0)
        )

    elif reference == "seasonal":
        month = observation.date().get("month")
        permanent_water = (
            ee.ImageCollection("JRC/GSW1_3/MonthlyRecurrence")
            .select("monthly_recurrence")
            .filter(ee.Filter.eq("month", month))
            .first()
            .gt(permanent_threshold)
        )

    elif reference == "occurrence":
        permanent_water = (
            ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
            .select("occurrence")
            .gt(permanent_threshold)
        )

    else:
        raise NotImplementedError(
            "the selected reference method, {reference}, is not implemeted. please select either 'yearly','seasonal', or 'occurrence'"
        )

    return discrete_difference(observation, permanent_water)


@decorators.keep_attrs
def discrete_difference(observation, reference):
    """Function to find the difference between an observation and reference dataset.
    The data should be discrete values (i.e. a classification). This function
    expects that water/floods = 1 and land = 0. Unknown values should be masked.

    args:
        observation (ee.Image): image object representing the observation of water to extract floods
        reference (ee.Image): image object representing reference/permanent water to extract floods from

    returns:
        ee.Image: extracted flood image where flood values = 1
    """
    og_mask = observation.mask()
    floods = observation.unmask(0).add(reference.unmask(0).multiply(2)).eq(1)
    return floods.updateMask(og_mask).rename("flood")


@decorators.keep_attrs
def lar_change_detection(
    observation, reference, band=None, in_units="dB", segmentation="edgeotsu", **kwargs
):
    """Log Amplitude Ratio change detection method. https://doi.org/10.1080/014311698215649
    Note: This method only works for SAR imagery.

    args:
        observation (ee.Image):
        reference (ee.Image):
        band (str | None,optional): band name to use for thresholding, if set to `None` will use first band in image. default = None
        in_units (str, optional): string specifying the input units for the imagery. Options are 'dB' or 'power'. If in_units = 'dB',
            then the data is converted to power units. default = dB
        segmentation (str | None, optional): segmentation method to use to . Options are 'edgeotsu', 'bmaxotsu', or None. If none
            is provided then no segmentation is applied and the raw log amplitude ratio is returned. default = edgeotsu
        **kwargs: optional keywords to pass to segmentation method, not required if segmentation = None


    returns:
        ee.Image:
    """
    if band is None:
        band = observation.bandNames().get(0)
    # convert the db data to amplitude units
    # amplitude = sqrt(power)
    # then divide post/pre and take the log
    if in_units == "dB":
        lar = (
            geeutils.db_to_power(observation).sqrt().select(band)
            .divide(geeutils.db_to_power(reference).sqrt().select(band))
            .log10()
        )

    elif in_units == "power":
        lar = observation.sqrt().select(band).divide(reference.sqrt().select(band)).log10()

    else:
        raise NotImplementedError(
            f"input units could not be infered from {in_units}, please select either 'dB' or 'power'"
        )

    if segmentation == "edgeotsu":
        floods_lar = thresholding.edge_otsu(lar, band=band, **kwargs)

    elif segmentation == "bmaxotsu":
        floods_lar = thresholding.bmax_otsu(lar, band=band, **kwargs)

    elif segmentation == None:
        floods_lar = lar

    else:
        raise NotImplementedError(
            f"selected segmentation method {segmentation} is not valid, please select 'edgeotsu', 'bmaxotsu', or None"
        )

    return floods_lar


# def temporal_zscore_change_detection(
#     collection,
#     reference_period,
# ):
#     """Flood dectection method from https://doi.org/10.1016/j.rse.2020.111664
#     Calculates the per-pixel z_score from an dryseason composite and determines
#     flooded pixels using a decision tree
#     """
#
#     return


def flood_duration(collection, is_masked=True):

    if not isinstance(collection, ee.ImageCollection):
        collection = collection.collection

    max_extent = collection.max().selfMask()

    if is_masked == True:
        f = lambda x: timeseries.add_time_band(
            x.selfMask(), apply_mask=True, offset="day"
        )

    else:
        f = lambda x: timeseries.add_time_band(x, apply_mask=True, offset="day")

    time_coll = collection.map(f)

    min_time = time_coll.select("time").min()
    max_time = time_coll.select("time").max()

    return max_time.subtract(min_time).updateMask(max_extent)
