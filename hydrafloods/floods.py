import ee

from hydrafloods import decorators


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
