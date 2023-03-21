from __future__ import absolute_import
import os
import ee
import math
import copy
import datetime
from pipetools import pipe
from pprint import pformat
from functools import partial
from ee.ee_exception import EEException
from hydrafloods import decorators, geeutils


class Dataset:
    """Base dataset class used to define an EE image collection by datetime and geographic region

    A dataset wraps an ee.ImageCollection by applying the spatial and temporal filtering upon init.
    Provides utility functionality to make working with and managing image collections less verbose

    Example:
        Create a dataset object for Sentinel-1 data over Alabama for 2019
        >>> ds = Dataset(
        ...    region = ee.Geometry.Rectangle([
                   -88.473227,
                   30.223334,
                   -84.88908,
                   35.008028
               ]),
        ...    start_time = "2019-01-01",
        ...    end_time = "2020-01-01",
        ...    asset_id = "COPERNICUS/S1_GRD"
        ... )
        >>> ds
        HYDRAFloods Dataset:
        {'asset_id': 'COPERNICUS/S1_GRD',
        'end_time': '2020-01-01',
        'name': 'Dataset',
        'region': [[[...], [...], [...], [...], [...]]],
        'start_time': '2019-01-01'}
    """

    def __init__(self, region, start_time, end_time, asset_id, use_qa=False):
        """Initialize Dataset class

        args:
            region (ee.Geometry): earth engine geometry object to filter image collection by
            start_time (str | datetime.datetime): start time used to filter image collection
            end_time (str | datetime.datetime): end time used to filter image collection
            asset_id (str): asset id of earth engine collection
            use_qa (bool, optional): boolean to determine to use an internal function qa().
                Used for definining custom dataset objects

        raises:
            AttributeError: if qa() method is not defined and use_qa is True
        """

        # TODO: add exceptions to check datatypes
        self.region = region  # dtype = ee.Geometry
        self.start_time = start_time
        self.end_time = end_time
        self.asset_id = asset_id
        self.use_qa = use_qa

        # dictionary mapping of band names used to harmonized optical datasets to same names
        self.BANDREMAP = ee.Dictionary(
            {
                "landsat7": ee.List(["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"]),
                "landsat8": ee.List(["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"]),
                "viirs": ee.List(["M2", "M4", "I1", "I2", "I3", "M11"]),
                "sen2": ee.List(["B2", "B3", "B4", "B8", "B11", "B12"]),
                "modis": ee.List(
                    [
                        "sur_refl_b03",
                        "sur_refl_b04",
                        "sur_refl_b01",
                        "sur_refl_b02",
                        "sur_refl_b06",
                        "sur_refl_b07",
                    ]
                ),
                "new": ee.List(["blue", "green", "red", "nir", "swir1", "swir2"]),
            }
        )

        # get the image collection and filter by geographic region and date time
        imgcollection = (
            ee.ImageCollection(self.asset_id)
            .filterBounds(self.region)
            .filterDate(self.start_time, self.end_time)
        )

        # check if to apply arbitrary qa process on collection
        # qa function can be defined in custom objects extending dataset
        if self.use_qa:
            try:
                imgcollection = imgcollection.map(self.qa)
            except AttributeError:
                raise AttributeError(
                    "qa() method is not defined...please define one or set `use_qa` to False"
                )

        self.collection = imgcollection

    def __repr__(self):
        # format datetime information
        if isinstance(self.start_time, datetime.datetime):
            ststr = self.start_time.strftime("%Y-%m-%d")
        else:
            ststr = self.start_time

        if isinstance(self.end_time, datetime.datetime):
            etstr = self.end_time.strftime("%Y-%m-%d")
        else:
            etstr = self.end_time

        # create dict of dataset information
        objDict = {
            "name": self.__class__.__name__,
            "asset_id": self.asset_id,
            "start_time": ststr,
            "end_time": etstr,
            "region": self.region.coordinates().getInfo(),
        }

        # pretty format dict and return information
        strRepr = pformat(objDict, depth=3)
        return f"HYDRAFloods Dataset:\n{strRepr}"

    @property
    def collection(self):
        """image collection object property wrapped by dataset"""
        return self._collection

    @collection.setter
    def collection(self, value):
        """setter function for collection property"""
        self._collection = value
        return

    @property
    def n_images(self):
        """Number of images contained in the dataset"""
        return self.collection.size().getInfo()

    @property
    def dates(self):
        """Dates of imagery contained in the image collection"""
        eeDates = self.collection.aggregate_array("system:time_start").map(
            lambda x: ee.Date(x).format("YYYY-MM-dd HH:mm:ss.SSS")
        )
        return eeDates.getInfo()

    @staticmethod
    def from_imgcollection(img_collection):
        """Static method to convert an ee.ImageCollection object to a hf.Dataset object.
        This method will take some time as it uses computed ee Objects from the image collection
        propeties to populate the Dataset object properties (passing info from server to client)

        args:
            img_collection (ee.ImageCollection): computed ee.ImageCollection object to create a hf.Dataset

        returns:
            hf.Dataset: dataset object with property information directly from the ee.ImageCollection
        """

        # get region and date information
        region = (
            img_collection.map(geeutils.get_geoms)
            .union(maxError=100)
            .geometry(maxError=100)
        )
        # convert ee.Date info to string format
        dates = (
            img_collection.aggregate_array("system:time_start")
            .sort()
            .map(lambda x: ee.Date(x).format("YYYY-MM-dd HH:mm:ss.S"))
        )

        # pull the date info to local strings
        start_time = dates.get(0).getInfo()
        end_time = dates.get(-1).getInfo()

        # get the collection id for `.asset_id` property
        collection_id = img_collection.get("system:id").getInfo()
        if collection_id == "None":
            collection_id = "Custom ImageCollection"

        # make a dummy dataset
        dummy_ds = Dataset(
            region,
            start_time,
            end_time,
            asset_id="NOAA/VIIRS/001/VNP09GA",
            use_qa=False,
        )

        # override the dummy dataset information with the correct data fr
        dummy_ds.asset_id = collection_id
        dummy_ds.collection = img_collection

        return dummy_ds

    def qa(
        self,
    ):
        return

    def _inplace_wrapper(self, collection, inplace):
        """Private helper function to replace the collection info for a class
        typically used when other class methods have inplace
        """
        if inplace:
            self.collection = collection
            return
        else:
            outCls = self.copy()
            outCls.collection = collection
            return outCls

    def copy(self):
        """returns a deep copy of the hydrafloods dataset class"""
        return copy.deepcopy(self)

    def apply_func(self, func, inplace=False, *args, **kwargs):
        """Wrapper method to apply a function to all of the image in the dataset.
        Makes a copy of the collection and reassigns the image collection property.
        Function must accept an ee.ImageCollection and return an ee.ImageCollection

        args:
            func (object): Function to map across image collection. Function must accept ee.Image as first argument
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False
            **kwargs: arbitrary keyword to pass to `func`

        returns:
            Dataset | None: copy of class with results from `func` as image within collection property
        """

        # get a partial function to map over imagery with the keywords applied
        # expects that the first positional arg is an
        func = partial(func, **kwargs)

        out_coll = self.collection.map(func)

        return self._inplace_wrapper(out_coll, inplace)

    def apply(self, func, inplace=False, *args, **kwargs):
        """Alias for the `apply_func` method

        args:
            func (object): Function to map across image collection. Function must accept ee.Image as first argument
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False
            **kwargs: arbitrary keyword to pass to `func`

        returns:
            Dataset | None: copy of class with results from `func` as image within collection property
        """

        return self.apply_func(func, inplace, *args, *kwargs)

    def filter(self, filter, inplace=False):
        """Wrapper method for applying a filter to a datset collection.

        args:
            filter (ee.Filter): an `ee.Filter` object to apply to the dataset collection.

        returns:
            Dataset | None: returns the dataset with the filtered collection.
        """

        filtered = self.collection.filter(filter)

        return self._inplace_wrapper(filtered, inplace)

    def select(self, *args, inplace=False):
        """Wrapper method for selecting bands from dataset collection

        args:
            *args: arbitrary arguments to pass to the `ee.ImageCollection.select()` method
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns the dataset with the collection with selected bands.

        """
        selected = self.collection.select(*args)

        return self._inplace_wrapper(selected, inplace)

    def clip_to_region(self, inplace=False):
        """Clips all of the images to the geographic extent defined by region.
        Useful for setting geometries on unbounded imagery in collection (e.g. MODIS or VIIRS imagery)

        args:
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset with imagery clipped to self.region or none depending on inplace
        """

        @decorators.keep_attrs
        def clip(img):
            """Closure function to perform the clipping while carrying metadata"""
            return ee.Image(img.clip(self.region))

        out_coll = self.collection.map(clip)

        return self._inplace_wrapper(out_coll, inplace)

    def merge(self, dataset, inplace=False):
        """Merge the collection of two datasets into one where self.collection will contain imagery from self and dataset arg.
        Results will be sorted by time

        args:
            dataset (Dataset): dataset object to merge
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset where collection is merged imagery or none depending on inplace

        """
        merged = self.collection.merge(dataset.collection).sort("system:time_start")

        return self._inplace_wrapper(merged, inplace)

    def join(self, dataset, inplace=False):
        """Performs spatiotemporal join between self.collection and dataset.collection.
        Result will be a dataset where the collection is colocated imagery in space and time

        args:
            dataset (Dataset): dataset object to apply join with. Used as right in join operations
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset where collection with joined imagery or none depending on inplace

        """

        def _merge(img):
            """Closure func to take results from the join and combine into one image with overlaping region"""
            join_coll = ee.ImageCollection.fromImages(img.get(key))

            img_geom = img.geometry(100)
            join_geom = join_coll.map(geeutils.get_geoms).union(100).geometry(100)
            overlap = img_geom.intersection(join_geom, 100)

            join_data = join_coll.mosaic()

            return img.addBands(join_data).clip(overlap)

        key = str(dataset.__class__.__name__)

        # get a time and space filter
        filter = ee.Filter.And(
            ee.Filter.maxDifference(
                **{
                    "difference": 1000 * 60 * 60 * 24,  # One day in milliseconds
                    "leftField": "system:time_start",
                    "rightField": "system:time_start",
                }
            ),
            ee.Filter.intersects(
                **{"leftField": ".geo", "rightField": ".geo", "maxError": 100}
            ),
        )
        # apply join on collections and save all results
        joined = ee.ImageCollection(
            ee.Join.saveAll(key).apply(
                primary=self.collection, secondary=dataset.collection, condition=filter
            )
        )

        # map over all filtered imagery, mosaic joined matches, and add bands to imagery
        joined = joined.map(_merge)

        return self._inplace_wrapper(joined, inplace)

    def aggregate_time(
        self,
        dates=None,
        period=1,
        period_unit="day",
        reducer="mean",
        rename=True,
        clip_to_area=False,
        inplace=False,
    ):
        """Aggregates multiple images into one based on time periods and a user defined reducer.
        Useful for mosaicing images from same date or time period.
        Expects the images in this dataset to be homogenous (same band names and types)

        args:
            dates (list[str], optional): list of dates defined as beginning time period of aggregatation. default = None,
                all available uniques dates in collection will be used
            period (int, optional): number of days to advance from dates for aggregation. default = 1
            period_unit (str, optional): time unit to advance period for aggregation. default = "day"
            reducer (str | ee.Reducer, optional): reducer to apply to images for aggregation, accepts string reducer name
                or ee.Reducer opbject, default = "mean"
            clip_to_area (bool): switch to clip imagery that has been merged to the overlaping region of imagery, default=False
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset.collection with aggregated imagery or none depending on inplace

        """

        def _aggregation(d):
            """Closure function to map through days and reduce data within a given time period"""
            t1 = ee.Date(d)
            t2 = t1.advance(period, period_unit)
            img = (
                self.collection.filterDate(t1, t2)
                .reduce(reducer)
                .set("system:time_start", t1.millis())
            )
            if rename:
                img = img.rename(band_names)

            geom = (
                ee.FeatureCollection(
                    self.collection.filterDate(t1, t2).map(geeutils.get_geoms)
                )
                .union(100)
                .geometry(100)
            )
            outimg = ee.Algorithms.If(clip_to_area, img.clip(geom), img)

            return outimg

        if dates is None:
            dates = (
                self.collection.aggregate_array("system:time_start")
                .map(lambda x: ee.Date(x).format("YYYY-MM-dd"))
                .distinct()
            )
        else:
            dates = ee.List(dates)

        band_names = ee.Image(self.collection.first()).bandNames()

        out_coll = ee.ImageCollection.fromImages(dates.map(_aggregation))

        return self._inplace_wrapper(out_coll, inplace)

    @decorators.keep_attrs
    def band_pass_adjustment(self, img):
        """Method to apply linear band transformation to dataset image collection.
        Expects that dataset has properties `self.gain` and `self.bias` set

        args:
            img (ee.Image): image to apply regression on

        """
        # linear regression coefficients for adjustment
        return (
            img.multiply(self.gain)
            .add(self.bias)
            .set("system:time_start", img.get("system:time_start"))
        )

    def pipe(self, steps, inplace=False, keep_attrs=True):
        """Method to pipe imagery within dataset through multiple functions at once.
        Assumes the first argument into piped functions are and ee.Image

        args:
            steps (list | tuple): iterable of functions/steps to apply to imagery.
                list must be in the form of (func,func) or with a tuple of function/keyword ((func,kwargs),func)
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset.collection with piped functions applied

        Example:
            ```python
            s1 = hf.Sentinel1(ee.Geometry.Point(105.03,11.72),"2019-10-03","2019-10-05")
            water = s1.pipe(
                (
                    hf.gamma_map, #apply speckle filter
                    (hf.egde_ostu,{'initial_threshold:-16}) # apply water mapping
                )
            )
            ```
        """

        def _piper(funcs):
            """Closure function to nest list of functions"""
            if len(funcs) > 1:
                one_shotter = funcs[0]
                for func in funcs[1:]:
                    one_shotter = pipe | one_shotter | func

            else:
                one_shotter = funcs[0]

            return one_shotter

        fs = []
        # loop through the steps and create partial funcs is kwargs are provided
        for step in steps:
            try:
                func, kwargs = step

            except TypeError:
                func = step
                kwargs = None

            if kwargs is not None:
                pfunc = partial(func, **kwargs)
            else:
                pfunc = func

            fs.append(pfunc)

        # get the piped function
        if keep_attrs:
            one_shot = decorators.keep_attrs(_piper(fs))
        else:
            one_shot = _piper(fs)

        # apply pipe to each image
        out_coll = self.collection.map(lambda img: one_shot(img))

        return self._inplace_wrapper(out_coll, inplace)


class Sentinel1(Dataset):
    """Class extending dataset for the Sentinel 1 collection
    This Sentinel 1 dataset is in backscatter units
    """

    def __init__(self, *args, asset_id="COPERNICUS/S1_GRD", use_qa=True, **kwargs):
        """Initialize Sentinel1 Dataset class

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Sentinel 1 earth engine collection. default="COPERNICUS/S1_GRD"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Sentinel1, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        self.collection = self.collection.filter(
            ee.Filter.listContains("transmitterReceiverPolarisation", "VH")
        )

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Sentinel1 backscatter based on view angle
        Angle threshold values taken from https://doi.org/10.3390/rs13101954
        """
        angle = img.select("angle")
        angle_mask = angle.lt(45.23993).And(angle.gt(30.63993))
        return img.updateMask(angle_mask)

    def add_orbit_band(self, inplace=False):
        """Method to add orbit band from S1 image metadata
        Useful for determining if pixels are from ascending or descending orbits

        args:
            inplace (bool, optional): define whether to return another dataset object or update inplace. default = False

        returns:
            Dataset | None: returns dataset.collection where imagery has the added bands
        """

        def _add_features(img):
            """Closure function to add features as bands to the images"""
            bounds = img.geometry(100)
            orbit = ee.String(img.get("orbitProperties_pass"))
            orbit_band = ee.Algorithms.If(
                orbit.compareTo("DESCENDING"), ee.Image(1), ee.Image(0)
            )

            extraFeatures = ee.Image(orbit_band).rename("orbit")

            return img.addBands(extraFeatures.clip(bounds))

        return self.apply_func(_add_features, inplace=inplace)

    def to_db(self, inplace=False):
        """Convience method to convert units from power to db"""
        out_coll = self.collection.map(geeutils.power_to_db)

        if inplace:
            self.collection = out_coll
            return
        else:
            outCls = self.copy()
            outCls.collection = out_coll
            return outCls

    def to_power(self, inplace=False):
        """Convience method to convert units from db to power"""
        out_coll = self.collection.map(geeutils.db_to_power)

        if inplace:
            self.collection = out_coll
            return
        else:
            outCls = self.copy()
            outCls.collection = out_coll
            return outCls


class Sentinel1Asc(Sentinel1):
    """Class extending dataset for the Sentinel 1 collection for ascending orbit
    This Sentinel 1 dataset is in power units
    """

    def __init__(self, *args, asset_id="COPERNICUS/S1_GRD", use_qa=True, **kwargs):
        """Initialize Sentinel1 Dataset class

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Sentinel 1 earth engine collection. default="COPERNICUS/S1_GRD"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Sentinel1Asc, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        self.filter(ee.Filter.eq("orbitProperties_pass", "ASCENDING"), inplace=True)

        return


class Sentinel1Desc(Sentinel1):
    """Class extending dataset for the Sentinel 1 collection for descending orbit
    This Sentinel dataset is in power units
    """

    def __init__(self, *args, asset_id="COPERNICUS/S1_GRD", use_qa=True, **kwargs):
        """Initialize Sentinel1 Dataset class

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Sentinel 1 earth engine collection. default="COPERNICUS/S1_GRD"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Sentinel1Desc, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        self.filter(ee.Filter.eq("orbitProperties_pass", "DESCENDING"), inplace=True)

        return


class Viirs(Dataset):
    def __init__(
        self,
        *args,
        asset_id="NOAA/VIIRS/001/VNP09GA",
        use_qa=True,
        apply_band_adjustment=False,
        **kwargs,
    ):
        """Initialize VIIRS Dataset class

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the VIIRS earth engine collection. default="NOAA/VIIRS/001/VNP09GA"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            apply_band_adjustment (bool, optional): boolean switch to apply linear band pass equation to convert values to Landsat8. default=False
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Viirs, self).__init__(*args, asset_id=asset_id, use_qa=use_qa, **kwargs)

        # get the bands and rename to common optical names
        coll = self.collection.select(
            self.BANDREMAP.get("viirs"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # band bass adjustment coefficients taken calculated from https://code.earthengine.google.com/876f53861690e483fb3e3439a3571f27
            # slope coefficients
            self.gain = ee.Image.constant(
                [0.68328, 0.66604, 0.78901, 0.95324, 0.98593, 0.88941]
            )
            # y-intercept coefficients
            self.bias = ee.Image.constant(
                [0.016728, 0.030814, 0.023199, 0.036571, 0.026923, 0.021615]
            ).multiply(10000)
            coll = coll.map(self.band_pass_adjustment)

        self.collection = coll

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for VIIRS VNP09GA dataset"""
        cloudMask = geeutils.extract_bits(
            img.select("QF1"), 2, end=3, new_name="cloud_qa"
        ).lt(1)
        shadowMask = geeutils.extract_bits(
            img.select("QF2"), 3, new_name="shadow_qa"
        ).Not()
        snowMask = geeutils.extract_bits(img.select("QF2"), 5, new_name="snow_qa").Not()
        sensorZenith = img.select("SensorZenith").abs().lt(6000)

        env_mask = cloudMask.And(shadowMask).And(sensorZenith)

        # internal pixel quality masks
        pixal_quality_3 = img.select("QF3")
        pixal_quality_4 = img.select("QF4")
        m2_qual = geeutils.extract_bits(pixal_quality_3, 1, new_name="m2_quality").eq(0)
        m4_qual = geeutils.extract_bits(pixal_quality_3, 3, new_name="m4_quality").eq(0)
        m11_qual = geeutils.extract_bits(pixal_quality_4, 0, new_name="m11_quality").eq(0)
        i1_qual = geeutils.extract_bits(pixal_quality_4, 1, new_name="i1_quality").eq(0)
        i2_qual = geeutils.extract_bits(pixal_quality_4, 2, new_name="i2_quality").eq(0)
        i3_qual = geeutils.extract_bits(pixal_quality_4, 3, new_name="i3_quality").eq(0)

        qual_mask = m2_qual.And(m4_qual).And(m11_qual).And(i1_qual).And(i2_qual).And(i3_qual)

        mask = env_mask.And(qual_mask)

        return geeutils.rescale(img).updateMask(mask)


class Modis(Dataset):
    def __init__(
        self, *args, asset_id="MODIS/006/MOD09GA", use_qa=True, **kwargs
    ):
        """Initialize MODIS Dataset class
        Can be used with MOD09GA and MYD09GA

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the MODIS earth engine collection. default="MODIS/006/MOD09GA"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Modis, self).__init__(*args, asset_id=asset_id, use_qa=use_qa, **kwargs)

        self.collection = self.collection.select(
            self.BANDREMAP.get("modis"), self.BANDREMAP.get("new")
        )

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for MODIS MXD09GA dataset"""
        # internal env masks
        qa = img.select("state_1km")
        cloudMask = geeutils.extract_bits(qa, 10, end=11, new_name="cloud_qa").lt(1)
        shadowMask = geeutils.extract_bits(qa, 2, new_name="shadow_qa").Not()
        snowMask = geeutils.extract_bits(qa, 12, new_name="snow_qa").Not()
        cloudAdjMask = geeutils.extract_bits(qa, 13, new_name="cloud_adjacency_qa").Not()
        sensorZenith = img.select("SensorZenith").abs().lt(6000)
        env_mask = cloudMask.And(shadowMask).And(snowMask).And(sensorZenith).And(cloudAdjMask)

        # internal pixel quality masks
        pixal_quality = img.select("QC_500m")
        b1_qual = geeutils.extract_bits(pixal_quality, 2, end=5, new_name="b1_quality").eq(0)
        b2_qual = geeutils.extract_bits(pixal_quality, 6, end=9, new_name="b1_quality").eq(0)
        b3_qual = geeutils.extract_bits(pixal_quality, 10, end=13, new_name="b1_quality").eq(0)
        b4_qual = geeutils.extract_bits(pixal_quality, 14, end=17, new_name="b1_quality").eq(0)
        b6_qual = geeutils.extract_bits(pixal_quality, 22, end=25, new_name="b1_quality").eq(0)
        b7_qual = geeutils.extract_bits(pixal_quality, 26, end=29, new_name="b1_quality").eq(0)

        qual_mask = b1_qual.And(b2_qual).And(b3_qual).And(b4_qual).And(b6_qual).And(b7_qual)

        mask = env_mask.And(qual_mask)
        return geeutils.rescale(img).updateMask(mask)

class Landsat9(Dataset):
    def __init__(
        self,
        *args,
        asset_id="LANDSAT/LC09/C02/T1_L2",
        use_qa=True,
        **kwargs,
    ):
        """Initialize Landsat9 Dataset class
        Can theoretically be useds with any Landsat collection 2 surface reflectance collection (e.g. LANDSAT/LT05/C01/T1_SR)

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Landsat earth engine collection. default="LANDSAT/LC08/C01/T1_SR"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Landsat9, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        self.collection = self.collection.select(
            self.BANDREMAP.get("landsat8"), self.BANDREMAP.get("new")
        )

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Landsat9 surface reflectance dataset"""
        qa_band = img.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = img.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return geeutils.rescale(img, scale = 0.0000275, offset = -0.2).updateMask(mask)


class Landsat8(Dataset):
    def __init__(
        self,
        *args,
        asset_id="LANDSAT/LC08/C02/T1_L2",
        use_qa=True,
        **kwargs,
    ):
        """Initialize Landsat8 Dataset class
        Can theoretically be useds with any Landsat surface reflectance collection (e.g. LANDSAT/LT05/C01/T1_SR)

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Landsat earth engine collection. default="LANDSAT/LC08/C01/T1_SR"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Landsat8, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        self.collection = self.collection.select(
            self.BANDREMAP.get("landsat8"), self.BANDREMAP.get("new")
        )

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Landsat8 surface reflectance dataset"""
        qa_band = img.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = img.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return geeutils.rescale(img, scale = 0.0000275, offset = -0.2).updateMask(mask)


class Landsat7(Dataset):
    def __init__(
        self,
        *args,
        asset_id="LANDSAT/LE07/C02/T1_L2",
        use_qa=True,
        apply_band_adjustment=False,
        **kwargs,
    ):
        """Initialize Landsat7 Dataset class
        Can theoretically be useds with any Landsat surface reflectance collection (e.g. LANDSAT/LT05/C01/T1_SR)

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Landsat7 earth engine collection. default="LANDSAT/LE07/C01/T1_SR"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            apply_band_adjustment (bool, optional): boolean switch to apply linear band pass equation to convert values to Landsat8. default=False
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Landsat7, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        coll = self.collection.select(
            self.BANDREMAP.get("landsat7"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # band bass adjustment coefficients taken from Roy et al., 2016 http://dx.doi.org/10.1016/j.rse.2015.12.024
            # slope coefficients
            self.gain = ee.Image.constant(
                [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]
            )
            # y-intercept coefficients
            self.bias = ee.Image.constant(
                [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]
            ).multiply(10000)
            coll = coll.map(self.band_pass_adjustment)

        self.collection = coll

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Landsat7 surface reflectance dataset"""
        qa_band = img.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = img.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return geeutils.rescale(img, scale = 0.0000275, offset = -0.2).updateMask(mask)


class Landsat5(Dataset):
    def __init__(
        self,
        *args,
        asset_id="LANDSAT/LT05/C02/T1_L2",
        use_qa=True,
        apply_band_adjustment=False,
        **kwargs,
    ):
        """Initialize Landsat7 Dataset class
        Can theoretically be useds with any Landsat surface reflectance collection (e.g. LANDSAT/LT05/C01/T1_SR)

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Landsat7 earth engine collection. default="LANDSAT/LE07/C01/T1_SR"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            apply_band_adjustment (bool, optional): boolean switch to apply linear band pass equation to convert values to Landsat8. default=False
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Landsat5, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        coll = self.collection.select(
            self.BANDREMAP.get("landsat7"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # TODO: get harmonization coefficients for TM to OLI
            # slope coefficients
            self.gain = ee.Image.constant(
                [
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                ]
            )
            # y-intercept coefficients
            self.bias = ee.Image.constant(
                [
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ]
            ).multiply(10000)
            coll = coll.map(self.band_pass_adjustment)

        self.collection = coll

        return

    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Landsat7 surface reflectance dataset"""
        qa_band = img.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = img.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return geeutils.rescale(img, scale = 0.0000275, offset = -0.2).updateMask(mask)


class Sentinel2(Dataset):
    def __init__(
        self,
        *args,
        asset_id="COPERNICUS/S2_SR_HARMONIZED",
        use_qa=True,
        apply_band_adjustment=False,
        **kwargs,
    ):
        """Initialize Sentinel2 Dataset class

        args:
            *args: positional arguments to pass to `Dataset` (i.e. `region`, `start_time`, `end_time`)
            asset_id (str): asset id of the Sentinel2 earth engine collection. default="COPERNICUS/S2_SR"
            use_qa (bool, optional): boolean to determine to use a private `self.qa()` function. default=True
            apply_band_adjustment (bool, optional): boolean switch to apply linear band pass equation to convert values to Landsat8. default=False
            rescale (bool, optional): boolean switch to convert units from scaled int (0-10000) to float (0-1). If false values will be scaled int. default = False
            **kwargs (optional): addtional arbitrary keywords to pass to `Dataset`
        """
        super(Sentinel2, self).__init__(
            *args, asset_id=asset_id, use_qa=use_qa, **kwargs
        )

        coll = self.collection.select(
            self.BANDREMAP.get("sen2"), self.BANDREMAP.get("new")
        )

        if use_qa:
            coll = coll.filter(ee.Filter.eq('cloud_mask_success', 1))

        if apply_band_adjustment:
            # band bass adjustment coefficients taken HLS project https://hls.gsfc.nasa.gov/algorithms/bandpass-adjustment/
            # slope coefficients
            self.gain = ee.Image.constant(
                [0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003]
            )
            # y-intercept coefficients
            self.bias = ee.Image.constant(
                [-0.00411, -0.00093, 0.00094, -0.0001, -0.0015, -0.0012]
            ).multiply(10000)
            coll = coll.map(self.band_pass_adjustment)

        self.collection = coll

        return


    def deduplicate(self, inplace=False):
        """Deduplicate Sentinel2 images by removing images with the same data take identifier"""
        def _dedupe(id):
            """Helper function to remove duplicate images"""
            id_collection = (
                self.collection
                .filter(
                    ee.Filter.stringStartsWith("DATATAKE_IDENTIFIER", ee.String(id))
                )
                .sort('PROCESSING_BASELINE')
            )
            image = (
                id_collection
                .mosaic()
            )
  
            imgbounds = id_collection.map(lambda x: x.geometry(1e3))
  
            start_time = id_collection.aggregate_array('system:time_start').reduce(ee.Reducer.min())
            end_time = id_collection.aggregate_array('system:time_start').reduce(ee.Reducer.max())
  
            return image.clipToCollection(imgbounds).set(
                {
                    'system:time_start': start_time, 
                    'system:time_end':end_time,
                    'system:index':id
                }
            )
        
        datatakes = (
            self.collection
            .aggregate_histogram("DATATAKE_IDENTIFIER")
            .keys()
            .map(lambda x: ee.String(x).slice(0,-7))
        )

        images = ee.ImageCollection.fromImages(
            datatakes.map(_dedupe)
        )

        return self._inplace_wrapper(images, inplace)


    @decorators.keep_attrs
    def qa(self, img):
        """Custom QA masking method for Sentinel2 surface reflectance dataset"""
        CLD_PRB_THRESH = 40
        NIR_DRK_THRESH = 0.175 * 1e4
        CLD_PRJ_DIST = 3
        BUFFER = 100
        CRS = img.select(0).projection()

        # Get s2cloudless image, subset the probability band.
        cld_prb = ee.Image(
            ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
            .filterDate(self.start_time, self.end_time)
            .filterBounds(self.region)
            .filter(ee.Filter.eq("system:index", img.get("system:index")))
            .first()
        ).select("probability")

        worked = cld_prb.bandNames().length().gt(0)

        # Condition s2cloudless by the probability threshold value.
        is_cloud = cld_prb.gt(CLD_PRB_THRESH)

        # Identify water pixels from the SCL band, invert.
        not_water = img.select("SCL").neq(6)

        # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
        dark_pixels = img.select("B8").lt(NIR_DRK_THRESH).multiply(not_water)

        # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
        shadow_azimuth = ee.Number(90).subtract(
            ee.Number(img.get("MEAN_SOLAR_AZIMUTH_ANGLE"))
        )

        # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
        cld_proj = (
            is_cloud.directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST * 10)
            .reproject(**{"crs": CRS, "scale": 120})
            .select("distance")
            .mask()
        )

        # Identify the intersection of dark pixels with cloud shadow projection.
        is_shadow = cld_proj.multiply(dark_pixels)

        # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
        is_cld_shdw = is_cloud.add(is_shadow).gt(0)

        # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
        # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
        is_cld_shdw = (
            is_cld_shdw.focal_min(2)
            .focal_max(BUFFER * 2 / 20)
            .reproject(**{"crs": CRS, "scale": 60})
            .rename("cloudmask")
        )

        # Subset reflectance bands and update their masks, return the result.
        return geeutils.rescale(img).select("B.*").updateMask(is_cld_shdw.Not()).set({'cloud_mask_success':worked})
