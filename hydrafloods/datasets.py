from __future__ import absolute_import
import os
import ee
import math
import copy
import datetime
from pprint import pformat
from functools import partial
from ee.ee_exception import EEException
from hydrafloods import (
    decorators,
    geeutils,
    thresholding,
    fusion,
    fetch,
    preprocess,
    utils,
    filtering,
)


class Dataset:
    """
    
    """

    def __init__(self, region, start_time, end_time, asset_id="", use_qa=True):
        # TODO: add exceptions to check datatypes
        self.region = region  # dtype = ee.Geometry
        self.start_time = start_time
        self.end_time = end_time
        self.asset_id = asset_id
        self.use_qa = use_qa

        self.BANDREMAP = ee.Dictionary(
            {
                "landsat7": ee.List(["B1", "B2", "B3", "B4", "B5", "B7"]),
                "landsat8": ee.List(["B2", "B3", "B4", "B5", "B6", "B7"]),
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

        imgcollection = (
            ee.ImageCollection(self.asset_id)
            .filterBounds(self.region)
            .filterDate(self.start_time, self.end_time)
        )

        if self.use_qa:
            self.collection = imgcollection.map(self._qa)
        else:
            self.collection = imgcollection

    def __repr__(self):
        if isinstance(self.start_time, datetime.datetime):
            ststr = self.start_time.strftime("%Y-%m-%d")
        else:
            ststr = self.start_time

        if isinstance(self.end_time, datetime.datetime):
            etstr = self.end_time.strftime("%Y-%m-%d")
        else:
            etstr = self.end_time
        objDict = {
            "name": self.__class__.__name__,
            "asset_id": self.asset_id,
            "start_time": ststr,
            "end_time": etstr,
            "region": self.region.coordinates().getInfo(),
        }
        strRepr = pformat(objDict, depth=3)
        return f"HYDRAFloods Collection:\n{strRepr}"

    @property
    def collection(self):
        return self._collection

    @collection.setter
    def collection(self, value):
        self._collection = value
        return

    @property
    def n_images(self):
        return self.collection.size().getInfo()

    @property
    def dates(self):
        eeDates = self.collection.aggregate_array("system:time_start").map(
            lambda x: ee.Date(x).format("YYYY-MM-dd HH:mm:ss.SSS")
        )
        return eeDates.getInfo()

    def copy(self):
        """
        Returns a deep copy of the hydrafloods dataset class
        """
        return copy.deepcopy(self)

    def apply_qa(self):
        if self.use_qa:
            coll = self.collection.map(self._qa)

        self.collection = coll.map(geeutils.addTimeBand)

        return

    def clip_to_region(self, inplace=False):
        """
        Clips all of the images to the geographic extent defined by region
        Useful for setting geometries on unbounded imagery in collection
        """

        @decorators.carry_metadata
        def clip(img):
            return ee.Image(img.clip(self.region))

        clipped = self.collection.map(clip)
        if inplace:
            self.collection = clipped
            return
        else:
            outCls = self.copy()
            outCls.collection = clipped
            return outCls

    def apply_func(self, func, inplace=False, **kwargs):
        """
        Wrapper method to apply a function to all of the image in the dataset
        Makes a copy of the collection and reassigns the image collection propety
        Function must accept an ee.ImageCollection and return an ee.ImageCollection

        Args:
            func: Function to map across image collection

        Keywords:
            **kwargs: arbitrary keyword to pass to func

        Returns:
            outCls: copy of class with results as image collection property
        """
        func = partial(func, **kwargs)
        if inplace:
            self.collection = self.collection.map(func)
            return
        else:
            outCls = self.copy()
            outCls.collection = self.collection.map(func)
            return outCls

    def merge(self, dataset, inplace=False):
        merged = self.collection.merge(dataset.collection).sort("system:time_start")

        if inplace:
            self.collection = merged
            return
        else:
            outCls = self.copy()
            outCls.collection = merged
            return outCls

    def join(self, dataset, inplace=False):
        key = str(dataset.__class__.__name__)
        filter = ee.Filter.And(
            ee.Filter.maxDifference(
                **{
                    "difference": 1000 * 60 * 60 * 24,  # One day in milliseconds
                    "leftField": "system:time_start",
                    "rightField": "system:time_start",
                }
            ),
            ee.Filter.intersects(**{"leftField": ".geo", "rightField": ".geo",}),
        )
        joined = ee.ImageCollection(
            ee.Join.saveAll(key).apply(
                primary=self.collection, secondary=dataset.collection, condition=filter
            )
        )

        joined = joined.map(
            lambda x: x.addBands(
                ee.ImageCollection.fromImages(x.get(key)).mosaic().clip(x.geometry())
            )
        )

        if inplace:
            self.collection = joined
            return
        else:
            outCls = self.copy()
            outCls.collection = joined
            return outCls

    def aggregate_time(
        self, dates=None, period=1, reducer="mean", inplace=False
    ):
        # expects the images in this dataset to be homogenous (same band names and types)
        def _aggregation(d):
            t1 = ee.Date(d)
            t2 = t1.advance(period, "day")
            return (
                self.collection.filterDate(t1, t2)
                .reduce(reducer)
                .rename(band_names)
                .set("system:time_start", t1.millis())
            )

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

        if inplace:
            self.collection = out_coll
            return
        else:
            outCls = self.copy()
            outCls.collection = out_coll
            return outCls


class Sentinel1(Dataset):
    def __init__(self, *args, asset_id="COPERNICUS/S1_GRD", **kwargs):
        super(Sentinel1, self).__init__(*args, asset_id=asset_id, **kwargs)

        self.collection = self.collection.filter(
            ee.Filter.listContains("transmitterReceiverPolarisation", "VH")
        )

        return

    @decorators.carry_metadata
    def _qa(self, img):
        angles = img.select("angle")
        return img.updateMask(angles.lt(45).And(angles.gt(30)))

    def add_fusion_features(self):
        def _add_fusion_features(img):
            bounds = img.geometry(100)
            orbit = ee.String(img.get("orbitProperties_pass"))
            orbit_band = ee.Algorithms.If(
                orbit.compareTo("DESCENDING"), ee.Image(1), ee.Image(0)
            )
            vv = img.select("VV")
            vh = img.select("VH")
            ratio = vv.divide(vh).rename("ratio")
            ndpi = vv.subtract(vh).divide(vv.add(vh)).rename("ndpi")

            extraFeatures = ee.Image.cat(
                [ee.Image(orbit_band).rename("orbit"), ratio, ndpi]
            )

            return img.addBands(extraFeatures.clip(bounds))

        self.collection = self.collection.map(_add_fusion_features)

        return self.copy()


class Atms(Dataset):
    def __init__(self, *args, **kwargs):
        super(Atms, self).__init__(*args, use_qa=False, **kwargs)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        return

    @staticmethod
    def extract(date, region, credentials, outDir="./", gridding_radius=50000):
        files = fetch.atms(
            credentials, startTime=date, endTime=None, region=region, outDir=outDir
        )
        geotiffs = list(map(lambda x: preprocess.atms(x, gridding_radius), files))
        return geotiffs

    @staticmethod
    def load(files, gcsBucket="", eeAsset=""):
        if gcsBucket[-1] != "/":
            gcsBucket += "/"
        if eeAsset[-1] != "/":
            eeAsset += "/"

        for f in files:
            fName = os.path.basename(f)
            utils.push_to_gcs(f, gcsBucket)
            bso = gcsBucket + fName
            tstr = fName.split(".")[3]
            t = "{0}-{1}-{2}:{3}:00".format(tstr[:4], tstr[4:6], tstr[6:11], tstr[11:])
            utils.push_to_gee(bso, eeAsset, properties={"time_start": t})

        return


class Viirs(Dataset):
    def __init__(
        self,
        *args,
        asset_id="NOAA/VIIRS/001/VNP09GA",
        apply_band_adjustment=False,
        **kwargs,
    ):
        super(Viirs, self).__init__(*args, asset_id=asset_id, **kwargs)

        coll = self.collection.select(
            self.BANDREMAP.get("viirs"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # band bass adjustment coefficients taken from Roy et al., 2016 http://dx.doi.org/10.1016/j.rse.2015.12.024
            self.gain = ee.Image.constant(
                [0.68328, 0.66604, 0.789009, 0.953243, 0.985927, 0.889414]
            )
            self.bias = ee.Image.constant(
                [0.016728, 0.03081, 0.02319, 0.03657, 0.026923, 0.021615]
            ).multiply(10000)
            coll = coll.map(self._band_pass_adjustment)

        self.collection = coll.map(geeutils.add_indices)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        cloudMask = geeutils.extract_bits(
            img.select("QF1"), 2, end=3, new_name="cloud_qa"
        ).lt(1)
        shadowMask = geeutils.extract_bits(
            img.select("QF2"), 3, new_name="shadow_qa"
        ).Not()
        snowMask = geeutils.extract_bits(img.select("QF2"), 5, new_name="snow_qa").Not()
        sensorZenith = img.select("SensorZenith").abs().lt(6000)

        mask = cloudMask.And(shadowMask).And(sensorZenith)
        return img.updateMask(mask)

    @decorators.carry_metadata
    def _band_pass_adjustment(self, img):
        # linear regression coefficients for adjustment
        return (
            img.multiply(self.gain)
            .add(self.bias)
            .set("system:time_start", img.get("system:time_start"))
        )

    @staticmethod
    def extract(date, region, outdir="./", creds=None):
        files = fetch.viirs(
            creds, startTime=date, endTime=None, region=region, outDir=outdir
        )

        return

    @staticmethod
    def load(self, files, gcsBucket="", eeAsset=""):

        return


class Modis(Dataset):
    def __init__(self, *args, asset_id="MODIS/006/MOD09GA", **kwargs):
        super(Modis, self).__init__(*args, asset_id=asset_id, **kwargs)

        self.collection = self.collection.select(
            self.BANDREMAP.get("modis"), self.BANDREMAP.get("new")
        ).map(geeutils.add_indices)

        self.clip_to_region(inplace=True)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        qa = img.select("state_1km")
        cloudMask = geeutils.extract_bits(qa, 10, end=11, new_name="cloud_qa").lt(1)
        shadowMask = geeutils.extract_bits(qa, 2, new_name="shadow_qa").Not()
        snowMask = geeutils.extract_bits(qa, 12, new_name="snow_qa").Not()
        sensorZenith = img.select("SensorZenith").abs().lt(6000)
        mask = cloudMask.And(shadowMask).And(snowMask).And(sensorZenith)
        return img.updateMask(mask)

    def extract(self, date, region, outdir="./", creds=None):

        return

    def load(self, files, gcsBucket="", eeAsset=""):

        return


class Landsat8(Dataset):
    def __init__(self, *args, asset_id="LANDSAT/LC08/C01/T1_SR", **kwargs):
        super(Landsat8, self).__init__(*args, asset_id=asset_id, **kwargs)

        self.collection = self.collection.select(
            self.BANDREMAP.get("landsat8"), self.BANDREMAP.get("new")
        ).map(geeutils.add_indices)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        qa_band = img.select("pixel_qa")
        qaCloud = geeutils.extract_bits(qa_band, start=5, new_name="cloud_mask").eq(0)
        qaShadow = geeutils.extract_bits(qa_band, start=3, new_name="shadow_mask").eq(0)
        qaSnow = geeutils.extract_bits(qa_band, start=4, new_name="snow_mask").eq(0)
        mask = qaCloud.And(qaShadow)  # .And(qaSnow)
        return img.updateMask(mask)


class Landsat7(Dataset):
    def __init__(
        self,
        *args,
        asset_id="LANDSAT/LE07/C01/T1_SR",
        apply_band_adjustment=False,
        **kwargs,
    ):
        super(Landsat7, self).__init__(*args, asset_id=asset_id, **kwargs)

        coll = self.collection.select(
            self.BANDREMAP.get("landsat7"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # band bass adjustment coefficients taken from Roy et al., 2016 http://dx.doi.org/10.1016/j.rse.2015.12.024
            self.gain = ee.Image.constant(
                [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]
            )
            self.bias = ee.Image.constant(
                [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]
            ).multiply(10000)
            coll = coll.map(self._band_pass_adjustment)

        self.collection = coll.map(geeutils.add_indices)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        qa_band = img.select("pixel_qa")
        qaCloud = geeutils.extract_bits(qa_band, start=5, new_name="cloud_mask").eq(0)
        qaShadow = geeutils.extract_bits(qa_band, start=3, new_name="shadow_mask").eq(0)
        qaSnow = geeutils.extract_bits(qa_band, start=4, new_name="snow_mask").eq(0)
        mask = qaCloud.And(qaShadow)  # .And(qaSnow)
        return img.updateMask(mask)

    @decorators.carry_metadata
    def _band_pass_adjustment(self, img):
        # linear regression coefficients for adjustment
        return (
            img.multiply(self.gain)
            .add(self.bias)
            .set("system:time_start", img.get("system:time_start"))
        )


class Sentinel2(Dataset):
    def __init__(
        self, *args, asset_id="COPERNICUS/S2_SR", apply_band_adjustment=False, **kwargs
    ):
        super(Sentinel2, self).__init__(*args, asset_id=asset_id, **kwargs)

        coll = self.collection.select(
            self.BANDREMAP.get("sen2"), self.BANDREMAP.get("new")
        )

        if apply_band_adjustment:
            # band bass adjustment coefficients taken HLS project https://hls.gsfc.nasa.gov/algorithms/bandpass-adjustment/
            self.gain = ee.Image.constant(
                [0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003]
            )
            self.bias = ee.Image.constant(
                [-0.00411, -0.00093, 0.00094, -0.0001, -0.0015, -0.0012]
            ).multiply(10000)
            coll = coll.map(self._band_pass_adjustment)

        self.collection = coll.map(geeutils.add_indices)

        return

    @decorators.carry_metadata
    def _qa(self, img):
        sclImg = img.select("SCL")  # Scene Classification Map
        mask = sclImg.gte(4).And(sclImg.lte(6))
        return img.updateMask(mask)

    @decorators.carry_metadata
    def _band_pass_adjustment(self, img):
        # linear regression coefficients for adjustment
        return (
            img.multiply(self.gain)
            .add(self.bias)
            .set("system:time_start", img.get("system:time_start"))
        )

