# daily surface water - fusion process

import os
import ee
import time
from ee.ee_exception import EEException
import gcsfs
import logging
import datetime
from functools import partial
import warnings
import hydrafloods as hf
from hydrafloods import (
    datasets,
    timeseries,
    ml,
    utils,
    geeutils,
    indices,
    thresholding,
    decorators,
    corrections,
)


def _get_output_names(band_names: ee.List, reducer_names: ee.List) -> ee.List:
    def name_outer_loop(prefix):
        """Closure function to loop over output band name prefixes"""

        def name_inner_loop(suffix):
            """Closure function to loop over output band name suffixes and combine with prefixes"""
            return ee.String(prefix).cat(ee.String("_")).cat(ee.String(suffix))

        # loop over suffixes
        return reducer_names.map(name_inner_loop)

    return band_names.map(name_outer_loop).flatten()


def _calc_sar_anomalies(
    years, s1_ds, sar_bands=["VH", "VV", "ndpi", "nvhi", "nvvi", "ratio"]
):
    def calc_anomalies_yr(yr):
        """Closure function to loop over each year and calculate stats then"""
        # @decorators.keep_attrs
        def calc_anomalies_img(img):
            """Closure function to calculate the difference of img from yearly stats"""
            anomalies = (
                ee.ImageCollection.fromImages(
                    ee.List(
                        # loop over each band prefix to make sure we calculate difference from correct combo of bands
                        sar_bands.map(
                            # calculate image difference from yearly statistics
                            lambda x: img.select(ee.String(x)).subtract(
                                stats_yr.select(ee.String(x).cat("_.*"))
                            )
                        )
                    )
                )
                .toBands()  # convert imagecollection to image
                .float()  # cast to float
                .rename(out_band_names)  # rename the output bands
            )

            return img.addBands(anomalies)

        # get the time range for the year
        t1 = ee.Date.fromYMD(yr, 1, 1)
        t2 = t1.advance(1, "year")
        # select the imagery to reduce for the year
        s1_yr = s1_ds.collection.filterDate(t1, t2)
        # apply reduction
        stats_yr = s1_yr.reduce(all_stats_reducer, 8)
        # loop over each image within year and calculate anomalies stats
        anomalies = s1_yr.map(calc_anomalies_img)

        return anomalies

    # cast band name list to ee.List and sort
    # need to sort to make sure order matches output
    sar_bands = ee.List(sar_bands).sort()

    # get all of the statistics in one reducer so it is more efficient
    all_stats_reducer = (
        ee.Reducer.mean()
        .combine(ee.Reducer.stdDev(), None, True)
        .combine(ee.Reducer.percentile([5, 95]), None, True)
    )

    # get the names of the reducer suffixes
    reducer_suffixes = all_stats_reducer.getOutputs().sort()

    # create a new list of expected band names based on
    out_band_names = _get_output_names(sar_bands, reducer_suffixes)

    # select the sar bands we want to process to trim any unneccessary computations
    s1_ds.collection = s1_ds.collection.select(sar_bands)

    # apply the anomaly process
    # outputs is weird formation List[<ImageCollection>, ..., <ImageCollection>]
    # need to convert to FeatureCollection -> flatten -> then ImageCollection
    s1_anomalies = ee.ImageCollection(
        ee.FeatureCollection(years.map(calc_anomalies_yr)).flatten()
    )

    # output is an ImageCollection...be careful if using with hf.Dataset
    return s1_anomalies


def _fuse_dataset(
    region,
    start_time,
    end_time,
    fusion_model,
    scaling_dict,
    feature_names,
    target_band="mndwi",
    use_viirs=False,
    use_modis=False,
):
    ds_kwargs = dict(region=region, start_time=start_time, end_time=end_time)
    dsa_kwargs = {**ds_kwargs, **{"apply_band_adjustment": True}}

    lc8 = datasets.Landsat8(**ds_kwargs)
    le7 = datasets.Landsat7(**dsa_kwargs)
    s2 = datasets.Sentinel2(**dsa_kwargs)

    optical = lc8.merge(le7).merge(s2)

    if use_viirs:
        viirs = datasets.Viirs(**dsa_kwargs)
        optical = optical.merge(viirs)

    if use_modis:
        modis = datasets.Modis(**ds_kwargs)
        optical = optical.merge(modis)

    optical = optical.apply_func(geeutils.add_indices, indices=[target_band])

    s1 = datasets.Sentinel1(**ds_kwargs)

    years = (
        s1.collection.aggregate_array("system:time_start")
        .map(lambda x: ee.Date(x).get("year"))
        .distinct()
    )

    dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation")

    s1_proc = (
        (corrections.slope_correction, dict(elevation=dem, buffer=30)),
        (geeutils.add_indices, dict(indices=["vv_vh_ratio", "ndpi", "nvvi", "nvhi"])),
        #     (
        #         ml.standard_image_scaling,
        #         dict(scaling_dict=scaling_dict, feature_names=feature_names),
        #     ),
        #     lambda x: x.classify(fusion_model, target_band),
    )
    s1 = s1.pipe(s1_proc)

    s1_anomalies = _calc_sar_anomalies(years, s1)

    s1.collection = s1_anomalies.map(
        decorators.keep_attrs(lambda x: x.classify(fusion_model, target_band))
    )

    fused_ds = optical.merge(s1)

    fused_ds.collection = (
        fused_ds.collection.select(target_band)
        .cast({target_band: "float"}, [target_band])
        .sort("system:time_start")
    )

    return fused_ds


def export_fusion_samples(
    region,
    start_time,
    end_time,
    output_asset_path,
    stratify_samples=True,
    sample_scale=30,
    n_samples=25,
    seed=0,
):
    """First step of the daily surface water fusion process.
    This procedure samples values from coincident optical and SAR data
    so that we can use ML for data fusion. This will calculate MNDWI, NWI,
    AEWInsh, and AEWIsh optical water indices and a few indices from SAR imagery
    (VV/VH, NDPI, NVVI, NVHI) to predict a water index.

    args:
        region (ee.Geometry): geographic region to look for coincident data and sample from
        start_time (str | datetime.datetime): start time used to look for coincident data
        end_time (str | datetime.datetime): end time used to look for coincident data
        output_asset_path (str): Earth Engine asset id to save sampled values too
        stratify_samples (bool, optional): boolean keyword to specify for sampling data stratified
            by a combination of the MODIS land cover and JRC surface water occurrence. If False,
            then a random sampling wil be used. default = False
        sample_scale (float, optional): resolution in meters to sample data at
        n_samples (int, optional): number of samples to collect per coincident image pair. If
            stratified_samples == True, this value be be samples per class. default = 25
        seed (int,optional): random number generator seed, used for setting random sampling. default = 0

    """

    dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation")

    optical_water_indices = ["mndwi", "nwi", "aewish", "aewinsh"]

    ds_kwargs = dict(
        region=region, start_time=start_time, end_time=end_time, rescale=True
    )
    dsa_kwargs = {**ds_kwargs, **{"apply_band_adjustment": True}}

    lc8 = datasets.Landsat8(**ds_kwargs)
    le7 = datasets.Landsat7(**dsa_kwargs)
    s2 = datasets.Sentinel2(**dsa_kwargs)

    _ = ds_kwargs.pop("rescale")

    s1a = datasets.Sentinel1Asc(**ds_kwargs)
    s1d = datasets.Sentinel1Desc(**ds_kwargs)

    years = (
        s1a.collection.aggregate_array("system:time_start")
        .map(lambda x: ee.Date(x).get("year"))
        .distinct()
    )

    sar_proc = (
        (
            corrections.slope_correction,
            dict(
                elevation=dem,
                buffer=50,
            ),
        ),
        hf.gamma_map,
        (geeutils.add_indices, dict(indices=["vv_vh_ratio", "ndpi", "nvvi", "nvhi"])),
    )

    s1a.pipe(sar_proc, inplace=True)
    s1d.pipe(sar_proc, inplace=True)

    s1a_anomalies = _calc_sar_anomalies(years, s1a)
    s1d_anomalies = _calc_sar_anomalies(years, s1d)

    s1a.collection = s1a_anomalies
    s1d.collection = s1d_anomalies

    s1 = s1a.merge(s1d)

    optical = lc8.merge(s2).merge(le7)

    optical = optical.apply_func(geeutils.add_indices, indices=optical_water_indices)

    ds = optical.join(s1)

    sample_region = (
        ds.collection.map(geeutils.get_geoms)
        .union(maxError=1000)
        .geometry(maxError=1000)
    ).intersection(region, maxError=1000)

    img_list = ds.collection.toList(ds.collection.size())

    output_features = ee.FeatureCollection([])
    aggregate_samples = 0

    if stratify_samples:

        class_band = "strata"
        interval = 20

        water_freq = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select("occurrence")

        class_intervals = ee.List.sequence(0, 80, interval)
        logging.info(f"Water intervals: {class_intervals.getInfo()}")
        n_water_classes = class_intervals.size()
        water_classes = ee.List.sequence(1, n_water_classes)
        logging.info(f"Water Classes: {water_classes.getInfo()}")

        water_img = (
            ee.ImageCollection.fromImages(
                class_intervals.map(lambda x: water_freq.gte(ee.Number(x)))
            )
            .reduce(ee.Reducer.sum())
            .uint8()
            .rename(class_band)
        )

        # class_band = "landcover"
        igbp_classes = ee.List(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        )
        ipcc_classes = ee.List([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 6, 3, 5, 3, 4, 4, 6])

        lc_img = (
            ee.ImageCollection("MODIS/006/MCD12Q1")
            .limit(1, "system:time_start", True)
            .first()
            .remap(igbp_classes, ipcc_classes)
            .rename(class_band)
        ).add(water_classes.size())
        lc_classes = (
            ipcc_classes.distinct()
            .map(lambda x: ee.Number(x).add(water_classes.size()))
            .sort()
        )
        logging.info(f"LC Classes: {lc_classes.getInfo()}")

        final_strata_img = water_img.unmask(lc_img).rename(class_band)

        half = ee.Number(n_samples).multiply(n_water_classes.subtract(1))
        n_lc = half.divide(lc_classes.size()).round()

        all_classes = (
            water_classes.slice(1)
            .cat(lc_classes)
            .map(lambda x: ee.Number(x).subtract(1))
        )
        n_water_samples = ee.List.repeat(n_samples, n_water_classes)
        n_lc_samples = ee.List.repeat(n_lc, lc_classes.size())
        logging.info(f"n Water Samples {n_water_samples.getInfo()}")
        logging.info(f"n LC Samples {n_lc_samples.getInfo()}")

        base_samples = water_img.select(class_band).stratifiedSample(
            region=sample_region,
            numPoints=n_samples,
            classBand=class_band,
            scale=sample_scale,
            seed=seed,
            classValues=water_classes,
            classPoints=n_water_samples,
            tileScale=16,
            geometries=True,
        )
    else:
        base_samples = ee.FeatureCollection.randomPoints(sample_region, n_samples, seed)

    def sample_img(img):
        geom = img.geometry()
        date = img.date()
        week = date.get("week")
        year = date.get("year")
        new_seed = week.add(year).add(seed)

        lc_samples = lc_img.select(class_band).stratifiedSample(
            region=sample_region,
            numPoints=n_samples,
            classBand=class_band,
            scale=sample_scale,
            seed=new_seed,
            classValues=lc_classes,
            classPoints=n_lc_samples,
            tileScale=16,
            geometries=True,
        )
        samples = (
            base_samples.merge(lc_samples)
            .filterBounds(geom)
            .randomColumn("random", seed)
        )

        features = img.sampleRegions(
            samples, scale=sample_scale, tileScale=16, geometries=True
        )

        features = features.map(lambda x: ee.Feature(x).set("timestamp", date.millis()))

        return features

    sample_features = ds.collection.map(sample_img).flatten()

    output_features = ee.FeatureCollection(
        sample_features.aggregate_array(class_band)
        .distinct()
        .map(
            lambda x: sample_features.filter(ee.Filter.eq(class_band, x)).randomColumn()
        )
    ).flatten()

    now = datetime.datetime.now()
    time_id = now.strftime("%Y%m%d%H%M%s")

    task = ee.batch.Export.table.toAsset(
        collection=output_features,
        assetId=output_asset_path,
        description=f"hydrafloods_fusion_samples_export_{time_id}",
    )
    task.start()
    logging.info(f"Started export task for {output_asset_path}")

    return


def export_surface_water_harmonics(
    region,
    start_time,
    end_time,
    output_asset_path,
    n_cycles=2,
    feature_names=None,
    label=None,
    fusion_samples=None,
    tile=False,
    tile_size=1.0,
    output_scale=30,
):
    """Second step of the daily surface water fusion process.
    This procedure uses samples from `export_fusion_samples` to build a
    random forest model to predict a water index from SAR imagery. This a
    time series of optical-SAR fused data is used to calculate a harmonic
    model for long-term surface water trend and is exported to an Earth Engine
    asset.

    args:
        region (ee.Geometry): geographic region to look for coincident data and sample from
        start_time (str | datetime.datetime): start time used to look for coincident data
        end_time (str | datetime.datetime): end time used to look for coincident data
        output_asset_path (str): Earth Engine asset id to save harmonic model weights to as image.
            If tile==True, then output_asset_path much be a precreated ImageCollection asset
        n_cycles (int, optional): number of interannual cycles to model. default = 2
        feature_names (list[str],):  names of feature columns used to calculate `label` from
        label (str): name of feature column to predict using `feature_names`
        fusion_samples (str): Earth Engine FeatureCollection asset id of samples to get a data
            fusion model from. Should be the asset output from `export_fusion_samples`
        tile (bool, optional): boolean keyword to tile exports. If false will try to calculate
            harmonic weights as image. If true, it will tile area and recusively call to export
            smaller areas. If true then expects that `output_asset_path` is an ImageCollection.
            default = False
        tile_size (float, optional): resolution in decimal degrees to create tiles over region
            for smaller exports. Only used if tile==True. default = 1.0
        output_scale (float, optional): output resolution of harmonic weight image. default = 30
    """

    if tile:
        land_area = (
            ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
            .filterBounds(region)
            .geometry(1000)
            .buffer(2500, maxError=1000)
        )
        grid = geeutils.tile_region(
            region, intersect_geom=land_area, grid_size=tile_size
        )

        n = grid.size().getInfo()
        grid_list = grid.toList(n)

        for i in range(n):
            grid_tile = ee.Feature(grid_list.get(i)).geometry()
            if output_asset_path is not None:
                output_tile_path = output_asset_path + f"harmonics_t{i:05d}"

            export_surface_water_harmonics(
                region=grid_tile,
                start_time=start_time,
                end_time=end_time,
                n_cycles=n_cycles,
                feature_names=feature_names,
                label=label,
                fusion_samples=fusion_samples,
                output_asset_path=output_tile_path,
                tile=False,
                tile_size=tile_size,
                output_scale=output_scale,
            )

    else:
        if fusion_samples is not None:
            fusion_model, scaling_dict = ml.random_forest_ee(
                30,
                fusion_samples,
                feature_names,
                label,
                scaling=None,
                mode="regression",
            )
        else:
            raise ValueError(
                "Either 'fusion_samples' or 'fusion_model_path' needs to be defined to run fusion process"
            )

        ds = _fuse_dataset(
            region,
            start_time,
            end_time,
            fusion_model,
            scaling_dict,
            feature_names,
            target_band=label,
        )

        now = datetime.datetime.now()
        time_id = now.strftime("%Y%m%d%H%M%s")
        time_str = now.strftime("%Y-%m-%d %H:%M:%s")

        scale_factor = 0.0001

        # create metadata dict
        metadata = ee.Dictionary(
            {
                "hf_version": hf.__version__,
                "scale_factor": scale_factor,
                "fit_time_start": start_time,
                "fit_time_end": end_time,
                "execution_time": time_str,
            }
        )

        harmonic_coefs = timeseries.fit_harmonic_trend(
            ds, dependent=label, n_cycles=n_cycles, output_err=True
        )
        harmonic_coefs = harmonic_coefs.divide(scale_factor).int32().set(metadata)

        if output_asset_path is not None:
            geeutils.export_image(
                harmonic_coefs,
                region,
                output_asset_path,
                description=f"hydrafloods_harmonic_coefficient_export_{time_id}",
                scale=output_scale,
                crs="EPSG:4326",
            )
        else:
            raise ValueError(
                "'output_asset_path' needs to be defined to run fusion export process"
            )

    return


def export_fusion_product(
    region,
    target_date,
    harmonic_image=None,
    harmonic_collection=None,
    feature_names=None,
    label=None,
    look_back=30,
    lag=4,
    n_cycles=2,
    fusion_samples=None,
    output_asset_path=None,
    output_bucket_path=None,
    tile=False,
    tile_size=1.0,
    tile_buffer=100000,
    output_scale=30,
):
    def get_residuals(i):
        """Closure function to calculate residuals of harmonic water estimate
        compared to observed data.
        """
        i = ee.Number(i)
        t_diff = (
            ee.Number(i).multiply(-1).subtract(lag)
        )  # calc how many days to adjust ini date
        new_date = target_date.advance(t_diff, "day")  # calculate new date

        corr_img = (
            ds.collection.select(label)
            .filterDate(new_date, new_date.advance(1, "day"))
            .median()
        )  # .select(f"^{label}.*",label)

        time_img = timeseries.get_dummy_img(new_date)

        harmon_pred = (
            timeseries.add_harmonic_coefs(time_img)
            .multiply(harmonic_coefs)
            .reduce("sum")
        )

        harmon_diff = harmon_pred.subtract(corr_img).rename("residual")

        return harmon_diff.set("system:time_start", new_date.millis())

    if tile:
        if tile:
            land_area = (
                ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
                .filterBounds(region)
                .geometry(100)
                .buffer(2500, maxError=100)
            )
            grid = geeutils.tile_region(
                region, intersect_geom=land_area, grid_size=tile_size
            )

            n = grid.size().getInfo()
            grid_list = grid.toList(n)

            for i in range(n):
                if output_asset_path is not None:
                    output_asset_tile = output_asset_path + f"daily_tile{i:05d}"
                else:
                    output_asset_tile = None
                if output_bucket_path is not None:
                    output_bucket_tile = output_bucket_path + f"_tile{i:05d}"
                else:
                    output_bucket_tile = None

                grid_tile = ee.Feature(grid_list.get(i)).geometry()
                export_fusion_product(
                    region=grid_tile,
                    target_date=target_date,
                    harmonic_image=harmonic_image,
                    harmonic_collection=harmonic_collection,
                    feature_names=feature_names,
                    label=label,
                    look_back=look_back,
                    lag=lag,
                    n_cycles=n_cycles,
                    fusion_samples=fusion_samples,
                    output_asset_path=output_asset_path,
                    output_bucket_path=output_bucket_path,
                    tile=tile,
                    tile_size=tile_size,
                    tile_buffer=tile_buffer,
                    output_scale=output_scale,
                )

    else:
        if not isinstance(target_date, ee.Date):
            target_date = ee.Date(target_date)

        end_time = target_date.advance(-(lag - 1), "day")
        start_time = end_time.advance(-look_back, "day")

        if fusion_samples is not None:
            fusion_model, scaling_dict = ml.random_forest_ee(
                30,
                fusion_samples,
                feature_names,
                label,
                scaling=None,
                mode="regression",
            )
        else:
            raise ValueError(
                "'fusion_samples' needs to be defined to run fusion process"
            )

        now = datetime.datetime.now()
        time_id = now.strftime("%Y%m%d%H%M%s")
        time_str = now.strftime("%Y-%m-%d %H:%M:%s")

        if harmonic_image is not None:
            harmonic_coefs = ee.Image(harmonic_image)
            harmonic_coefs = harmonic_coefs.multiply(
                ee.Image(ee.Number(harmonic_coefs.get("scale_factor")))
            )
        elif harmonic_collection is not None:
            harmonic_collection = ee.ImageCollection(harmonic_collection)
            first = ee.Image(harmonic_collection.first())
            harmonic_coefs = harmonic_collection.mosaic().multiply(
                ee.Image(ee.Number(first.get("scale_factor")))
            )
        else:
            raise ValueError(
                "Either 'harmonic_image' or 'harmonic_collection' needs to be defined to run fusion process"
            )

        harmonic_coefs = harmonic_coefs.select("^(c|t|s).*")

        prod_region = region.buffer(tile_buffer, 100)

        ds = _fuse_dataset(
            region,
            start_time,
            end_time,
            fusion_model,
            scaling_dict,
            feature_names,
            target_band=label,
            use_viirs=True,
        )

        dummy_target = timeseries.get_dummy_img(target_date)

        weights = ee.ImageCollection.fromImages(
            ee.List.sequence(0, look_back - 1).map(get_residuals)
        ).sort("system:time_start")

        weights_lr = timeseries.fit_linear_trend(
            weights, dependent="residual", output_err=False
        )

        weights_coefs = weights_lr.select("^(c|t).*")

        lin_pred = (
            dummy_target.multiply(weights_coefs).reduce("sum").rename("residual_est")
        )

        har_pred = (
            timeseries.add_harmonic_coefs(dummy_target, n_cycles=n_cycles)
            .multiply(harmonic_coefs)
            .reduce("sum")
        )

        fused_pred = (har_pred.subtract(lin_pred)).rename("fused_product")

        fused_pred = fused_pred.multiply(10000).int16()

    if output_asset_path is not None:
        # create metadata dict
        metadata = ee.Dictionary(
            {
                "hf_version": hf.__version__,
                "system:time_start": target_date.millis(),
                "system:time_end": target_date.advance(86399, "seconds").millis(),
                "execution_time": time_str,
                "lag": lag,
                "look_back": look_back,
                "scale_factor": 0.0001,
            }
        )

        geeutils.export_image(
            fused_pred.set(metadata.combine({"product": "fusion"})),
            region,
            output_asset_path + "_fusion",
            description=f"hydrafloods_fusion_ee_export_{time_id}",
            scale=output_scale,
            crs="EPSG:4326",
        )

    elif output_bucket_path is not None:
        export_region = region.bounds(maxError=100).getInfo()["coordinates"]
        bucket_path, ext = os.path.splitext(output_bucket_path)
        fcomponents = bucket_path.split("/")
        bucket = fcomponents[2]
        fpath = fcomponents[3:-1]

        # TODO: remove extension from string formulation
        f_water = "/".join(fpath + [fcomponents[-1] + "_water" + ext])
        f_fusion = "/".join(fpath + [fcomponents[-1] + "_fusion" + ext])

        fusion_task = ee.batch.Export.image.toCloudStorage(
            image=fused_pred,
            description=f"hydrafloods_fusion_gcp_export_{time_id}",
            bucket=bucket,
            fileNamePrefix=f_fusion,
            region=export_region,
            scale=output_scale,
            crs="EPSG:4326",
            maxPixels=1e13,
            fileFormat="GeoTIFF",
            formatOptions={"cloudOptimized": True},
        )
        fusion_task.start()

    else:
        raise ValueError(
            "Either 'output_asset_path' or 'output_bucket_path' needs to be defined to run fusion export process"
        )

    return


def export_daily_surface_water(
    region,
    target_date,
    harmonic_image=None,
    harmonic_collection=None,
    feature_names=None,
    label=None,
    look_back=30,
    lag=4,
    n_cycles=2,
    include_confidence=False,
    include_flood=False,
    fusion_samples=None,
    output_asset_path=None,
    output_bucket_path=None,
    initial_threshold=0.1,
    thresh_no_data=None,
    tile=False,
    tile_size=1.0,
    tile_buffer=100000,
    output_scale=30,
):
    """Last and repeated step of the daily surface water fusion process.
    This procedure uses the results from `export_fusion_samples` and
    `export_surface_water_harmonics` to build a random forest model to predict
    a water index from SAR imagery and predict water using the harmonic model.
    This process will correct the harmonic estimate using observed data and export
    the resulting imagery.

    args:
        region (ee.Geometry): geographic region to look for coincident data and sample from
        target_date (str | datetime.datetime): date to estimate surface water extent for
        harmonic_image (str, optional): Earth Engine Image asset id of the harmonic model weights
            exported by `export_surface_water_harmonics`. If left as None then `harmonic_collection`
            must be defined. default = None
        harmonic_collection (str, optional): Earth Engine ImageCollection asset id of the harmonic
            model weights from tile `export_surface_water_harmonics`. If left as None then
            `harmonic_image` must be defined. default = None
        feature_names (list[str],):  names of feature columns used to calculate `label` from
        label (str): name of feature column to predict using `feature_names`
        look_back (int,optional): number of days used to estimate short-term trend in water. default = 30
        lag (int, optional): number of days after `target_date` to begin `look_back`. default=4
        n_cycles (int, optional): number of interannual cycles to model. default = 2
        include_confidence (bool, optional): boolean keyword to specify if a confidence band will
            be exported with surface water image. If True then confidence will be calculated. default = False
        include_flood (bool, optional): boolean keyword to specify if a flood band will
            be exported with surface water image. If True then flood will be calculated based on JRC
             permanent water data. default = False
        export_fusion (bool, optional): boolean keyword to specify if the fusion image used to calculate
            water should be exported as a seperated task. If True then run fusion export task. default = False
        fusion_samples (str): Earth Engine FeatureCollection asset id of samples to get a data
            fusion model from. Should be the asset output from `export_fusion_samples`
        output_asset_path (str): Earth Engine asset id to save estimate water and fusion results to as image.
            If tile==True, then output_asset_path much be a precreated ImageCollection asset. If left
            as None then `output_bucket_path` must be specified. default = None
        output_bucket_path (str): GCP cloud bucket path to save estimate water and fusion results to
            cloud optimized geotiffs. If tile==True, then multiple file will be created. If left
            as None then `output_asset_path` must be specified. default = None
        initial_threshold (float, optional): initial threshold value used in `edge_otsu` thresholding
            algorithm to segment water from fusion image. default = 0.1
        tile (bool, optional): boolean keyword to tile exports. If false will try to calculate
            harmonic weights as image. If true, it will tile area and recusively call to export
            smaller areas. If true then expects that `output_asset_path` is an ImageCollection.
            default = False
        tile_size (float, optional): resolution in decimal degrees to create tiles over region
            for smaller exports. Only used if tile==True. default = 1.0
        tile_buffer (float,optional): buffer size in meters to buffer tiles to calculate threshold. This
            is used to ensure running tiled exports produces consistent results at tile seams. default = 100000
        output_scale (float, optional): output resolution of harmonic weight image. default = 30

    raises:
        ValueError: if `fusion_samples` is None
        ValueError: if both`harmonic_image` and `harmonic_collection` is None
        ValueError: if both 'output_asset_path' and 'output_bucket_path' is None
    """

    def get_residuals(i):
        """Closure function to calculate residuals of harmonic water estimate
        compared to observed data.
        """
        i = ee.Number(i)
        t_diff = (
            ee.Number(i).multiply(-1).subtract(lag)
        )  # calc how many days to adjust ini date
        new_date = target_date.advance(t_diff, "day")  # calculate new date

        corr_img = (
            ds.collection.select(label)
            .filterDate(new_date, new_date.advance(1, "day"))
            .median()
        )

        time_img = timeseries.get_dummy_img(new_date)

        harmon_pred = (
            timeseries.add_harmonic_coefs(time_img)
            .multiply(harmonic_coefs)
            .reduce("sum")
        )

        harmon_diff = harmon_pred.subtract(corr_img).rename("residual")

        return harmon_diff.set("system:time_start", new_date.millis())

    def calc_confidence(i):
        """Closure function to calculate confidence in water estimate using
        monte carlo methods and simulating errors in long- and short-term water
        dynamics
        """
        i = ee.Number(i)
        # uniform sampling of std dev at 95% confidence interval
        long_term_seed = i.add(500)
        short_term_seed = i.add(1000)
        long_term_random = ee.Image.random(long_term_seed).multiply(3.92).subtract(1.96)
        short_term_random = (
            ee.Image.random(short_term_seed).multiply(3.92).subtract(1.96)
        )

        lin_sim = lin_pred.add(short_term_random.multiply(linCi))
        har_sim = har_pred.add(long_term_random.multiply(harCi))

        sim_pred = har_sim.subtract(lin_sim)
        # random_water = thresholding.bmax_otsu(random_combination,invert=True)
        # naive estimate of water (>0)
        return sim_pred.gt(ci_threshold).uint8()

    if tile:
        if tile:
            land_area = (
                ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
                .filterBounds(region)
                .geometry(100)
                .buffer(2500, maxError=100)
            )
            grid = geeutils.tile_region(
                region, intersect_geom=land_area, grid_size=tile_size
            )

            n = grid.size().getInfo()
            grid_list = grid.toList(n)

            for i in range(n):
                if output_asset_path is not None:
                    output_asset_tile = output_asset_path + f"daily_tile{i:05d}"
                else:
                    output_asset_tile = None
                if output_bucket_path is not None:
                    output_bucket_tile = output_bucket_path + f"_tile{i:05d}"
                else:
                    output_bucket_tile = None

                grid_tile = ee.Feature(grid_list.get(i)).geometry()
                export_daily_surface_water(
                    region=grid_tile,
                    target_date=target_date,
                    harmonic_image=harmonic_image,
                    harmonic_collection=harmonic_collection,
                    feature_names=feature_names,
                    label=label,
                    look_back=look_back,
                    lag=lag,
                    n_cycles=n_cycles,
                    include_confidence=include_confidence,
                    include_flood=include_flood,
                    fusion_samples=fusion_samples,
                    output_asset_path=output_asset_tile,
                    output_bucket_path=output_bucket_tile,
                    initial_threshold=initial_threshold,
                    thresh_no_data=thresh_no_data,
                    tile=False,
                    tile_buffer=tile_buffer,
                    output_scale=output_scale,
                )

    else:
        if not isinstance(target_date, ee.Date):
            target_date = ee.Date(target_date)

        end_time = target_date.advance(-(lag - 1), "day")
        start_time = end_time.advance(-look_back, "day")

        if fusion_samples is not None:
            fusion_model, scaling_dict = ml.random_forest_ee(
                30,
                fusion_samples,
                feature_names,
                label,
                scaling=None,
                mode="regression",
            )
        else:
            raise ValueError(
                "'fusion_samples' needs to be defined to run fusion process"
            )

        now = datetime.datetime.now()
        time_id = now.strftime("%Y%m%d%H%M%s")
        time_str = now.strftime("%Y-%m-%d %H:%M:%s")

        if harmonic_image is not None:
            harmonic_coefs = ee.Image(harmonic_image)
            harmonic_coefs = harmonic_coefs.multiply(
                ee.Image(ee.Number(harmonic_coefs.get("scale_factor")))
            )
        elif harmonic_collection is not None:
            harmonic_collection = ee.ImageCollection(harmonic_collection)
            first = ee.Image(harmonic_collection.first())
            harmonic_coefs = ee.Image(harmonic_collection.mosaic()).multiply(
                ee.Image(ee.Number(first.get("scale_factor")))
            )
        else:
            raise ValueError(
                "Either 'harmonic_image' or 'harmonic_collection' needs to be defined to run fusion process"
            )

        if include_confidence:
            harmonic_err = harmonic_coefs.select(".*(x|y|n)$")
            harmonic_coefs = harmonic_coefs.select("^(c|t|s).*")
        else:
            harmonic_coefs = harmonic_coefs.select("^(c|t|s).*")

        prod_region = region.buffer(tile_buffer, 100)

        ds = _fuse_dataset(
            region,
            start_time,
            end_time,
            fusion_model,
            scaling_dict,
            feature_names,
            target_band=label,
            use_viirs=True,
        )

        dummy_target = timeseries.get_dummy_img(target_date)

        weights = ee.ImageCollection.fromImages(
            ee.List.sequence(0, look_back - 1).map(get_residuals)
        ).sort("system:time_start")

        weights_lr = timeseries.fit_linear_trend(
            weights, dependent="residual", output_err=include_confidence
        )

        weights_coefs = weights_lr.select("^(c|t).*")

        lin_pred = (
            dummy_target.multiply(weights_coefs).reduce("sum").rename("residual_est")
        )

        har_pred = (
            timeseries.add_harmonic_coefs(dummy_target, n_cycles=n_cycles)
            .multiply(harmonic_coefs)
            .reduce("sum")
        )

        fused_pred = (har_pred.subtract(lin_pred)).rename("fused_product")

        ci_threshold = thresholding.edge_otsu(
            fused_pred,
            initial_threshold=initial_threshold,
            thresh_no_data=thresh_no_data,
            edge_buffer=300,
            region=prod_region,
            invert=True,
            scale=150,
            return_threshold=True,
        )

        permanent_water = (
            ee.ImageCollection("JRC/GSW1_2/YearlyHistory")
            .filterDate("1985-01-01", end_time)
            .limit(5, "system:time_start", False)
            .map(lambda x: x.select("waterClass").eq(3))
            .sum()
            .unmask(0)
            .gt(0)
        )

        water = fused_pred.gt(ci_threshold).Or(permanent_water).rename("water").uint8()

        if include_flood:
            flood = water.select("water").And(permanent_water.Not()).rename("flood")
            water = water.addBands(flood)

        if include_confidence:
            weights_err = weights_lr.select(".*(x|y|n)$")

            linCi = weights_err.expression(
                "mse * (1 + (1/n) + ((t-xmean)**2/xr))**(1/2)",
                {
                    "mse": weights_err.select("residual_y"),
                    "n": weights_err.select("n"),
                    "xmean": weights_err.select("mean_x"),
                    "xr": weights_err.select("residual_x"),
                    "t": dummy_target.select("time"),
                },
            )

            harCi = harmonic_err.expression(
                "mse * (1 + (1/n) + ((t-xmean)**2/xr))**(1/2)",
                {
                    "mse": harmonic_err.select("residual_y"),
                    "n": harmonic_err.select("n"),
                    "xmean": harmonic_err.select("mean_x"),
                    "xr": harmonic_err.select("residual_x"),
                    "t": dummy_target.select("time"),
                },
            )

            confidence = (
                ee.ImageCollection.fromImages(
                    ee.List.sequence(0, 99).map(calc_confidence)
                )
                .reduce(ee.Reducer.mean(), 16)
                .multiply(100)
                .uint8()
                .rename("confidence")
            )

            out_water = water.addBands(confidence)
        else:
            out_water = water

        if output_asset_path is not None:
            # create metadata dict
            metadata = ee.Dictionary(
                {
                    "hf_version": hf.__version__,
                    "system:time_start": target_date.millis(),
                    "system:time_end": target_date.advance(86399, "seconds").millis(),
                    "execution_time": time_str,
                    "lag": lag,
                    "look_back": look_back,
                }
            )
            geeutils.export_image(
                out_water.set(metadata.combine({"product": "water"})),
                region,
                output_asset_path + "_water",
                description=f"hydrafloods_water_ee_export_{time_id}",
                scale=output_scale,
                crs="EPSG:4326",
            )

        elif output_bucket_path is not None:
            export_region = region.bounds(maxError=100).getInfo()["coordinates"]
            bucket_path, ext = os.path.splitext(output_bucket_path)
            fcomponents = bucket_path.split("/")
            bucket = fcomponents[2]
            fpath = fcomponents[3:-1]

            # TODO: remove extension from string formulation
            f_water = "/".join(fpath + [fcomponents[-1] + "_water" + ext])
            f_fusion = "/".join(fpath + [fcomponents[-1] + "_fusion" + ext])

            water_task = ee.batch.Export.image.toCloudStorage(
                image=out_water,
                description=f"hydrafloods_water_gcp_export_{time_id}",
                bucket=bucket,
                fileNamePrefix=f_water,
                region=export_region,
                scale=output_scale,
                crs="EPSG:4326",
                maxPixels=1e13,
                fileFormat="GeoTIFF",
                formatOptions={"cloudOptimized": True},
            )
            water_task.start()

        else:
            raise ValueError(
                "Either 'output_asset_path' or 'output_bucket_path' needs to be defined to run fusion export process"
            )

    return


def merge_gcp_tiled_results(
    bucket_path,
    pattern,
    region,
    retries=-1,
    clean_up=False,
    cloud_project=None,
    file_dims=None,
    output_scale=30,
):
    """Helper function to merge tiled surface water estimates from `export_daily_surface_water`
    into on cloud optimized geotiff on Google Cloud Platform

    args:
        bucket_path (str): GCP cloud bucket path to read tiled cloud optimized geotiffs. Will
            also export merged file to this path.
        pattern (str): regex string to search for specific files. Useful for selecting files from
            a specific bucket subdirectory or date in filename
        region (ee.Geometry): geographic region to export merged data over. region must align with
            region from the tiled export
        retries (int, optional): number of retries to search for tiled data. Useful is runtime is unknown
            and the merge function is run at a set time each day. If less than or equal to zero, no retries
            will be used. default =  -1
        clean_up (bool, optional): boolean keyword to delete tile geotiffs after merge is complete. Only use
            if you are sure the merging is/will be successful. default = False
        cloud_project (str, optional): name of GCP cloud project name that the cloud storage bucket is part
            of. If nothing is provided, it will try to use default system GCP credentials. default = None
        file_dims (int | list[int], optional): the dimensions in pixels of each image file, if the image
            is too large to fit in a single file. May specify a single number to indicate a square shape,
            or a list of two dimensions to indicate (width,height). Note that the image will still be
            clipped to the overall image dimensions. Must be a multiple of shardSize (256). If none, then
            Earth Engine will automatically estimate dimensions. default = None
        output_scale (float, optional): output resolution of harmonic weight image. default = 30
    """
    land_area = (
        ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
        .filterBounds(region)
        .geometry(100)
        .buffer(2500, maxError=100)
    )
    grid = geeutils.tile_region(region, intersect_geom=land_area, grid_size=1.0)

    expected_n = grid.size().getInfo()

    fcomponents = bucket_path.split("/")
    bucket = fcomponents[2]
    fpath = pattern.replace("*", "")

    files = utils.list_gcs_objs(bucket, pattern=pattern, project=cloud_project)

    if len(files) == expected_n:
        images = [ee.Image.loadGeoTIFF(file) for file in files]

        merged = ee.ImageCollection.fromImages(images).mosaic()

        export_region = region.bounds(maxError=100).getInfo()["coordinates"]
        # bucket_path, ext = os.path.splitext(output_bucket_path)

        now = datetime.datetime.now()
        time_id = now.strftime("%Y%m%d%H%M%s")

        task = ee.batch.Export.image.toCloudStorage(
            image=merged,
            description=f"hydrafloods_merge_gcp_export_{time_id}",
            bucket=bucket,
            fileNamePrefix=fpath,
            region=export_region,
            scale=output_scale,
            crs="EPSG:4326",
            fileDimensions=file_dims,
            maxPixels=1e13,
            fileFormat="GeoTIFF",
            formatOptions={"cloudOptimized": True},
        )
        task.start()

        if clean_up:
            gcsfs.GCSFileSystem.rm(files)

    elif retries > 0:
        time.sleep(60 * 10)

        merge_gcp_tiled_results(bucket_path, pattern, region, retries=(retries - 1))

    else:
        raise RuntimeError(
            f"could not find all expected tiles to merge...found {len(files)} tiles but expected {expected_n}"
        )

    return


if __name__ == "__main__":
    raise NotImplementedError(
        "Worklow is currently not implemented for CLI use, please check back later"
    )
