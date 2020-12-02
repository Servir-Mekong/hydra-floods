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
)


def export_fusion_samples(
    region,
    start_time,
    end_time,
    stratify_samples=False,
    sample_scale=30,
    n_samples=25,
    output_asset_path=None,
    seed=0,
):
    """
    """

    optical_water_indices = ["mndwi", "nwi", "aewish", "aewinsh"]

    ds_kwargs = dict(
        region=region, start_time=start_time, end_time=end_time, rescale=True
    )
    dsa_kwargs = {**ds_kwargs, **{"apply_band_adjustment": True}}

    lc8 = datasets.Landsat8(**ds_kwargs)
    le7 = datasets.Landsat7(**dsa_kwargs)
    s2 = datasets.Sentinel2(**dsa_kwargs)

    _ = ds_kwargs.pop("rescale")

    s1 = datasets.Sentinel1(**ds_kwargs)
    # s1.collection = timeseries.temporal_iqr_filter(s1.collection)
    s1 = s1.apply_func(
        geeutils.add_indices, indices=["vv_vh_ratio", "ndpi", "vv_vh_abs_sum","nvvi","nvhi"]
    )

    optical = lc8.merge(s2).merge(le7)
    optical = optical.apply_func(geeutils.add_indices, indices=optical_water_indices)

    optical.collection = optical.collection.select(optical_water_indices)

    ds = optical.join(s1)

    n = ds.n_images
    logging.info(f"Found {n} images to sample from")

    img_list = ds.collection.toList(n)

    output_features = ee.FeatureCollection([])

    if stratify_samples:
        class_band = "landcover"
        igbp_classes = ee.List(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        )
        ipcc_classes = ee.List([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 6, 3, 5, 3, 4, 4, 6])

        stratification_img = (
            ee.ImageCollection("MODIS/006/MCD12Q1")
            .limit(1, "system:time_start", True)
            .first()
            .remap(igbp_classes, ipcc_classes)
            .rename(class_band)
        )
        classes = ipcc_classes.distinct()
    else:
        stratification_img = None

    for i in range(n):
        try:
            sample_img = ee.Image(img_list.get(i))

            sample_region = sample_img.geometry(sample_scale).bounds(sample_scale)

            if stratification_img is not None:

                samples = (
                    sample_img.addBands(stratification_img.select(class_band))
                    .stratifiedSample(
                        region=sample_region,
                        numPoints=n_samples,
                        classBand=class_band,
                        scale=sample_scale,
                        seed=seed + i,
                        classValues=classes,
                        classPoints=ee.List.repeat(
                            n_samples, classes.size().subtract(1)
                        ).add(n_samples * 4),
                        tileScale=16,
                        geometries=True,
                    )
                )

            else:
                samples = sample_img.sample(
                    region=sample_region,
                    scale=sample_scale,
                    numPixels=n_samples,
                    seed=seed + i,
                    tileScale=16,
                    geometries=True,
                )

            output_features = (
                output_features.merge(samples.randomColumn())
                if samples.size().getInfo() > 0
                else output_features
            )

        except EEException as e:
            warnings.warn(f"Sampling process ran into an error: {str(e)}")
            break

    task = ee.batch.Export.table.toAsset(
        collection=output_features, assetId=output_asset_path
    )
    task.start()
    logging.info(f"Started export task for {output_asset_path}")

    return


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
    @decorators.carry_metadata
    def _apply_fusion(img):
        img_norm = ml.standard_image_scaling(img, scaling_dict, feature_names)
        return img_norm.classify(fusion_model, target_band)

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
    s1 = s1.apply_func(
        geeutils.add_indices, indices=["vv_vh_ratio", "ndpi", "nvvi", "nvhi"]
    )

    s1 = s1.apply_func(_apply_fusion)

    fused_ds = optical.merge(s1)

    fused_ds.collection = (
        fused_ds.collection.select(target_band)
        .cast({target_band: "float"}, [target_band])
        .sort("system:time_start")
    )

    return fused_ds


def export_surface_water_harmonics(
    region,
    start_time,
    end_time,
    n_cycles=2,
    feature_names=None,
    label=None,
    fusion_samples=None,
    output_asset_path=None,
    output_bucket=None,
    tile=False,
    tile_size=1.0,
    output_scale=30,
):

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
                25,
                fusion_samples.limit(10000),
                feature_names,
                label,
                scaling="standard",
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
        elif output_bucket is not None:
            raise NotImplementedError()
        else:
            raise ValueError(
                "Either 'output_asset_path' or 'output_bucket' needs to be defined to run fusion export process"
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
    export_fusion=False,
    fusion_samples=None,
    output_asset_path=None,
    output_bucket_path=None,
    initial_threshold=0.1,
    tile=False,
    tile_size=1.0,
    tile_buffer=100000,
    output_scale=30,
):
    def get_weights(i):
        i = ee.Number(i)
        t_diff = (
            ee.Number(i).multiply(-1).subtract(lag)
        )  # calc how many days to adjust ini date
        new_date = target_date.advance(t_diff, "day")  # calculate new date

        corr_img = (
            ds.collection.select(label)
            .filterDate(new_date, new_date.advance(1, "day"))
            .qualityMosaic(label)
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
                    export_fusion=export_fusion,
                    fusion_samples=fusion_samples,
                    output_asset_path=output_asset_tile,
                    output_bucket_path=output_bucket_tile,
                    initial_threshold=initial_threshold,
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
                25,
                fusion_samples.limit(10000),
                feature_names,
                label,
                scaling="standard",
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
            use_viirs=True
        )

        dummy_target = timeseries.get_dummy_img(target_date)

        weights = ee.ImageCollection.fromImages(
            ee.List.sequence(0, look_back - 1).map(get_weights)
        ).sort("system:time_start")

        weights_lr = timeseries.fit_linear_trend(
            weights, dependent="residual", output_err=include_confidence
        )

        weights_coefs = weights_lr.select("^(c|t).*")

        lin_pred = (
            dummy_target.multiply(weights_coefs).reduce("sum").rename("residual_est")
        )

        har_pred = (
            timeseries.add_harmonic_coefs(dummy_target,n_cycles=n_cycles)
            .multiply(harmonic_coefs)
            .reduce("sum")
        )

        fused_pred = hf.filtering.p_median((har_pred.subtract(lin_pred))).rename(
            "fused_product"
        )

        ci_threshold = thresholding.edge_otsu(
            fused_pred,
            initial_threshold=initial_threshold,
            edge_buffer=300,
            region=prod_region,
            invert=True,
            scale=200,
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

            out_water = ee.Image.cat([confidence, water,])
        else:
            out_water = water

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
            if export_fusion:
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

            if export_fusion:
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


def merge_gcp_tiled_results(
    bucket_path,
    pattern,
    region,
    retries=-1,
    clean_up=False,
    cloud_project=None,
    file_dims=None,
):
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
            scale=10,
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
