# daily surface water - fusion process

import os
import ee
from ee.ee_exception import EEException
import gcsfs
import logging
import datetime
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing as mp
from functools import partial
from pprint import pformat
from sklearn import metrics, model_selection, preprocessing
import hydrafloods as hf
from hydrafloods import (
    datasets,
    timeseries,
    ml,
    utils,
    geeutils,
    thresholding,
    decorators,
)

# temporary to see what is going on
logging.basicConfig(level=logging.INFO)


def export_fusion_samples(
    region,
    start_time,
    end_time,
    stratify_samples=False,
    sample_scale=30,
    n_samples=100,
    img_limit=1000,
    export_to="asset",
    output_asset_path=None,
    export_kwargs=None,
    skip_empty=True,
):
    """
    """

    export_opts = dict(
        cloud=ee.batch.Export.table.toCloudStorage, asset=ee.batch.Export.table.toAsset,
    )
    export_func = export_opts[export_to]

    ds_kwargs = dict(region=region, start_time=start_time, end_time=end_time)
    dsa_kwargs = {**ds_kwargs, **{"apply_band_adjustment": True}}

    lc8 = datasets.Landsat8(**ds_kwargs)
    le7 = datasets.Landsat7(**dsa_kwargs)
    s2 = datasets.Sentinel2(**dsa_kwargs)

    s1 = datasets.Sentinel1(**ds_kwargs)
    s1 = s1.add_fusion_features()

    optical = lc8.merge(s2).merge(le7)

    ds = optical.join(s1)

    n = img_limit if img_limit is not None else ds.n_images
    img_list = ds.collection.toList(n)

    output_features = ee.FeatureCollection([])

    if stratify_samples:
        stratification_img = (
            ee.Image("JRC/GSW1_2/GlobalSurfaceWater")
            .select("occurrence")
            .unmask(0)
            .gt(60)
            .rename("water")
        )
    else:
        stratification_img = None

    for i in range(n):
        try:
            sample_img = ee.Image(img_list.get(i))

            sample_region = sample_img.geometry().bounds()

            if stratification_img is not None:
                class_band = "water"
                classes = ee.List([0, 1])

                samples = sample_img.addBands(
                    stratification_img.select(class_band)
                ).stratifiedSample(
                    region=sample_region,
                    numPoints=n_samples,
                    classBand=class_band,
                    scale=sample_scale,
                    seed=i,
                    classValues=classes,
                    classPoints=ee.List.repeat(n_samples, classes.size()),
                    tileScale=16,
                    geometries=True,
                )

            else:
                samples = sample_img.sample(
                    region=sample_region,
                    scale=sample_scale,
                    numPixels=n_samples,
                    seed=i,
                    tileScale=16,
                    geometries=True,
                )

            if skip_empty:
                output_features = (
                    output_features.merge(samples)
                    if samples.size().getInfo() > 0
                    else output_features
                )
            else:
                output_features = output_features.merge(samples)

        except EEException as e:
            break

    export_info = dict(collection=output_features, assetId=output_asset_path)
    print(export_info)
    if export_kwargs is not None:
        export_info = {**export_info, **export_kwargs}
        # if "fileNamePrefix" in export_info.keys():
        #     prefix = export_kwargs["fileNamePrefix"]
        #     true_prefix = (
        #         prefix + desc
        #         if prefix.endswith("/")
        #         else prefix + f"_{dstr}_{i}"
        #     )
        #     export_info["fileNamePrefix"] = true_prefix

    task = export_func(collection=output_features, assetId=output_asset_path)
    task.start()
    logging.info(f"Started task")

    return


def build_fusion_model(
    sample_path,
    features,
    label,
    export_scaler_to_ee=True,
    ee_asset_path=None,
    output_bucket=None,
    framework="sklearn",
    output_training_report=True,
    filter_outliers=False,
    seed=0,
):

    # framework options = [sklearn, xgboost, lightgbm]

    if sample_path.startswith("gs://"):
        tables = utils.list_gcs_objs(sample_path.replace("gs://", ""), pattern="*.csv")
    else:
        raise NotImplementedError(
            "Currently only fetching data from Google Cloud Storage is supported"
        )

    df_list = (pd.read_csv(table) for table in tables)
    df = pd.concat(df_list, axis=0, ignore_index=True)

    X = df[features]
    y = df[label]

    feature_names = list(X.columns)
    n_features = len(feature_names)

    if filter_outliers:
        logging.info(f"applying filters on feature coluns")
        α = 0.05
        z_threshold = 3
        masks = np.ones(df.shape)
        for i, feature in enumerate(list(X.columns)):
            values = X[feature]
            # ratio of number of unique values to the total number of unique values
            # probability is less than α thereshold then do not filter outliers
            is_categorical = 1.0 * np.unique(values).size / values.size < α
            if is_categorical:
                logging.info(f"{feature} is categorial, passing filter")
                continue

            k2, p = stats.normaltest(values)
            if p < α:
                logging.info(f"{feature} is gaussian, using z-score filter")
                z = np.abs(stats.zscore(values))
                masks[:, i] = z < z_threshold

            else:
                logging.info(f"{feature} is not gaussian, using iqr filter")
                q1, q3 = stats.iqr(values)
                iqr = q3 - q1
                masks[:, i] = (values > (q1 - 1.5 * iqr)) & (values < (q3 + 1.5 * iqr))

        X = X[masks.any(axis=1)]
        y = y[masks.any(axis=1)]

    if output_training_report:
        X_train, X_test, y_train, y_test = model_selection.train_test_split(
            X, y, train_size=0.80, random_state=seed
        )
    else:
        X_train, y_train = X, y

    scaler = preprocessing.MinMaxScaler().fit(X_train)
    X_train = scaler.transform(X_train)

    fmin = {f"{feature_names[i]}_min": scaler.data_min_[i] for i in range(n_features)}
    fmax = {f"{feature_names[i]}_max": scaler.data_max_[i] for i in range(n_features)}
    scaler_fc = ee.FeatureCollection(
        ee.Feature(ee.Geometry.Point([0, 0]), {**fmin, **fmax})
    )

    now = datetime.datetime.now()
    time_id = now.strftime("%Y%m%d%H%M%s")

    if export_scaler_to_ee:

        desc = f"fusion_model_feature_scaling_{time_id}"
        task = ee.batch.Export.table.toAsset(
            collection=scaler_fc, description=desc, assetId=ee_asset_path + desc
        )
        task.start()

    logging.info(f"training model...")
    if framework is "sklearn":
        from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor

        # model = RandomForestRegressor(n_estimators=100, n_jobs=-1, random_state=seed,)
        model = ExtraTreesRegressor(
            n_estimators=50, n_jobs=-1, max_depth=30, random_state=0
        )

    else:
        raise NotImplementedError(
            "only fusion modeling using the scikit-learn framework is avaiable now"
        )

    t1 = datetime.datetime.now()
    model.fit(X_train, y_train)
    training_time = datetime.datetime.now() - t1

    if output_training_report:
        X_test = scaler.transform(X_test)
        y_pred = model.predict(X_test)

        mae = metrics.mean_absolute_error(y_test, y_pred)
        me = metrics.max_error(y_test, y_pred)
        r2 = metrics.r2_score(y_test, y_pred)
        bias = np.mean((y_test - y_pred))
        rmse = np.mean(np.sqrt((y_test - y_pred) ** 2))

        fs = gcsfs.GCSFileSystem()

        with fs.open(f"{output_bucket}/training_report_{time_id}.txt", "w") as report:
            content = [
                f"Model framework: {framework}",
                f"Execution time: {now}",
                f"Model parameters: {{\n {pformat(model.get_params(),indent=8)[1:]}",
                f"Training information:",
                f"\trandom seed: {seed}",
                f"\ttraining time: {training_time}",
                f"\tn training examples: {X_train.shape[0]}",
                f"\tn testing examples: {X_test.shape[0]}",
                f"\tbias: {bias}",
                f"\trmse: {rmse}",
                f"\tmae: {mae}",
                f"\tr2: {r2}",
                f"\tmax error: {me}",
                f"Feature scaling:",
                f"\tfeature minimum: {{\n {pformat(fmin,indent=16)[1:]}",
                f"\tfeature maximum: {{\n {pformat(fmax,indent=16)[1:]}",
            ]
            if export_scaler_to_ee:
                content.append(f"\tOutput scaling feature collection: {desc}")

            report.write("\n".join(content))

    estimators = model.estimators_

    # fs = gcsfs.GCSFileSystem()
    # for i, estimator in enumerate(estimators):
    #     string = ml.sklearn_tree_to_string(estimator, features)
    #     with fs.open(f"{output_bucket}/estimator_{time_id}_{i:04d}.txt", "w") as f:
    #         f.write(string)

    tree_properties = {}
    # args = [(est,features) for est in estimators]
    with mp.Pool(7) as pool:
        proc = pool.map_async(
            partial(ml.sklearn_tree_to_string, feature_names=features), estimators
        )
        trees = list(proc.get())

    df = pd.DataFrame(
        {
            "lat": np.zeros(len(trees)),
            "lon": np.zeros(len(trees)),
            "tree": [i.replace("\n", "#") for i in trees],
        }
    )

    bucket_obj = f"gs://{output_bucket}/rf_model_{time_id}.csv"
    df.to_csv(bucket_obj, index=False)

    os.system(
        f"earthengine upload table {bucket_obj} --asset_id users/kelmarkert/rf_tree_test --x_column lon --y_column lat"
    )

    # tree_properties = {f"tree_{i:30d}": trees[i] for i in range(len(trees))}

    # features = [ee.Feature(ee.Geometry.Point([0,0])).set('tree',tree_properties[k]) for k in tree_properties.keys()]

    # ee_tree_obj= ee.FeatureCollection(features)

    # task = ee.batch.Export.table.toAsset(
    #     collection=ee_tree_obj,
    #     description="rf_tree_export_test",
    #     assetId="users/kelmarkert/rf_tree_test"
    # )
    # task.start()

    return


def _fuse_dataset(
    region,
    start_time,
    end_time,
    fusion_model,
    scaling_dict=None,
    target_band="mndwi",
    use_viirs=False,
    use_modis=False,
):
    @decorators.carry_metadata
    def _apply_scaling(img):
        return img.subtract(min_img).divide(max_img.subtract(min_img)).float()

    @decorators.carry_metadata
    def _apply_fusion(img):
        return img.classify(fusion_model, target_band)

    ds_kwargs = dict(region=region, start_time=start_time, end_time=end_time)
    dsa_kwargs = {**ds_kwargs, **{"apply_band_adjustment": True}}

    lc8 = datasets.Landsat8(**ds_kwargs)
    le7 = datasets.Landsat7(**dsa_kwargs)
    s2 = datasets.Sentinel2(**dsa_kwargs)

    optical = lc8.merge(le7).merge(s2)

    if use_viirs:
        viirs = datasets.Viirs(**ds_kwargs)
        optical = optical.merge(viirs)

    if use_modis:
        modis = datasets.Modis(**ds_kwargs)
        optical = optical.merge(modis)

    s1 = datasets.Sentinel1(**ds_kwargs)
    s1 = s1.add_fusion_features()

    if scaling_dict is not None:
        scaling_img = scaling_dict.toImage()
        min_img = scaling_img.select(".*_min")
        max_img = scaling_img.select(".*_max")
        s1 = s1.apply_func(_apply_scaling)

    s1.collection = s1.collection.map(_apply_fusion)

    fused_ds = optical.merge(s1)

    fused_ds.collection = (
        fused_ds.collection.select(target_band)
        .cast({target_band: "float"}, [target_band])
        .sort("system:time_start")
    )

    return fused_ds, target_band


def export_harmonics(
    region,
    start_time,
    end_time,
    feature_names=None,
    label=None,
    fusion_samples=None,
    fusion_model_path=None,
    output_asset_path=None,
    output_bucket=None,
    tile=False,
):
    def constuctGrid(i):
        def contructXGrid(j):
            j = ee.Number(j)
            box = ee.Feature(
                ee.Geometry.Rectangle(j, i, j.add(grid_size), i.add(grid_size))
            )
            out = ee.Algorithms.If(
                land_area.intersects(box.geometry(), maxError=100), box, None
            )
            return ee.Feature(out)

        i = ee.Number(i)
        out = ee.List.sequence(west, east.subtract(grid_size), grid_size).map(
            contructXGrid
        )
        return out

    land_area = (
        ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
        .filterBounds(region)
        .geometry(100)
        .buffer(2500, maxError=100)
    )

    if fusion_samples is not None:
        fusion_model, scaling_dict = ml.random_forest_ee(
            25, fusion_samples, feature_names, label, mode="regression"
        )
    elif model_estimator_path is not None:
        raise NotImplementedError()
    else:
        raise ValueError(
            "Either 'fusion_samples' or 'fusion_model_path' needs to be defined to run fusion process"
        )

    ds, label = _fuse_dataset(
        region, start_time, end_time, fusion_model, scaling_dict, target_band="mndwi"
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
        ds, dependent="mndwi", output_err=True
    )
    harmonic_coefs = harmonic_coefs.divide(scale_factor).int32().set(metadata)

    if tile:
        bounds = region.bounds(maxError=100)
        coords = ee.List(bounds.coordinates().get(0))
        grid_size = ee.Number(1.0)

        west = ee.Number(ee.List(coords.get(0)).get(0))
        south = ee.Number(ee.List(coords.get(0)).get(1))
        east = ee.Number(ee.List(coords.get(2)).get(0))
        north = ee.Number(ee.List(coords.get(2)).get(1))

        west = west.subtract(west.mod(grid_size))
        south = south.subtract(south.mod(grid_size))
        east = east.add(grid_size.subtract(east.mod(grid_size)))
        north = north.add(grid_size.subtract(north.mod(grid_size)))

        grid = ee.FeatureCollection(
            ee.List.sequence(south, north.subtract(grid_size), grid_size)
            .map(constuctGrid)
            .flatten()
        )

        n = grid.size().getInfo()
        grid_list = grid.toList(n)

        for i in range(n):
            out_region = ee.Feature(grid_list.get(i)).geometry()
            output_asset_img = f"{output_asset_path}/img_{i}"
            geeutils.export_image(
                harmonic_coefs,
                out_region,
                output_asset_img,
                description=f"hydrafloods_harmonic_coefficient_export_{i}_{time_id}",
                scale=10,
                crs="EPSG:4326",
            )

    else:
        if output_asset_path is not None:
            geeutils.export_image(
                harmonic_coefs,
                region,
                output_asset_path,
                description=f"hydrafloods_harmonic_coefficient_export_{time_id}",
                scale=10,
                crs="EPSG:4326",
            )
        elif output_bucket is not None:
            raise NotImplementedError()
        else:
            raise ValueError(
                "Either 'output_asset_path' or 'output_bucket' needs to be defined to run fusion export process"
            )

    return


def export_harmonic_daily_surface_water(
    region,
    target_date,
    harmonic_coefs,
    feature_names=None,
    label=None,
    look_back=20,
    decay_factor=0.25,
    lag=1,
    fusion_samples=None,
    fusion_model_path=None,
    output_asset_path=None,
    output_bucket=None,
):
    def _exponential_decay_correction(t, look_back, decay_factor, lag):
        # lambda =  decay rate for time weights, higher means weighting resent info more
        # lookBack = how many days to use for weights
        # lag = how many days back to start the lookBack

        def _get_weight(i):
            i = ee.Number(i)
            tDiff = (
                ee.Number(i).multiply(-1).subtract(lag)
            )  # calc how many days to adjust ini date
            new_date = t.advance(tDiff, "day")  # calculate new date

            # get the imagery for the day
            correction_imgs = ds.collection.filterDate(
                new_date, new_date.advance(1, "day")
            ).median()

            # get dummy image to predict harmonics from
            dummy = timeseries.get_dummy_img(new_date)

            # predict harmonics for lag date
            harmonic_est = (
                timeseries.add_harmonic_coefs(dummy)
                .multiply(harmonic_coefs)
                .reduce("sum")
                .rename("estimate")
            )

            # unmask any data to prevent no data gaps
            # unmasked area will have correction of 0 (no adjustment)
            correction = ee.Image(
                ee.Algorithms.If(
                    correction.bandNames().length().lt(1), harmonic_est, correction,
                )
            )

            # calculate the correction factor and apply weight
            residual = harmonic_est.subtract(correction).multiply(
                ee.Image(decay_weights.get([i.int()]))
            )

            return residual.set("system:time_start", new_date.millis())

        t = ee.Date(t)  # initial time

        # calculate exponential decay weights
        decay_weights = ee.Array(
            ee.List.sequence(0, (look_back - 1))
            .map(lambda x: ee.Number(decay_factor).multiply(ee.Number(x).add(1)).exp())
            .reverse()
        )
        # reverse so recent data = higher

        # force weights to sum to 1
        decay_sum = decay_weights.reduce(ee.Reducer.sum(), [0]).get([0])
        decay_weights = decay_weights.divide(decay_sum)

        # calculate composite images and apply weights
        weights = ee.ImageCollection(
            ee.List.sequence(0, (look_back - 1)).map(_get_weight)
        )

        return weights  # image collection so will need to sum later

    harmonic_coefs = ee.Image(harmonic_coefs)
    harmonic_coefs = harmonic_coefs.multiply(
        ee.Image(ee.Number(harmonic_coefs.get("scale_factor")))
    )

    end_time = ee.Date(target_date)
    start_time = end_time.advance(-look_back, "day")

    if fusion_samples is not None:
        fusion_model, scaling_dict = ml.random_forest_ee(
            25, fusion_samples, feature_names, label, mode="regression"
        )
    elif model_estimator_path is not None:
        raise NotImplementedError()
    else:
        raise ValueError(
            "Either 'fusion_samples' or 'fusion_model_path' needs to be defined to run fusion process"
        )

    ds, label = _fuse_dataset(
        region,
        start_time,
        end_time,
        fusion_model,
        scaling_dict,
        target_band=label,
        use_viirs=False,
    )

    now = datetime.datetime.now()
    time_id = now.strftime("%Y%m%d%H%M%s")
    time_str = now.strftime("%Y-%m-%d %H:%M:%s")

    dummy = timeseries.get_dummy_img(end_time)

    harmonic_est = (
        timeseries.add_harmonic_coefs(dummy)
        .multiply(harmonic_coefs)
        .reduce("sum")
        .rename("estimate")
    )
    weights = _exponential_decay_correction(end_time, look_back, decay_factor, lag)
    correction = weights.sum()

    fused_index = harmonic_est.subtract(correction)

    # create metadata dict
    metadata = ee.Dictionary(
        {
            "hf_version": hf.__version__,
            "system:time_start": end_time.millis(),
            "execution_time": time_str,
        }
    )

    # water_est = thresholding.edge_otsu(
    #     fused_index, initialThreshold=0, smooth_edges=300, invert=True
    # ).set(metadata)

    if output_asset_path is not None:
        geeutils.export_image(
            fused_index.set(metadata),
            region,
            output_asset_path,
            description=f"hydrafloods_fused_water_export",
            scale=10,
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
    harmonic_coefs=None,
    harmonic_collection=None,
    feature_names=None,
    label=None,
    look_back=30,
    lag=4,
    output_confidence=False,
    fusion_samples=None,
    fusion_model_path=None,
    output_asset_path=None,
    output_bucket_path=None,
    output_dims=None,
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
        random = (
            ee.Image.random(i).multiply(3.92).subtract(1.96)
        )  # uniform sampling of std dev at 95% confidence interval

        lin_sim = lin_pred.add(random.multiply(linCi))
        har_sim = har_pred.add(random.multiply(harCi))

        sim_pred = har_sim.subtract(lin_sim)
        # random_water = thresholding.bmax_otsu(random_combination,invert=True)
        # naive estimate of water (>0)
        return sim_pred.gt(0).uint8()

    end_time = ee.Date(target_date).advance(-(lag - 1), "day")
    start_time = end_time.advance(-look_back, "day")

    if fusion_samples is not None:
        fusion_model, scaling_dict = ml.random_forest_ee(
            25, fusion_samples, feature_names, label, mode="regression"
        )
    elif fusion_model_path is not None:
        raise NotImplementedError()
    else:
        raise ValueError(
            "Either 'fusion_samples' or 'fusion_model_path' needs to be defined to run fusion process"
        )

    if not isinstance(target_date, ee.Date):
        target_date = ee.Date(target_date)

    now = datetime.datetime.now()
    time_id = now.strftime("%Y%m%d%H%M%s")
    time_str = now.strftime("%Y-%m-%d %H:%M:%s")

    if harmonic_coefs is not None:
        harmonic_coefs = ee.Image(harmonic_coefs)
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
            "Either 'harmonic_coefs' or 'harmonic_collection' needs to be defined to run fusion process"
        )

    if output_confidence:
        harmonic_err = harmonic_coefs.select(".*(x|y|n)$")
        harmonic_coefs = harmonic_coefs.select("^(c|t|s).*")
    else:
        harmonic_coefs = harmonic_coefs.select("^(c|t|s).*")

    ds, label = _fuse_dataset(
        region,
        start_time,
        end_time,
        fusion_model,
        scaling_dict,
        target_band=label,
        use_viirs=True,
        use_modis=False,
    )

    dummy_target = timeseries.get_dummy_img(target_date)

    weights = ee.ImageCollection.fromImages(
        ee.List.sequence(0, look_back - 1).map(get_weights)
    ).sort("system:time_start")

    weights_lr = timeseries.fit_linear_trend(
        weights, dependent="residual", output_err=output_confidence
    )

    weights_coefs = weights_lr.select("^(c|t).*")

    lin_pred = dummy_target.multiply(weights_coefs).reduce("sum").rename("residual_est")

    har_pred = (
        timeseries.add_harmonic_coefs(dummy_target)
        .multiply(harmonic_coefs)
        .reduce("sum")
    )

    fused_pred = (har_pred.subtract(lin_pred)).convolve(ee.Kernel.gaussian(2.5))

    # water = thresholding.bmax_otsu(
    #     fused_pred,
    #     initial_threshold=0,
    #     grid_size=0.2,
    #     region=region,
    #     invert=True,
    #     reduction_scale=150,
    # )
    water = thresholding.edge_otsu(
        fused_pred,
        initial_threshold=0,
        edge_buffer=300,
        region=region,
        invert=True,
        reduction_scale=100,
    )

    if output_confidence:
        weights_err = weights_lr.select(".*(x|y|n)$")

        linCi = weights_err.expression(
            "mse * ((1/n) + ((t-xmean)**2/xr))**(1/2)",
            {
                "mse": weights_err.select("residual_y"),
                "n": weights_err.select("n"),
                "xmean": weights_err.select("mean_x"),
                "xr": weights_err.select("residual_x"),
                "t": dummy_target.select("time"),
            },
        )

        harCi = harmonic_err.expression(
            "mse * ((1/n) + ((t-xmean)**2/xr))**(1/2)",
            {
                "mse": harmonic_err.select("residual_y"),
                "n": harmonic_err.select("n"),
                "xmean": harmonic_err.select("mean_x"),
                "xr": harmonic_err.select("residual_x"),
                "t": dummy_target.select("time"),
            },
        )

        confidence = (
            ee.ImageCollection.fromImages(ee.List.sequence(0, 99).map(calc_confidence))
            .reduce(ee.Reducer.mean(), 16)
            .multiply(100)
            .uint8()
            .rename("confidence")
        )

        out_water = ee.Image.cat([confidence, water,])
    else:
        out_water = water
    # water.uint8().rename("water"),

    fused_pred = fused_pred.multiply(10000).int16().rename("fused_product")

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
            scale=10,
            crs="EPSG:4326",
        )
        geeutils.export_image(
            fused_pred.set(metadata.combine({"product": "fusion"})),
            region,
            output_asset_path + "_fusion",
            description=f"hydrafloods_fusion_ee_export_{time_id}",
            scale=10,
            crs="EPSG:4326",
        )

    elif output_bucket_path is not None:
        export_region = region.bounds(maxError=100).getInfo()["coordinates"]
        bucket_path, ext = os.path.splitext(output_bucket_path)
        fcomponents = bucket_path.split("/")
        bucket = fcomponents[2]
        fpath = fcomponents[3:-1]

        f_water = "/".join(fpath + [fcomponents[-1] + "_water" + ext])
        f_fusion = "/".join(fpath + [fcomponents[-1] + "_fusion" + ext])

        water_task = ee.batch.Export.image.toCloudStorage(
            image=out_water,
            description=f"hydrafloods_water_gcp_export_{time_id}",
            bucket=bucket,
            fileNamePrefix=f_water,
            region=export_region,
            scale=10,
            crs="EPSG:4326",
            maxPixels=1e13,
            fileFormat="GeoTIFF",
            fileDimensions=output_dims,
            formatOptions={"cloudOptimized": True},
        )
        water_task.start()

        fusion_task = ee.batch.Export.image.toCloudStorage(
            image=fused_pred,
            description=f"hydrafloods_fusion_gcp_export_{time_id}",
            bucket=bucket,
            fileNamePrefix=f_fusion,
            region=export_region,
            scale=10,
            crs="EPSG:4326",
            maxPixels=1e13,
            fileFormat="GeoTIFF",
            fileDimensions=output_dims,
            formatOptions={"cloudOptimized": True},
        )
        fusion_task.start()

    else:
        raise ValueError(
            "Either 'output_asset_path' or 'output_bucket_path' needs to be defined to run fusion export process"
        )

    return


if __name__ == "__main__":
    raise NotImplementedError(
        "Application is currently not implemented, please check back later"
    )
