import ee
from collections import OrderedDict
from hydrafloods import decorators


@decorators.keep_attrs
def apply_fcnn(
    image,
    project_name,
    model_name,
    model_kwargs=None,
    output_probas=False,
    output_names=None,
):
    """
    args:
        image (ee.Image): input image for FCNN model, must have all of the features as bands
        project_name (str): cloud project name to reference the model
        model_name (str): ai platform model name
        model_kwargs (dict, optional): dictionary of keyword arguments to pass to ee.Model. default = None
        output_probas (bool, optional): flag to set the output image as class probabilities. If False
            then the ouput will be a one band output of the classification. default = False
        output_bands (Iterable, optional): list of band names to set for the output image
    """

    if model_kwargs is None:
        model_kwargs = OrderedDict(projectName=project_name, modelName=model_name)
    else:
        positional = OrderedDict(projectName=project_name, modelName=model_name)

        model_kwargs = OrderedDict({**positional, **model_kwargs})

    # Load the trained model and use it for prediction.
    model = ee.Model.fromAiPlatformPredictor(**model_kwargs)

    # run the predictions
    predictions = model.predictImage(image.toFloat().toArray())

    if output_probas:
        if output_names is None:
            raise ValueError(
                "please provide `output_names` when `ouput_probas` is set to True"
            )

        output = predictions.arrayFlatten([output_names])

    else:
        output_names = "classification" if output_names is None else output_names
        # find highest probability class
        output = predictions.arrayArgmax().arrayFlatten([[output_names]])

    return output


def minmax_scaling_dict(fc, feature_names):
    """Function to calculate the minimum and maximum values of feautures in a collection
    Expects that fc has all feature names

    args:
        fc (ee.FeatureCollection): feature collection with the features used to calculate min/max value
        feature_names (list[str]):  names of feature columns to calculat min/max values from

    returns
        ee.Dictionary: dictionary of minimum and maximum values for each feature name
    """

    # force ee types
    fc = ee.FeatureCollection(fc)
    ee_feature_names = ee.List(feature_names)

    # apply reducer on each feature column
    feature_min_max = fc.reduceColumns(
        ee.Reducer.minMax().repeat(ee_feature_names.length()), ee_feature_names
    )

    # min/max feature names
    names = ee_feature_names.map(
        lambda x: ee.List([ee.String(x).cat("_min"), ee.String(x).cat("_max")])
    ).flatten()

    # get the min/max values for each feature
    # used to scale values from 0-1
    min_max_dict = ee.Dictionary.fromLists(
        names,
        ee.List(feature_min_max.get("min")).zip(feature_min_max.get("max")).flatten(),
    )

    return min_max_dict


def standard_scaling_dict(fc, feature_names):
    """Function to calculate the mean and standard deviation values of feautures in a collection
    Expects that fc has all feature names

    args:
        fc (ee.FeatureCollection): feature collection with the features used to calculate mean/std dev value
        feature_names (list[str]): names of feature columns to calculat mean/std dev values from

    returns
        ee.Dictionary: dictionary of mean and standard deviation values for each feature name
    """
    # force ee types
    fc = ee.FeatureCollection(fc)
    ee_feature_names = ee.List(feature_names)

    # get a combined reducer for caluclating mean and standard dev
    mean_stddev = ee.Reducer.mean().combine(ee.Reducer.stdDev(), None, True)

    # apply reducer on each feature column
    feature_mean_stddev = fc.reduceColumns(
        mean_stddev.repeat(ee_feature_names.length()), ee_feature_names
    )

    # mean / std dev feature names
    names = ee_feature_names.map(
        lambda x: ee.List([ee.String(x).cat("_mean"), ee.String(x).cat("_stdDev")])
    ).flatten()

    # get the mean / std dev values for each feature
    # used to scale values from ~ -3 to 3
    mean_stddev_dict = ee.Dictionary.fromLists(
        names,
        ee.List(feature_mean_stddev.get("mean"))
        .zip(feature_mean_stddev.get("stdDev"))
        .flatten(),
    )

    return mean_stddev_dict


def minmax_feature_scaling(fc, scaling_dict, feature_names):
    """Function to apply min/max scaling to feature collection

    args:
        fc (ee.FeatureCollection): feature collection to scale
        scaling_dict (ee.Dictionary): dictionary of min/max values to scale to
        feature_names (list[str]):  names of feature columns to calculate apply scaling to

    returns:
        ee.FeatureCollection: scaled feature collection
    """

    def feature_scaling(feature):
        """Nested closure function to apply scaling on each column in each feature"""

        def iter_cols(i):
            """Loops through feature columns"""
            i = ee.String(i)
            v = ee.Number(feature.get(i))
            minv = ee.Number(scaling_dict.get(i.cat("_min")))
            maxv = ee.Number(scaling_dict.get(i.cat("_max")))
            return v.subtract(minv).divide(maxv.subtract(minv))

        # apply scaling on each column of feature
        scaled = ee_feature_names.map(iter_cols)
        # get a dictionary of new values with old feature names
        newVals = ee.Dictionary.fromLists(ee_feature_names, scaled)
        # set feature columns new values
        return feature.set(newVals)

    # force ee types
    fc = ee.FeatureCollection(fc)
    ee_feature_names = ee.List(feature_names)

    # normalize the features in the entire featureCollection
    fc_norm = fc.map(feature_scaling)

    return fc_norm


@decorators.keep_attrs
def minmax_image_scaling(image, scaling_dict, feature_names):
    """Function to scale image between min/max values
    Expects that scaling_dict keys match bands

    args:
        image (ee.Image): image to scale
        scaling_dict (ee.Dictionary): dictionary of min/max values to scale to

    returns
        ee.Image: scaled image
    """
    # get dict as image
    scaling_img = scaling_dict.toImage()
    # extract the min/max values per band
    min_img = scaling_img.select(".*_min")
    max_img = scaling_img.select(".*_max")
    # apply scaling
    return (
        image.select(sorted(feature_names))
        .subtract(min_img)
        .divide(max_img.subtract(min_img))
        .float()
    )


def standard_feature_scaling(fc, scaling_dict, feature_names):
    """Function to apply standard (Z-score) scaling to feature collection

    args:
        fc (ee.FeatureCollection): feature collection to scale
        scaling_dict (ee.Dictionary): dictionary of mean/std dev values for scaling
        feature_names (list[str]):  names of feature columns to calculate apply scaling to

    returns:
        ee.FeatureCollection: scaled feature collection
    """

    def feature_scaling(feature):
        """Nested closure function to apply scaling on each column in each feature"""

        def iter_cols(i):
            """Loops through feature columns"""
            i = ee.String(i)
            v = ee.Number(feature.get(i))
            mean = ee.Number(scaling_dict.get(i.cat("_mean")))
            stddev = ee.Number(scaling_dict.get(i.cat("_stdDev")))
            return v.subtract(mean).divide(stddev)

        # apply scaling on each column of feature
        scaled = ee_feature_names.map(iter_cols)
        # get a dictionary of new values with old feature names
        newVals = ee.Dictionary.fromLists(ee_feature_names, scaled)
        # set feature columns new values
        return feature.set(newVals)

    # force ee types
    fc = ee.FeatureCollection(fc)
    ee_feature_names = ee.List(feature_names)

    # normalize the features in the entire featureCollection
    fc_norm = fc.map(feature_scaling)

    return fc_norm


@decorators.keep_attrs
def standard_image_scaling(image, scaling_dict, feature_names):
    """Function to apply z-score scaling to image
    Expects that scaling_dict keys match bands

    args:
        image (ee.Image): image to scale
        scaling_dict (ee.Dictionary): dictionary of mean/std dev values to scale to

    returns
        ee.Image: scaled image
    """
    # get dict as image
    scaling_img = scaling_dict.toImage()
    # extract the min/max values per band
    mean_img = scaling_img.select(".*_mean")
    stddev_img = scaling_img.select(".*_stdDev")
    # apply scaling
    return (
        image.select(sorted(feature_names))
        .subtract(mean_img)
        .divide(stddev_img)
        .float()
    )


def onehot_feature_encoding(fc, column_name, classes, class_names=None):
    """Function to calculate one-hot encoded columns from categorial columns
    where each new column equals 1 where the class value is the column index

    args:
        fc (ee.FeatureCollection): Feature collection with categorial data to encode
        column_name (str | ee.String): name of column that is categorial to encode
        classes (list[int] | ee.List): list of class values to encode

    kwargs:
        class_names (ee.List, optional): list of names to rename output bands.
            if None then bands will be named b0, b1, ...,bn. default = None

    returns:
        ee.Image: Feature collection with one-hot encoded columns with n new column as n classes
    """

    def feature_encoding(feature):
        c = ee.Number(feature.get(column_name))
        encoded = classes.map(lambda x: c.eq(ee.Number(x)))
        new_cols = ee.Dictionary.fromLists(class_names, encoded)
        return feature.set(new_cols)

    if class_names is None:
        class_names = ee.List.sequence(0, classes.length()).map(
            lambda x: ee.String("b").cat(ee.String(x))
        )

    fc_encoded = fc.map(feature_encoding)

    return fc_encoded


@decorators.keep_attrs
def onehot_image_encoding(img, classes, class_names=None, band=None):
    """Function to convert an categorial image image to one-hot encoded image
    where each new band equals 1 where the class value is the band index

    args:
        img (ee.Image): categorical image to encode
        classes (list[int] | ee.List): list of class values to encode

    kwargs:
        class_names (ee.List, optional): list of names to rename output bands.
            if None then bands will be named b0, b1, ...,bn. default = None
        band (str | ee.String, optional): name of band from input image to endcode.
            if None then the first band is used. default = None

    returns:
        ee.Image: one-hot encoded image with n bands as n classes
    """

    classes = ee.List(classes)

    if class_names is None:
        class_names = ee.List.sequence(0, classes.length()).map(
            lambda x: ee.String("b").cat(ee.String(x))
        )

    if band is None:
        img = img.select([0])
    else:
        img = img.select(band)

    encoded_imgs = classes.map(lambda x: img.eq(ee.Number(x)))

    return (
        ee.ImageCollection.fromImages(ee.List(encoded_imgs))
        .toBands()
        .rename(class_names)
    )


def random_forest_ee(
    n_trees,
    feature_collection,
    feature_names,
    label,
    scaling=None,
    mode="classification",
    min_samples_leaf=1,
):
    """Helper function to scale feature collection and train random forest model

    args:
        n_trees (int): number of trees for random forest model
        feature_collection (ee.FeatureCollection): features to train random forest model
        feature_names (list[str]): names of feature columns to use in random forest model (x values)
        label (str): name of feature column to fit random forest model (y value)
        scaling (str | None, optional): name of scaling to apply before training. One of: "minmax", "standard", `None`.
            default = `None`
        mode (str, optional): The output mode of the random forest model. One of: "classification", "regression",
            "probability". default = "classification"
    """

    if scaling == "minmax":
        scaling_dict = minmax_scaling_dict(feature_collection, feature_names)
        fc_norm = minmax_feature_scaling(
            feature_collection, scaling_dict, feature_names
        )

    elif scaling == "standard":
        scaling_dict = standard_scaling_dict(feature_collection, feature_names)
        fc_norm = standard_feature_scaling(
            feature_collection, scaling_dict, feature_names
        )

    elif scaling is None:
        scaling_dict = None
        fc_norm = feature_collection

    else:
        raise ValueError(
            "Could not determine scaling option. Options are ['minmax', 'standard', or None]"
        )

    classifier = (
        ee.Classifier.smileRandomForest(n_trees, minLeafPopulation=min_samples_leaf)
        .setOutputMode(mode.upper())
        .train(fc_norm, label, feature_names)
    )

    return classifier, scaling_dict


def gradient_boosting_ee(
    n_trees,
    feature_collection,
    feature_names,
    label,
    scaling=None,
    mode="classification",
    shrinkage=0.01,
    loss="LeastAbsoluteDeviation",
):
    """Helper function to scale feature collection and train gradient tree boosting model

    args:
        n_trees (int): number of trees for gradient boosting model
        feature_collection (ee.FeatureCollection): features to train random forest model
        feature_names (list[str]): names of feature columns to use in random forest model (x values)
        label (str): name of feature column to fit random forest model (y value)
        scaling (str | None, optional): name of scaling to apply before training. One of: "minmax", "standard", `None`.
            default = `None`
        mode (str, optional): The output mode of the random forest model. One of: "classification", "regression",
            "probability". default = "classification"
        learning_rate (float,optional): The shrinkage parameter in (0, 1] controls the learning rate of procedure. default = 0.01
        loss (str, optional): Loss function to be optimized. default = "LeastAbsoluteDeviation"
    """

    if scaling == "minmax":
        scaling_dict = minmax_scaling_dict(feature_collection, feature_names)
        fc_norm = minmax_feature_scaling(
            feature_collection, scaling_dict, feature_names
        )

    elif scaling == "standard":
        scaling_dict = standard_scaling_dict(feature_collection, feature_names)
        fc_norm = standard_feature_scaling(
            feature_collection, scaling_dict, feature_names
        )

    elif scaling is None:
        scaling_dict = None
        fc_norm = feature_collection

    else:
        raise ValueError(
            "Could not determine scaling option. Options are ['minmax', 'standard', or None]"
        )

    classifier = (
        ee.Classifier.smileGradientTreeBoost(
            numberOfTrees=n_trees, shrinkage=shrinkage, loss=loss
        )
        .setOutputMode(mode.upper())
        .train(fc_norm, label, feature_names)
    )

    return classifier, scaling_dict


def unsupervised_rf(
    n_trees,
    samples,
    features=None,
    rank_feature=None,
    ranking="min",
):
    """Unserpersived machine learning workflow to classify water
    Methods similar to: https://doi.org/10.1016/j.rse.2020.112209

    args:
        n_trees (int): number of trees to creat random forest model for class generalization
        samples (ee.FeatureCollection): input samples to create water classifier for
        features (list | ee.List): property names from samples to use for the semi supervised classification, If none then all properties are used. default = None
        rank_feature (str, optional): property name used to rank which unserpervised class is water. If None then first band name in `bands` is used. default = None
        ranking (str, optional): method to rank the classes by `rank_band`. Options are 'min' or 'max'. If 'min', then the lowest class mean is considered water. default = 'min'

    returns:
        ee.Classifier.RandomForest: random forest classifier to estimate probability that a pixel is water
    """

    def _cluster_center(x):
        return (
            feature_arr.mask(classes.eq(ee.Number(x))).reduce(ee.Reducer.mean(), [0])
        ).get([0])

    if features is None:
        bands = samples.first().propertyNames()
    else:
        features = ee.List(features)

    if rank_feature is None:
        rank_feature = ee.String(features.get(0))

    clusterer = ee.Clusterer.wekaXMeans(3, 12, 5).train(samples, features)

    samples = samples.cluster(clusterer, "init_classes")

    classes = samples.aggregate_array("init_classes")
    unique = classes.distinct().sort()
    classes = ee.Array(classes)
    feature_arr = ee.Array(samples.aggregate_array(rank_feature))

    class_means = unique.map(_cluster_center)

    if ranking == "min":
        ranker = ee.Reducer.min()
    elif ranking == "max":
        ranker = ee.Reducer.max()
    else:
        raise NotImplementedError(
            "ranking selection is not implemented. options are 'min' or 'max'"
        )

    ranked_mean = class_means.reduce(ranker)

    water_class = class_means.indexOf(ranked_mean)

    binary_samples = samples.map(
        lambda x: (
            ee.Feature(x).set(
                "init_classes", ee.Number(x.get("init_classes")).eq(water_class)
            )
        )
    )

    classifier = (
        ee.Classifier.smileRandomForest(numberOfTrees=n_trees)
        .setOutputMode("PROBABILITY")
        .train(binary_samples, "init_classes", features)
    )

    return classifier


def calc_image_pca(image, region=None, scale=90, max_pixels=1e9, method="svd"):
    """Principal component analysis decomposition of image bands

    args:
        image (ee.Image): image to apply pca to
        region (ee.Geometry | None, optional): region to sample values for covariance matrix,
            if set to `None` will use img.geometry(). default = None
        scale (int, optional): scale at which to perform reduction operations, setting higher will prevent OOM errors. default = 90
        max_pixels (int, optional): maximum number of pixels to use in reduction operations. default = 1e9
        method (str, optional): the decomposition method for obtaining the eigen vectors and values
            options are 'svd' or 'eigendecomp'. note: svd is usually faster as is does not need to
            compute the covariance matrix of input features. default = 'svd'

    returns:
        ee.Image: principal components scaled by eigen values
    """

    bandNames = image.bandNames()

    out_band_names = ee.List.sequence(1, bandNames.length()).map(
        lambda x: ee.String("pc_").cat(ee.Number(x).int())
    )

    # Mean center the data to enable a faster covariance reducer
    # and an SD stretch of the principal components.
    meanDict = image.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=region, scale=scale, maxPixels=max_pixels
    )
    means = ee.Image.constant(meanDict.values(bandNames))
    centered = image.subtract(means)

    # Collapse the bands of the image into a 1D array per pixel.
    arrays = centered.toArray()

    if method == "svd":
        svd = arrays.toArray(1).matrixSingularValueDecomposition()

        eigen_vecs = svd.select("V")  # .arrayTranspose()
        eigen_vals = svd.select("S").matrixDiagonal()

    elif method == "eigendecomp":
        # Compute the covariance of the bands within the region.
        covar = arrays.reduceRegion(
            reducer=ee.Reducer.centeredCovariance(),
            geometry=region,
            scale=scale,
            maxPixels=max_pixels,
        )

        # Get the 'array' covariance result and cast to an array.
        # This represents the band-to-band covariance within the region.
        covarArray = ee.Array(covar.get("array"))

        # Perform an eigen analysis and slice apart the values and vectors.
        eigens = covarArray.eigen()

        # This is a P-length vector of Eigenvalues.
        eigenValues = eigens.slice(1, 0, 1)
        # This is a PxP matrix with eigenvectors in rows.
        eigenVectors = eigens.slice(1, 1)

    # Convert the array image to 2D arrays for matrix computations.
    arrayImage = arrays.toArray(1)

    # Left multiply the image array by the matrix of eigenvectors.
    principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage)

    # Turn the square roots of the Eigenvalues into a P-band image.
    sdImage = (
        ee.Image(eigenValues.sqrt()).arrayProject([0]).arrayFlatten([out_band_names])
    )

    # Turn the PCs into a P-band image, normalized by SD.
    return (
        principalComponents
        # Throw out an an unneeded dimension, [[]] -> [].
        .arrayProject([0])
        # Make the one band array image a multi-band image, [] -> image.
        .arrayFlatten([out_band_names])
        # Normalize the PCs by their SDs.
        .divide(sdImage)
    )


def calc_feature_pca(fc, names, is_centered=False, method="svd"):
    """Principal component decomposition of features

    args:
        fc (ee.FeatureCollection): feature collection to caluculate PCA from
        names (list[str]): property names to uses as features in PCA
        is_centered (bool, optional): boolean to identify if features need to be centered before PCA.
            False means apply centering. default = False
        method (str, optional): the decomposition method for obtaining the eigen vectors and values
            options are 'svd' or 'eigendecomp'. note: svd is usually faster as is does not need to
            compute the covariance matrix of input features. default = 'svd'

    returns:
        ee.Array: eigen vectors of PCA
        ee.Array: eigen values of PCA
        ee.Array: mean values of each feature
    """
    array_ = ee.Array(
        ee.List(
            fc.makeArray(names)
            .aggregate_array("array")
            .map(lambda x: ee.Array(x).toList())
        )
    )
    center = array_.reduce(ee.Reducer.mean(), [0]).repeat(0, array_.length().get([0]))
    if not is_centered:
        centered = array_.subtract(center)
    else:
        centered = array_

    if method == "svd":
        svd = centered.matrixSingularValueDecomposition()

        eigen_vecs = ee.Array(svd.get("V")).transpose()
        eigen_vals = ee.Array(svd.get("S")).matrixDiagonal()

    elif method == "eigendecomp":
        # Compute the covariance of the bands within the region.
        covar = centered.transpose().matrixMultiply(centered)
        # Perform an eigen analysis and slice apart the values and vectors.
        eigens = covar.eigen()

        eigen_vecs = eigens.slice(1, 1)
        eigen_vals = eigens.slice(1, 0, 1)

    else:
        raise ValueError(
            "could not understand provided method keyword. Options are 'svd' or 'eigendecomp'"
        )

    out_band_names = [f"pc_{i}" for i in range(len(names))]

    return eigen_vecs, eigen_vals, center.slice(0, 0, 1).project([1])


def apply_feature_pca(fc, eigen_vecs, names, center=None):
    """Applies Principal component decomposition on features

    args:
        fc (ee.FeatureCollection): feature collection to caluculate pricipal components from
        eigen_vecs (ee.Array): eigen vectors of PCA to transform features
        names (list[str]): property names to uses as features in PCA
        center (ee.Array | None, optional): Array of mean values to center features. If None then no
            centering is applies. default = None

    returns:
        ee.FeatureCollection: feacture collection with new properties within each feature being the principal components
    """
    array_ = ee.Array(
        ee.List(
            fc.makeArray(names)
            .aggregate_array("array")
            .map(lambda x: ee.Array(x).toList())
        )
    )
    if center is not None:
        centered = array_.subtract(
            ee.Array.cat([center], 1).transpose().repeat(0, array_.length().get([0]))
        )
    else:
        centered = array_

    pca_arr = eigen_vecs.matrixMultiply(centered.transpose()).transpose()

    out_band_names = [f"pc_{i}" for i in range(len(names))]

    fc_size = fc.size()
    fc_list = fc.toList(fc_size)
    fc_pca = ee.FeatureCollection(
        ee.List.sequence(0, fc_size.subtract(1)).map(
            lambda x: ee.Feature(fc_list.get(x)).set(
                ee.Dictionary.fromLists(
                    out_band_names,
                    pca_arr.slice(0, x, ee.Number(x).add(1), 1).project([1]).toList(),
                )
            )
        )
    )

    return fc_pca


@decorators.keep_attrs
def apply_image_pca(img, eigen_vecs, names, center=None):
    """Applies Principal component decomposition on image

    args:
        img (ee.Image): image to caluculate pricipal components from
        eigen_vecs (ee.Array): eigen vectors of PCA to transform features
        names (list[str]): band names to uses as features in PCA
        center (ee.Array | None, optional): Array of mean values to center features. If None then no
            centering is applies. default = None

    returns:
        ee.Image: principal components calculated from image
    """
    if center is not None:
        arrayImage = (
            img.select(names)
            .subtract(ee.Image.constant(center.toList()))
            .toArray()
            .toArray(1)
        )
    else:
        arrayImage = img.select(names).toArray().toArray(1)

    principalComponents = ee.Image(eigen_vecs).matrixMultiply(arrayImage)

    out_band_names = [f"pc_{i}" for i in range(len(names))]

    pcaImage = (
        principalComponents
        # Throw out an an unneeded dimension, [[]] -> [].
        .arrayProject([0])
        # Make the one band array image a multi-band image, [] -> image.
        .arrayFlatten([out_band_names])
    )

    return pcaImage


def hist_matching(samples, predictor, target, n_estimators=50):
    """Trains classifiers to perform histogram matching

    args:
        samples (ee.FeatureCollection): feature collection with samples for histogram matching
        predictor (str): column name of values to transform
        target (str): column name of values to match
        n_estimators (int, optional): number of trees to create random forest models from. default = 50

    returns:
        list[ee.Classifier]: list of classifiers with first element being the val to proba and second being proba to val classifiers

    """

    def get_cdf(fc, column):
        def array_to_features(l):
            return ee.Feature(
                None, {column: ee.List(l).get(0), "probability": ee.List(l).get(1)}
            )

        # Histogram equalization start:
        histo = ee.Dictionary(
            fc.reduceColumns(
                ee.Reducer.histogram(
                    maxBuckets=2 ** 12,
                ),
                [column],
            ).get("histogram")
        )

        valsList = ee.List(histo.get("bucketMeans"))
        freqsList = ee.List(histo.get("histogram"))
        cdfArray = ee.Array(freqsList).accum(0)
        total = cdfArray.get([-1])
        normalizedCdf = cdfArray.divide(total)

        array = ee.Array.cat([valsList, normalizedCdf], 1)

        return ee.FeatureCollection(array.toList().map(array_to_features))

    pred_cdf = get_cdf(samples, predictor)
    target_cdf = get_cdf(samples, target)

    proba_to_val = (
        ee.Classifier.smileRandomForest(n_estimators)
        .setOutputMode("REGRESSION")
        .train(
            features=target_cdf, classProperty=target, inputProperties=["probability"]
        )
    )

    val_to_proba = (
        ee.Classifier.smileRandomForest(n_estimators)
        .setOutputMode("REGRESSION")
        .train(
            features=pred_cdf, classProperty="probability", inputProperties=[predictor]
        )
    )

    return val_to_proba, proba_to_val


@decorators.keep_attrs
def apply_image_matching(image, matching_classifiers, output_name="dn"):
    return image.classify(matching_classifiers[0], "probability").classify(
        matching_classifiers[1], output_name
    )
