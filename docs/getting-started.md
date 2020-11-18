Here are some quick examples of what you can do with `hydrafloods`. It is expected that the code is run in an interactive python session such as IPython or in a Jupyter Notebook as later code blocks will use variables from previous ones.

To get started, first import the `ee`, trigger the Earth Engine authentication flow, and import the `hydrafloods` package:

```python
import ee
ee.Initialize()
import hydrafloods as hf
```

## Get a `hf.Dataset`

You can access commonly used image collections on Earth Engine as a [`hydrafloods.Dataset`](/datasets/) to quickly filter by space and time as well as apply pre-written QA masking functions.

```python
# define a geographic region
region = hf.country_bbox("Cambodia")

# define start and end times
start_time = "2019-09-15"
end_time = "2019-09-20"

# get the Sentinel 1 collection as a Dataset
s1 = hf.Sentinel1(region,start_time,end_time)

# print dataset info
print(s1)

# print the number of images in Dataset
print(s1.n_images)

# print a list of image acquisition dates for the Dataset
print(s1.dates)

# access the actual ee.ImageCollection object
ee_collection = s1.collection
```

The `hydrafloods.Dataset` object is a wrapper around an `ee.ImageCollection` by applying the spatial and temporal filtering upon initialization.
This provides a quick and consistent access to imagery. The `Dataset` class also provides utility functionality to make working with and managing multiple image collections less verbose.

There are many ways to interface with datasets (i.e. ImageCollections) using `hydrafloods`, more examples on merging or joining datasets can be found on the [Using Dataset class](/using-datasets/) page.


## Image processing

The main purpose of `hydrafloods` is to lower the barrier to creating high-quality surface water maps, this requires image processing. Although the Dataset class wraps an Earth Engine image collection we can apply image processing functions using [`apply_func()`](/datasets/#hydrafloods.datasets.Dataset.apply_func) by passing a function object. 

This method wraps a function that accepts an image as the first argument (which most `hydrafloods` image processing algorithms do) and maps it over the collection. For example, we would like to apply a speckle filter algorithm on SAR imagery. We can easily do this with the following code.

```python
# apply a speckle filtering algorithm on SAR imagery
# here we will use the Gamma Map filter
filtered = s1.apply_func(hf.gamma_map)
```

The previous example is synonymous with using `s1.collection = s1.collection.map(hf.gamma_map)` which access the image collection, applies the function, and sets the results to the `s1.collection` property. Although this technically works, using the `apply_func()` method is advantageous and preferred as it allows us to pass arbitrary keyword parameters to functions which we want to apply. For example, the water mapping algorithms found in [`hydrafloods.thresholding`](/thresholding/) take many keyword parameters and we can customize function as in the following example.

```python
# apply the edge otsu surface water mapping 
# we apply this on the speckle filtered SAR data
water_maps = filtered.apply_func(hf.edge_otsu, 
    initial_threshold=-16,
    edge_buffer=300,
    scale=250
)
```

It should be noted that using the `apply_func()` method will return a `hydrafloods.Dataset` where the collection property has the results of the function. One can access the `ee.ImageCollection` object and reduce to and image using the following code:

```python
# reduce the Dataset.collection property to an ee.Image
water_img = water_maps.collection.reduce("mode")
```

There are a variety of image processing functions available in `hydrafloods`, more information on specific algorithms can be found on the [Algorithms](/algorithms/) page.


## Time series processing

In addition to image processing, processing data in time is valuable. Therefore, `hydrafloods` has a specific module for time series processing, [`hydrafloods.timeseries`](/timeseries/), specifically for processing stacks of imagery in time.

```python
# import in the timeseries module
from hydrafloods import timeseries
```

Here we are going to take a longer time series of SAR imagery for 2019 so we have more data for our model:

```python
# define start and end times for one year in 2019
start_time = "2019-01-01"
end_time = "2020-01-01"

# get the Sentinel 1 collection as a Dataset
# using Cambodia still as the region
s1 = hf.Sentinel1(region,start_time,end_time)
```

Now that we have our time series of data, we can begin to model some temporal information. For example, we want to model a harmonic trend of SAR imagery to predict an image using time information.

```python
# fit a harmonic trend model on the VV band in time
# use two harmonic cycles
harmonics_weights = timeseries.fit_harmonic_trend(
    s1, dependent='VV', n_cycles=2
)
```

The result from [`fit_harmonic_trend()`](/timeseries/#hydrafloods.timeseries.fit_harmonic_trend) will be an image with many bands. Some bands are the coefficeint weights for prediction (i.e. cos_1, sin_1, or time), others can be awareness information such as number of valid observations used (i.e. n). So we will filter out the coefficient weight bands we need which are cos_n, sin_n and time which start with either "c", "t", or "s". Then get a dummy image with time information and apply the prediction.

```python
# extract bands needed for prediction
harmonic_weights = harmonics.select("^(c|t|s).*")

# get a dummy image with just time information for prediction
# for flooding date in Oct 2019
dummy_img = timeseries.get_dummy_img("2019-10-05")

# predict the VV values using the dummy img and computed harmonic coeffients
prediction = (
    timeseries.add_harmonic_coefs(dummy_img)
    .multiply(harmonic_weights)
    .reduce("sum")
)
```

Time series functionality in `hydrafloods` is focused around modeling data in time, more information on the functions can be found in the [timeseries module API reference](/timeseries/)


## Machine Learning

`hydrafloods` also has a specific module for machine learning workflows with Earth Engine, [`hydrafloods.ml`](/ml/). 

```python
# import in the ml module
from hydrafloods import ml
```

The aim with this module is to make high-quality machine learning workflows easier and less verbose when working with Earth Engine. However, we will will still need to access feature collection information to train models. Here we will define some parameters for an example:

```python
# define some parameters for the ml workflow we will use
# define feature columns to use for training/prediction
feature_names = ['VV','VH','ratio','ndpi'] 

# get a feature collection for training/testing
fc = (
    ee.FeatureCollection("projects/servir-ee/assets/sar_wi_samples_20200825104400")
    .randomColumn("random") 
)

# split the feature collection into training/testing datasets
# good practice to do this
training = fc.filter(ee.Filter.lte("random",0.7)) 
testing = fc.filter(ee.Filter.gt("random",0.7))
```

Now that we have our datasets and features defined, we can begin the workflow. Typically with machine learning one would perform a feature scaling to get all features within the same range, this helps with numerical stability, training speed, and model accuracy. This can be quite cumbersome to do efficiently in Earth Engine...however, the `ml` module allows us to scale and train in one shot.

```python
# scale training dataset and train random forest model
# note this returns a trained ee.Classifier and dictionary to scale data later
rf, scaling_dict = ml.random_forest_ee(
    25,
    training,
    feature_names,
    'mndwi',
    scaling="standard",
    mode="regression"
)
```

This function returns a train random forest `ee.Classifier` and a `ee.Dictionary` with per band values needed to scale other feature collections or imagery as seen in the next example. Here we scale the testing dataset using the values from the training dataset and apply the model on the feature collection.

```python
# scale the testing dataset using the scaling values from training
testing_norm = ml.standard_feature_scaling(testing,scaling_dict,feature_names)

# apply random forest model on scaled test feature collection dataset
y_test = testing_norm.classify(rf,"predicted")
```

Majority of the time we would like to apply the predictions on imagery and the `ml` module has functionality to perfrom the scaling for imagery easily. First we have to add the bands to dataset collection which the `Sentinel1` dataset class has a custom method to add the VV/VH ratio and VV-VH normalized difference bands.

```python
# add bands for features used in RF model
s1_features = s1.add_fusion_features()

# scale the bands using the scaling_dict
s1_norm = s1_features.apply_func(
    ml.standard_image_scaling, 
    scaling_dict=scaling_dict,
)

# apply the prediction on the dataset
predicted = s1_norm.apply_func(lambda x: x.classify(rf))
```

Again, most of the functionality around the `hydrafloods.ml` module is to make end-to-end machine learning work flows more straightforward. Please see the [ml module](/ml/) documentation for information on functions.
