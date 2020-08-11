# Installation

`hydrafloods` itself is a pure Python package, but its dependencies are not. Furthermore, the package relies on Google Cloud and Google Earth Engine to stage and process data relying on specific software and even more importantly authentication. Here we describe two 

## Docker Image

The easiest way to get up and started using the `hydrafloods` packages is via a Docker Image. 


## Using Conda

Another convient way to install the package and its dependencies is using [Anaconda](https://www.anaconda.com/). It is recommend using the community maintained [conda-forge](https://conda-forge.github.io/) channel to handle difficult-to-build dependencies.urthermore, it is good practice to use a virtual environment within conda. To create a new environment,
install dependencies, and activate the environment:

```
  $ conda create -n hydra -c conda-forge python=3.7 \
    numpy \
    scipy \
    pandas \
    requests \
    yaml \
    xmltodict \
    gdal \
    shapely \
    pyproj \
    netCDF4 \
    xarray \
    pyresample \
    geopandas \
    earthengine-api \
    fire -y
  $ conda activate hydra
```


We will now also need to install the [Google Cloud SDK](https://cloud.google.com/sdk/docs/downloads-versioned-archives) to interface to with the Google cloud. Follow the directions provided by the website. 
After successful installation,
the `gcloud` environment will need to be initialized:

```
  $ gcloud init
```

`hydrafloods` also relies on using [Google Earth Engine](https://earthengine.google.com) for most of the geospatial processing so we will now have to authenticate the EE Python API:

```
  $ earthengine authenticate
```

<span style="color:red">
**Warning:** *Make sure you initialize the `earthengine` and `gcloud` APIs with Google accounts that have permissions to read and write to Google Cloud Storage and Google Earth Engine assets.*
</span>


## Testing Installation