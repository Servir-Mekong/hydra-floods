# hydra-floods

## Introduction
HYDrologic Remote sensing Analysis for Floods (HydraFloods) is a Python application for downloading, processing, and delivering flood maps derived from remote sensing data.
The bases behind the tool is to use as many remote sesnsing dataset as possible to provide daily (sometime twice daily) flood maps.

### Installation

1. Create a new [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) with the required dependencies and activate the environment

```
$ conda create -n hydra -c conda-forge python=3.6 numpy scipy pandas requests yaml xmltodict gdal shapely pyproj netCDF4 xarray pyresample geopandas earthengine-api fire

$ conda activate hydra
```

2. Authenticate the Earth Engine API

```
$ earthengine authenticate
```

3. Install [Google Storage Utility command-line interface (gsutil)](https://cloud.google.com/storage/docs/gsutil_install) and the setup the `gsutil` environment

```
$ gcloud init
```

---
**NOTE:**

Make sure your initialize the `ee` and `gsutil` APIs with Google accounts that have permissions to read and write to Google Cloud Storage and Google Earth Engine assets.

---

4. Download the hydrafloods source code and install the package
 - COMING SOON: the hydrafloods package will be on PyPI in the near future to prevent installing from source

```
$ git clone https://github.com/servir-mekong/hydra-floods.git
$ cd hydra-floods/py pkg
$ python setup.py install
```

5. Test to see if the installation was successful

```
$ hydrafloods run_tests
```


### How-to
Coming soon...
