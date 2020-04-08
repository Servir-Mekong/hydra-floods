# hydra-floods

## Introduction
HYDrologic Remote sensing Analysis for Floods (HydraFloods) is a Python application for downloading, processing, and delivering flood maps derived from remote sensing data.
This is backend behind the [HYDRAViewer](https://github.com/Servir-Mekong/hydrafloodviewer) Geospatial tool which utilizes many remote sesnsing dataset as possible to provide daily (sometime twice daily) flood maps.

![hydrafloodviewer](https://user-images.githubusercontent.com/39947610/65659571-5044b480-e056-11e9-8b1b-d3e82cace3c8.PNG)

[![IMAGE ALT TEXT HERE](https://i.ytimg.com/vi/Iqh_uwtg6yI/hqdefault.jpg)](https://www.youtube.com/watch?v=Iqh_uwtg6yI&feature=youtu.be)

### Installation

1. Create a new [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) with the required dependencies and activate the environment

```
$ conda create -n hydra -c conda-forge python=3.6 numpy scipy pandas requests yaml xmltodict gdal shapely pyproj netCDF4 xarray pyresample geopandas earthengine-api fire

$ conda activate hydra
```

2. Install additional packages via `pip` that are not on conda-forge

```
$ pip install git+git://github.com/KMarkert/simple-cmr.git
```

3. Authenticate the Earth Engine API

```
$ earthengine authenticate
```

4. Install [Google Storage Utility command-line interface (gsutil)](https://cloud.google.com/storage/docs/gsutil_install) and the setup the `gsutil` environment

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

To initiate the HydraFloods process use the command below.

```
$hydrafloods process atms 2015-07-19 --configuration myconf.yaml
```
**NOTE:**
Details regarding the command

* hydrafloods      --- call the hydrafloods process to begin

* process           --- start the hydraflood process

* atms              --- specify sensor (see processing.py for details options include: Sentinel2, Landsat, Modis, Viirs, Atms, Sentinel1)

* 2015-07-19       --- specify the date of intrest

* --configuration  --- use the configuration file

* myconf.yaml      --- call the .yaml file that specifies the ROI etc. (see the exampl.yaml for details)


<img src="exampleyaml.JPG" width="50%">
