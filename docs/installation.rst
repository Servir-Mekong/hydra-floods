.. _installation:

Installation
============

Instructions
------------

hydrafloods itself is a pure Python package, but its dependencies are not. The
easiest way to get most dependencies installed is to use conda_. It is recommend
using the community maintained `conda-forge <https://conda-forge.github.io/>`_
channel to handle difficult\-to\-build dependencies.Furthermore, it is good
practice to use a virtual environment within conda. To create a new environment,
install dependencies, and activate the environment::

  $ conda create -n hydra -c conda-forge python=3.7 numpy scipy pandas requests yaml xmltodict gdal shapely pyproj netCDF4 xarray pyresample geopandas earthengine-api fire
  $ conda activate hydra

.. _conda: http://conda.io/

We will now also need to install the `Google Storage Utility command-line interface (gsutil) <https://cloud.google.com/storage/docs/gsutil_install>`_.
Follow the directions provided by the website. After successful installation,
the ``gcloud`` environment will need to be initialized::

  $ gcloud init

hydrafloods also relies on using `Google Earth Engine <https://earthengine.google.com>`_
so we will now have to authenticate the EE Python API::

  $ earthengine authenticate


.. warning::

    Make sure you initialize the ``earthengine`` and ``gcloud`` APIs with Google
    accounts that have permissions to read and write to Google Cloud Storage and
    Google Earth Engine assets.
