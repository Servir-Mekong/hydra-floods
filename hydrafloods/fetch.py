import os
import datetime
import requests
import simplecmr as scmr
from pathlib import Path

# requires simplecmr to be installed


def fetching(
    conceptid,
    start_time,
    region,
    credentials,
    out_directory,
    max_results=500,
    end_time=None,
):
    """ Function to download data from NASA CMR by specifying a dataset, time and region.
    Uses CMR to handle spatio-temporal query and extracts data download urls

    args:
        conceptid (str): String of dataset concept id to search for
        start_time (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
        region (tuple[float] | list[float]): Bounding box of region to search as iterable in W,S,E,N order
        credentials (tuple | list): EarthData username and password login credentials as iterable
        out_directory (str|pathlib.Path): Local directory to downaload data to
        max_results (int, optional): Maximum number of items to search and download.
            default = 500
        end_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd).
            default = start_time + 1day

    returns:
        list[pathlib.Path]: List of local paths that data was downloaded to
    """
    if type(start_time) != datetime.datetime:
        start_time = scmr.utils.decode_date(start_time)
    if end_time is None:
        end_time = start_time + datetime.timedelta(seconds=86399)
    else:
        end_time = scmr.utils.decode_date(end_time)

    # construct query
    query = scmr.Query(
        conceptid=conceptid,
        startTime=start_time,
        endTime=end_time,
        spatialExtent=region,
        maxResults=max_results,
    )

    # fetch datasets from query
    query.granules.fetch(
        credentials=credentials,
        directory=out_directory,
        limit=max_results,
        maxWorkers=4,
    )

    # return a list of the granules for later processing
    return query.granules.getLocalPaths(directory=out_directory)


def viirs(
    credentials,
    start_time="2000-01-01",
    end_time=None,
    region=[-180, 60, 180, 85],
    out_directory="./",
):
    """Function to download Suomi-NPP VIIRS surface reflectance data for specified time and region,
    wraps `fetching()`

    args:
        credentials (tuple | list): EarthData username and password login credentials as iterable
        start_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        end_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = start_time + 1day
        region (tuple | list, optional): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        out_directory (str | pathlib.Path, optioanl): Local directory to downaload data to
            default = "./"

    returns:
        list[pathlib.Path]: List of local paths that data was downloaded to
    """

    # check if date requested is within science-quality production time
    now = datetime.datetime.now()
    offset = now - scmr.utils.decode_date(start_time)
    if offset.days > 4:
        # use science-quality collection
        CONCEPTID = "C1373412034-LPDAAC_ECS"
    else:
        # use LANCE-NRT collection
        CONCEPTID = "C1344293643-LANCEMODIS"

    return fetching(
        conceptid=CONCEPTID,
        start_time=start_time,
        region=region,
        credentials=credentials,
        out_directory=out_directory,
        max_results=500,
        end_time=end_time,
    )


def modis(
    credentials,
    start_time="2000-01-01",
    end_time=None,
    region=[-180, 60, 180, 85],
    out_directory="./",
):
    """Function to download MODIS surface reflectance data for specified time and region,
    wraps `fetching()`

    args:
        credentials (tuple | list): EarthData username and password login credentials as iterable
        start_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        end_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = start_time + 1day
        region (tuple[float] | list[float], optional): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        out_directory (str | pathlib.Path, optional): Local directory to downaload data to
            default = "./"

    returns:
        list[pathlib.Path]: List of local paths that data was downloaded to
    """

    # check if date requested is within science-quality production time
    now = datetime.datetime.now()
    offset = now - scmr.utils.decode_date(start_time)
    if offset.days > 4:
        # use science-quality collection
        CONCEPTID = "C193529902-LPDAAC_ECS"
    else:
        # use LANCE-NRT collection
        CONCEPTID = "C1219249711-LANCEMODIS"

    return fetching(
        conceptid=CONCEPTID,
        start_time=start_time,
        region=region,
        credentials=credentials,
        out_directory=out_directory,
        max_results=500,
        end_time=end_time,
    )


def atms(
    credentials,
    start_time="2000-01-01",
    end_time=None,
    region=[-180, 60, 180, 85],
    out_directory="./",
):
    """ Function to download Suomi-NPP ATMS passive microwave data for specified time and region,
    wraps `fetching()`

    args:
        credentials (tuple | list): EarthData username and password login credentials as iterable
        start_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        end_time (str, optional): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = start_time + 1day
        region (tuple | list, optional): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        out_directory (str | pathlib.Path, optioanl): Local directory to downaload data to
            default = "./"

    returns:
        list[pathlib.Path]: List of local paths that data was downloaded to
    """

    CONCEPTID = "C1442068516-GES_DISC"
    return fetching(
        conceptid=CONCEPTID,
        start_time=start_time,
        region=region,
        credentials=credentials,
        out_directory=out_directory,
        max_results=500,
        end_time=end_time,
    )
