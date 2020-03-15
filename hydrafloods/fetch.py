import os
import datetime
import requests
import simplecmr as scmr
from pathlib import Path

def fetching(conceptid,startTime,region,credentials,outDir,maxResults=500,endTime=None):
    """
    Function to download data from NASA by specifying a datase, time and region.
    Uses CMR to handle spatio-temporal query and extracts data download urls

    Args:
        conceptid (str): String of dataset concept id to search for
        startTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
        region (tuple|list): Bounding box of region to search as iterable in W,S,E,N order
        credentials (tuple|list): EarthData username and password login credentials as iterable
        outDir (str|pathlib.Path): Local directory to downaload data to
    Kwargs:
        maxResults (int): Maximum number of items to search and download
            default = 500
        endTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = None (endTime = startTime + 1day)

    Returns:
        List of local paths that data was downloaded to
    """
    
    startTime = scmr.utils.decode_date(startTime)
    if endTime is None:
        endTime = startTime + datetime.timedelta(seconds=86399)
    else:
        endTime = scmr.utils.decode_date(endTime)

    # construct query
    query = scmr.Query(conceptid=conceptid,
        startTime=startTime,
        endTime=endTime,
        spatialExtent=region,
        maxResults=maxResults,
    )
    print(query.granules)

    query.granules.fetch(credentials=credentials,
        directory=outDir,
        limit=maxResults,
        maxWorkers=4
    )

    return query.granules.getLocalPaths()


def viirs(credentials,startTime='2000-01-01',endTime=None,region=[-180,60,180,85],outDir='./'):
    """
    Function to download Suomi-NPP VIIRS surface reflectance data for specified time and region,
    wraps fetching()

    Args:
        credentials (tuple|list): EarthData username and password login credentials as iterable

    Kwargs:
        startTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        endTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = None (endTime = startTime + 1day)
        region (tuple|list): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        outDir (str|pathlib.Path): Local directory to downaload data to
            default = './' (current working directory)

    Returns:
        List of local paths that data was downloaded to
    """

    # check if date requested is within science-quality production time
    now = datetime.datetime.now()
    offset = now - scmr.utils.decode_date(startTime)
    if offset.days > 4:
        # use science-quality collection
        CONCEPTID = "C1373412034-LPDAAC_ECS"
    else:
        # use science-quality collection
        CONCEPTID = "C1344293643-LANCEMODIS"

    return fetching(conceptid=CONCEPTID,
        startTime=startTime,
        region=region,
        credentials=credentials,
        outDir=outDir,
        maxResults=500,
        endTime=endTime
    )

def modis(credentials,startTime='2000-01-01',endTime=None,region=[-180,60,180,85],outDir='./'):
    """
    Function to download Terra MODIS surface reflectance data for specified time and region,
    wraps fetching()

    Args:
        credentials (tuple|list): EarthData username and password login credentials as iterable

    Kwargs:
        startTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        endTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = None (endTime = startTime + 1day)
        region (tuple|list): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        outDir (str|pathlib.Path): Local directory to downaload data to
            default = './' (current working directory)

    Returns:
        List of local paths that data was downloaded to
    """

    # check if date requested is within science-quality production time
    now = datetime.datetime.now()
    offset = now - scmr.utils.decode_date(startTime)
    if offset.days > 4:
        # use science-quality collection
        CONCEPTID = "C193529902-LPDAAC_ECS"
    else:
        # use science-quality collection
        CONCEPTID = "C1219249711-LANCEMODIS"

    return fetching(conceptid=CONCEPTID,
        startTime=startTime,
        region=region,
        credentials=credentials,
        outDir=outDir,
        maxResults=500,
        endTime=endTime
    )


def atms(credentials,startTime='2000-01-01',endTime=None,region=[-180,60,180,85],outDir='./'):
    """
    Function to download Suomi-NPP ATMS passive microwave data for specified time and region,
    wraps fetching()

    Args:
        credentials (tuple|list): EarthData username and password login credentials as iterable

    Kwargs:
        startTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = 2000-01-01
        endTime (str): Date as string preferrably as ISO8601 format (YYYY-MM-dd)
            default = None (endTime = startTime + 1day)
        region (tuple|list): Bounding box of region to search as iterable in W,S,E,N order
            default = [-180,60,180,85]
        outDir (str|pathlib.Path): Local directory to downaload data to
            default = './' (current working directory)

    Returns:
        List of local paths that data was downloaded to
    """

    CONCEPTID = "C1442068516-GES_DISC"
    return fetching(conceptid=CONCEPTID,
        startTime=startTime,
        region=region,
        credentials=credentials,
        outDir=outDir,
        maxResults=500,
        endTime=endTime
    )
