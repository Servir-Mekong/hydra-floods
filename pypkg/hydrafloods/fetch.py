import os
import netrc
import requests
import datetime
import xmltodict
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

def viirs(date,h,v,outdir='./',creds=None,product=None):
    """Function to download VIIRS NRT data for specified time and tile

    Args:
        date (datetime.datetime): Datetime object specifying which date the data of interest was acquired.
        h (int): horizontal tile grid to fetch
        v (int): vertical tile grid to fetch
        outdir (str, optional): out directory to dump retrieved data to
        default = './' or current working directory
        creds (str, optional): path to .netrc file with NASA EarthData login in credentials
        default = None

    Returns:
        None
    """

    if outdir[-1] != '/':
        outdir = outdir+'/'

    acct = netrc.netrc(creds)
    usr,_,pswrd = acct.hosts['https://urs.earthdata.nasa.gov']

    today = datetime.datetime.now()

    if outdir[-1] != '/':
        outdir = outdir+'/'

    basename = '{0}.A{1}{2:03d}.h{3:02d}v{4:02d}.001.h5'

    if (today - date).days > 8:
        yr = date.year
        dt = (date-datetime.datetime(yr,1,1)).days + 1
        url = 'https://e4ftl01.cr.usgs.gov/DP102/VIIRS/{0}.001/{1}.{2:02d}.{3:02d}/'\
                .format(product,yr,date.month,date.day)
        with requests.Session() as s:
            s.auth = (usr, pswrd)

            r1 = s.request('get', url)
            r = s.get(r1.url, auth=(usr, pswrd))

            if r.ok:
                result = r.content.split(' '.encode())
                filtered = []
                for i in result:
                    if 'href'.encode() in i:
                        this = i.split('"'.encode())[1]
                        if this[-3:] =='.h5'.encode():
                            filtered.append(this.decode("utf-8"))

                for f in filtered:
                    if 'h{:02d}'.format(h) in f\
                    and 'v{:02d}'.format(v)  in f:
                        filename = basename.format(product,yr,dt,h,v)
                        outFile = outdir + filename
                        if os.path.exists(outFile) != True:
                            newurl = url+f
                            r3 = s.request('get', newurl)
                            r4 = s.get(r3.url, auth=(usr, pswrd))

                            with open(outFile, 'wb') as this:
                                this.write(r4.content) # Say

    else:
        yr = date.year
        dt = (date-datetime.datetime(yr,1,1)).days + 1

        url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/5000/{0}/Recent/'.format(product+'_NRT')
        filename = basename.format(product+'_NRT',yr,dt,h,v)
        fileUrl = url + filename

        outFile = outdir + filename
        if os.path.exists(outFile) != True:
            with requests.Session() as s:
                s.auth = (usr, pswrd)

                r1 = s.request('get', fileUrl)
                r2 = s.get(r1.url, auth=(usr, pswrd))

                with open(outFile, 'wb') as f:
                   f.write(r2.content)

    return outFile

def modis(date,h,v,outdir='./',creds=None,product=None):
    if outdir[-1] != '/':
        outdir = outdir+'/'

    acct = netrc.netrc(creds)
    usr,_,pswrd = acct.hosts['https://urs.earthdata.nasa.gov']

    today = datetime.datetime.now()

    if product:
        platform = product[:3]
        if platform.upper() == 'MOD':
            sensor = 'MOLT'
        else:
            sensor = 'MOLA'
    else:
        raise ValueError('product keyword was provided as None please specify product to fetch')

    basename = '{0}.A{1}{2:03d}.h{3:02d}v{4:02d}.006.hdf'

    if (today - date).days > 8:
        yr = date.year
        dt = (date-datetime.datetime(yr,1,1)).days + 1
        url = 'https://e4ftl01.cr.usgs.gov/{0}/{1}.006/{2}.{3:02d}.{4:02d}/'\
                .format(sensor,product,yr,date.month,date.day)

        with requests.Session() as s:
            s.auth = (usr, pswrd)

            r1 = s.request('get', url)
            r = s.get(r1.url, auth=(usr, pswrd))

            if r.ok:
                result = r.content.split(' '.encode())
                filtered = []
                for i in result:
                    if 'href'.encode() in i:
                        this = i.split('"'.encode())[1]
                        if this[-4:] =='.hdf'.encode():
                            filtered.append(this.decode("utf-8"))

                for f in filtered:
                    if 'h{:02d}'.format(h) in f\
                    and 'v{:02d}'.format(v)  in f:
                        filename = basename.format(product,yr,dt,h,v)
                        outFile = outdir + filename
                        if os.path.exists(outFile) != True:
                            newurl = url+f
                            r3 = s.request('get', newurl)
                            r4 = s.get(r3.url, auth=(usr, pswrd))

                            with open(outFile, 'wb') as this:
                                this.write(r4.content) # Say

    else:
        yr = date.year
        dt = (date-datetime.datetime(yr,1,1)).days + 1

        url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/6/{0}/Recent/'.format(product)
        filename = basename.format(product,yr,dt,h,v)
        fileUrl = url + filename[:-4] + '.NRT.hdf'

        outFile = outdir + filename
        if os.path.exists(outFile) != True:
            with requests.Session() as s:
                s.auth = (usr, pswrd)

                r1 = s.request('get', fileUrl)
                r2 = s.get(r1.url, auth=(usr, pswrd))

                with open(outFile, 'wb') as f:
                   f.write(r2.content)

    return outFile


def atms(date,region,outdir='./',creds=None):
    if outdir[-1] != '/':
        outdir = outdir+'/'

    acct = netrc.netrc(creds)
    usr,_,pswrd = acct.hosts['https://urs.earthdata.nasa.gov']

    foy = datetime.datetime(date.year,1,1)

    jday = (date-foy).days + 1

    url = 'https://sounder.gesdisc.eosdis.nasa.gov/data/SNPP_Sounder_Level1/SNPPATMSL1B.2/{0}/{1:03d}/'.format(
    date.year,jday
    )

    r = requests.get(url)

    if r.ok:
        result = r.content.split(' '.encode())
        filtered = []
        for i in result:
            if 'href'.encode() in i:
                this = i.split('"'.encode())[1]
                if this[-4:] =='.xml'.encode():
                    filtered.append(this.decode("utf-8"))

    xmls = set([url+xml for xml in filtered])
    sdrfiles = swathFilter(region,xmls)

    fileList = []

    with requests.Session() as s:
        s.auth = (usr, pswrd)
        for sdr in sdrfiles:
            outFile = os.path.join(outdir,sdr)
            if os.path.exists(outFile) != True:
                newurl = url+sdr
                r1 = s.request('get', newurl)
                r2 = s.get(r1.url, auth=(usr, pswrd))

                with open(outFile, 'wb') as this:
                    this.write(r2.content)

            fileList.append(outFile)

    return fileList

def swathFilter(region,xmls):
    geoms = []
    sdrnames = []
    for xml in xmls:
        xmlStr = requests.get(xml).content
        data = xmltodict.parse(xmlStr)

        ptList= data["S4PAGranuleMetaDataFile"]['SpatialDomainContainer']['HorizontalSpatialDomainContainer']['GPolygon']['Boundary']['Point']

        verts = [(float(pt['PointLongitude']),float(pt['PointLatitude'])) for pt in ptList]
        verts.append(verts[0])
        x,y = list(zip(*verts))

        maxDist = max([abs(x[-1]-x[i]) for i in range(len(x)-1)])

        if maxDist < 60:
            geoms.append(Polygon(verts))
            sdrnames.append(data["S4PAGranuleMetaDataFile"]['DataGranule']['GranuleID'])

    swathGeo = gpd.GeoDataFrame(pd.DataFrame({'sdr':sdrnames,'geometry':geoms}),geometry=geoms)

    swathGeo.crs = {'init':'epsg:4326'}

    intersection = gpd.overlay(region,swathGeo,how='intersection')

    return list(intersection.sdr)


def findTiles(region, tiles):
    """Returns the tile IDs that need to be downloaded for
    a given region bounded by *region*."""

    if region is None:
        raise ValueError("No bounding box provided for study area. Aborting download!")
        ids = None
    else:
        intersection = gpd.overlay(region,tiles,how='intersection')

        if 'PATH' in intersection.columns:
            h,v = 'PATH','ROW'
        elif 'h' in intersection.columns:
            h,v = 'h','v'
        else:
            raise AttributeError('cannot parse the needed tile information from provided geopadas dataframe')

        ids = [(intersection.iloc[i][h],intersection.iloc[i][v]) for i in range(len(intersection))]

    return ids
