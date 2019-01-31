import os
import glob
import gzip
import netrc
import shutil
import ftplib
import tarfile
import requests
import datetime
import subprocess
import numpy as np
import pandas as pd
from osgeo import gdal
import geopandas as gpd
from shapely import geometry
from itertools import groupby
from landsat.google_download import GoogleDownload

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

def modis(date,h,v,outdir='./',creds=None,platform='terra',product=None):
    if outdir[-1] != '/':
        outdir = outdir+'/'

    acct = netrc.netrc(creds)
    usr,_,pswrd = acct.hosts['https://urs.earthdata.nasa.gov']

    today = datetime.datetime.now()

    if platform.lower() in ['terra','aqua']:
        if platform.lower() == 'terra':
            sensor = 'MOLT'
        else:
            sensor = 'MOLA'
    else:
        raise ValueError('platform options are "terra"|"aqua" entered values is : {}'.format(platform))

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


def landsat(date,p,r,outdir='./',updateScenes=False,maxClouds=100):
    if outdir[-1] != '/':
        outdir = outdir+'/'

    if updateScenes == True:
        # get the latest scene metadata
        subprocess.run(['landsat', '--update-scenes'],shell=True)

    sDate = (date - datetime.timedelta(1)).strftime('%Y-%m-%d')
    eDate = (date + datetime.timedelta(1)).strftime('%Y-%m-%d')

    downloader = GoogleDownload(sDate,eDate,8,path=p, row=r,
                                max_cloud_percent=maxClouds,
                                output_path=outdir)

    search = list(downloader.scenes_low_cloud.SCENE_ID)

    if len(search) == 1:
        fileList = search[0]
        downloader.download()
    else: fileList = None

    return None

def atms(date,outdir='./',creds=None):
    if outdir[-1] != '/':
        outdir = outdir+'/'

    today = datetime.datetime.now()

    tDiff = (today-date).days

    if tDiff > 95:
        raise ValueError('date argument is outside of range, must be within 95days of today')

    url = 'ftp-npp.bou.class.noaa.gov'

    ftp = ftplib.FTP(url)
    ftp.login('anonymous','anonymous')

    dateStr = date.strftime('%Y%m%d')

    products = ['ATMS-SDR-Ellipsoid-Geo','ATMS-SDR']

    basepath = '~/{0}/ATMS-SDR/{1}/NPP/'

    tarballs = []
    for i in products:
        ftp.cwd(basepath.format(dateStr,i))
        files = ftp.nlst()
        for file in files:
            if file[-4:] == '.tar':
                outFile = outdir + file
                tarballs.append(outFile)
                if os.path.exists(outFile) != True:
                    with open(outFile, 'wb') as f:
                       ftp.retrbinary('RETR ' + file, f.write)

    ftp.close()

    for t in tarballs:
        with tarfile.open(t,'r') as tar:
            tar.extractall(path=outdir)

    swathballs = glob.glob(os.path.join(outdir,'*.gz'))
    for sb in swathballs:
        outFile = sb[:-3]

        with gzip.open(sb,'rb') as f_in, open(outFile,"wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(sb)

    fileList = glob.glob(os.path.join(outdir,'*.h5'))

    return fileList

def groupSwathFiles(h5files):
    sortKey = lambda x: x.split('_')[-6]
    swathGroups, keys = [],[]
    for k,g in groupby(sorted(h5files,key=sortKey),sortKey):
        swathGroups.append(sorted(list(g)))
        keys.append(k)
    return swathGroups


def spatialSwathFilter(region,h5files,subgroup='//All_Data/ATMS-SDR-GEO_All/'):
    swathGroups = groupSwathFiles(h5files)

    geoms,sdrnames,geonames = [],[],[]
    for sg in swathGroups:
        geoFile,sdrFile = sg

        geoBase = 'HDF5:"{0}":{1}{2}'

        lats = gdal.Open(geoBase.format(geoFile,subgroup,'Latitude')).ReadAsArray()
        lons = gdal.Open(geoBase.format(geoFile,subgroup,'Longitude')).ReadAsArray()
        view = gdal.Open(geoBase.format(geoFile,subgroup,'SatelliteZenithAngle')).ReadAsArray()

        yp,xp = np.where(view<50)
        minindex,maxindex = xp.min(),xp.max()
        lons = lons[:,minindex:maxindex]
        lats = lats[:,minindex:maxindex]

        wVerts = [(lons[i,0],lats[i,0]) for i in range(lats.shape[0]) if (lons[i,0]>-200)and(lats[i,0]>-200)][::-1]
        nVerts = [(lons[0,i],lats[0,i]) for i in range(lats.shape[1]) if (lons[0,i]>-200)and(lats[0,i]>-200)]
        eVerts = [(lons[i,-1],lats[i,-1]) for i in range(lats.shape[0]) if (lons[i,-1]>-200)and(lats[i,-1]>-200)]
        sVerts = [(lons[-1,i],lats[-1,i]) for i in range(lats.shape[1]) if (lons[-1,i]>-200)and(lats[-1,i]>-200)]
        sVerts.append(wVerts[0])

        verts = wVerts + nVerts + eVerts + sVerts
        geoms.append(geometry.Polygon(verts))
        sdrnames.append(sdrFile)
        geonames.append(geoFile)

    swathGeo = gpd.GeoDataFrame(pd.DataFrame({'sdr':sdrnames,'geo':geonames,'geometry':geoms}),geometry=geoms)
    swathGeo.crs = {'init':'epsg:4326'}

    intersection = gpd.overlay(region,swathGeo,how='intersection')

    return list(intersection.sdr),list(intersection.geo)

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
