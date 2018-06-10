from __future__ import print_function,division
import copy
import glob
import math
import datetime
import itertools
import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib.pyplot as plt

import geopandas as gpd
from affine import Affine
from osgeo import gdal,osr
from gdalconst import GA_ReadOnly
from pyproj import Proj, transform
from shapely.geometry import Point, Polygon

import cartopy
import cartopy.crs as ccrs

from . import utils


class raster(object):
    
    def __init__(self,path,sensor):
        self.src = path
        self.sensor = sensor
        self.crs = {'init':'epsg:6933'}
    
        return
    
    
    def _copy(self):
        return copy.deepcopy(self)
    
#    def __str__(self):
#        print()
    
    
    def _makeShape(self,linearRing,shapeid=None):
        if shapeid == None:
            shapeid = 'box'
            
        inProj = Proj(init='epsg:4326')
        outProj = Proj(init='epsg:6933')
                
        geom = []
        linRing = []
        for n in linearRing:
            x,y = transform(inProj,outProj,n[1],n[0])
            linRing.append((x,y))
            vertex = [shapeid,Point((x,y))]
            geom.append(vertex)
            
        self.linearRing = linRing
                                
        df = gpd.GeoDataFrame(geom, columns = ['shape_id', 'geometry'], 
                              geometry='geometry')
        
        df['geometry'] = df['geometry'].apply(lambda x: x.coords[0])

        df = df.groupby('shape_id')['geometry'].apply(lambda x: Polygon(x.tolist())).reset_index()

        # Declare the result as a new a GeoDataFrame
        df = gpd.GeoDataFrame(df, geometry = 'geometry')
        df.crs = {'init':'epsg:6933'}
                        
        return df
    

    def _extractBits(self,image,start,end):
        """Helper function to convert Quality Assurance band bit information to flag values

        Args:
            image (ndarray): Quality assurance image as a numpy array
            start (int): Bit position to start value conversion
            end (int): Bit position to end value conversion

        Returns:
            out (ndarray): Output quality assurance in values from bit range
        """  

        pattern = 0;
        for i in range(start,end+1):
            pattern += math.pow(2, i)

        bits = image.astype(np.uint16) & int(pattern)
        out = bits >> start

        return out

    
    def _geoCoords(self,north,south,east,west,dims,projStr=None):
        
        outProj = Proj(init='epsg:6933')
        gcsProj = Proj(init='epsg:4326')
        
        if self.sensor == 'viirs':
            projStr = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
            native = Proj(projStr)

            llx,lly = transform(gcsProj,native,west,south)
            urx,ury = transform(gcsProj,native,east,north)

            yCoords = np.linspace(lly,ury,dims[0],endpoint=False)[::-1]
            xCoords = np.linspace(llx,urx,dims[1],endpoint=False)

            xx,yy = np.meshgrid(xCoords,yCoords)

        elif self.sensor == 'landsat':        
            native = Proj(projStr)

            llx,lly = west,south
            urx,ury = east,north

            yCoords = np.linspace(lly,ury,dims[0],endpoint=False)[::-1]
            xCoords = np.linspace(llx,urx,dims[1],endpoint=False)

            xx,yy = np.meshgrid(xCoords,yCoords)
            
        else:
            raise NotImplementedError("Raster from sensor {} is not yet implemented".format(self.sensor))

        lons,lats = transform(native,outProj,xx,yy)

        return lons,lats

    
    def normalizedDifference(self,band1,band2):

        nd = (self.bands[band1] - self.bands[band2]) / (self.bands[band1] + self.bands[band2])
        
        self.bands['nd'] = nd
        
        clone = self._copy()

        for n in clone.bandNames:
            del clone.bands[n]
            
        clone.bandNames = ['nd']
        
        return clone


class viirs(raster):
    
    def __init__(self,infile):
        raster.__init__(self,infile,'viirs')

        tree = '//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/'
        field = 'SurfReflect_{0}{1}_1'
        base = 'HDF5:"{0}":{1}{2}'

        m = [i for i in range(12) if i not in [0,6,9]]
        i = [i for i in range(1,4)]
        bands = [m,i]

        res = ['1km','500m']
        mode = ['M','I']

        band = gdal.Open(base.format(self.src,tree.format('1km'),field.format('QF',1)))
        self.metadata = band.GetMetadata()
        cloudQA = self._extractBits(band.ReadAsArray(),2,3)
        hiresCloudQA = ndimage.zoom(cloudQA,2,order=0)
        band = None

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',2)))
        shadowQA = self._extractBits(band.ReadAsArray(),3,3)
        hiresShadowQA = ndimage.zoom(shadowQA,2,order=0)

        qa = (cloudQA>0)&(shadowQA<1)
        hiresqa = (hiresCloudQA>0)&(hiresShadowQA<1)
        
        yList = self.metadata['GRingLatitude'].strip().split(' ')
        xList = self.metadata['GRingLongitude'].strip().split(' ')
        
        xx,yy = [],[]
        for i in range(len(yList)):
            xx.append(float(xList[i]))
            yy.append(float(yList[i]))
            
        linRing =  list(zip(*[yy,xx]))
                        
        self.shape = self._makeShape(linRing,shapeid=self.sensor)

        east,west = float(self.metadata['EastBoundingCoord']), float(self.metadata['WestBoundingCoord'])
        north,south = float(self.metadata['NorthBoundingCoord']), float(self.metadata['SouthBoundingCoord'])
        
        self.extent = [west,south,east,north]

        databands = {'QA':qa}
        
        bandNames = ['QA']

        for i in range(2):
            if i == 0:
                mask = qa
            else:
                mask = hiresqa
            for j in range(len(bands[i])):

                subdataset = base.format(infile,tree.format(res[i]),field.format(mode[i],bands[i][j]))

                band = gdal.Open(subdataset)
                data = band.ReadAsArray()
                data = np.ma.masked_where(data<0,data)
                data = np.ma.masked_where(data>10000,data)
                data = np.ma.masked_where(mask!=0,data)
                bName = '{0}{1}'.format(mode[i],bands[i][j])
                databands[bName] = data.astype(np.int16)
                bandNames.append(bName)

                band = None
                data = None
                
        self.bands = databands
        self.bandNames = bandNames
        
        coords = {}

        coords['MLon'],coords['MLat'] = self._geoCoords(north,south,east,west,self.bands['M1'].shape)
        coords['ILon'],coords['ILat'] = self._geoCoords(north,south,east,west,self.bands['I1'].shape)
        
        self.coords = coords
        
        self.crs = {'init':'epsg:6974'}
        self.proj = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'

        mgt, igt = self._getViirsGt(north,west)
        
        transform = {}
        transform['Mtransform'] = Affine.from_gdal(*mgt)
        transform['Itransform'] = Affine.from_gdal(*igt)
        self.transform = transform
        
        date = '{0}{1}{2}'.format(self.metadata['RangeBeginningDate'],self.metadata['RangeBeginningTime'],' UTC')
        
        self.date = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')

        return

    def _getViirsGt(self,north,west):
        projStr = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
        outProj = Proj(projStr)
        inProj = Proj(init='epsg:6933')

        ulx,uly = transform(inProj,outProj,west,north)

        #(originX, pixelWidth, 0, originY, 0, pixelHeight)

        igt = (ulx,500,0,uly,0,-500)
        mgt = (ulx,1000,0,uly,0,-1000)

        return mgt,igt

    
    def plot(self,bandKey,latKey,lonKey,vmin=0,vmax=1,cmap='gray'):
        inProj = Proj(init='epsg:6933')
        outProj = Proj(init='epsg:4326')
        xx,yy = transform(inProj,outProj,self.coords[lonKey],self.coords[latKey])
        
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.pcolormesh(xx, yy,self.bands[bandKey]\
                       ,transform=ccrs.PlateCarree(),
                       vmin=vmin,vmax=vmax,cmap=cmap)

        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        
        return

    

class atms(raster):

    def __init__(infile, ingeo):
        raster.__init__(self,infile,'atms')
    
        trees = ['//All_Data/ATMS-SDR_All/',
                 '//All_Data/ATMS-SDR-GEO_All/'
                ]

        fields = [['GainCalibration','BrightnessTemperature',],
                  ['Latitude','Longitude',]
                 ]

        base = 'HDF5:"{0}":{1}{2}'    

        band = gdal.Open(base.format(infile,trees[0],fields[0][3]))
        metadata = band.GetMetadata()
        qa = _extractBits(band.ReadAsArray(),2,3)
        hiresqa = ndimage.zoom(qa,2,order=0)
        band = None

        out = {'QA':qa}
        
        btrBands = ['C1','C2','C3','C4','C5','C16','C17','C18','C19','C20','C21','C22']
        toKeep = [0,1,2,3,4,15,16,17,18,19,20,21]

        for i in range(len(trees)):
            if i == 0:
                f = infile
                iters = len(fields[i])-1
            else:
                f = ingeo
                iters = len(fields[i])

            for j in range(iters):
                subdataset = base.format(infile,tree.format(res[i]),field.format(mode[i],bands[i][j]))

                band = gdal.Open(subdataset)
                data = band.ReadAsArray()
                if i == 0:
                    if j == 1:
                        for b in range(toKeep):
                           #out[btrBands[b]] = 
                            data = np.ma.masked_where(data<0,data)
                            data = np.ma.masked_where(mask!=0,data)
                out['{0}{1}'.format(mode[i],bands[i][j])] = data/10000.

                band = None
                data = None


        return out
    
    
    def _project():
        
        return
    
    def _calibrate():
        
        
        return calibrated
    
    
    def plot():
        
        return


    
class landsat(raster):

    def __init__(self,infile):
        raster.__init__(self,infile,'landsat')
        
        self.metadata = self._parseMetadata()
        
        path = '/'.join(infile.split('/')[0:-1])+'/'

        sZenith = 90 - float(self.metadata['SUN_ELEVATION'])
        dSun = float(self.metadata['EARTH_SUN_DISTANCE'])

        north = max(float(self.metadata['CORNER_UR_PROJECTION_Y_PRODUCT']),
                    float(self.metadata['CORNER_UL_PROJECTION_Y_PRODUCT']))
        south = min(float(self.metadata['CORNER_LR_PROJECTION_Y_PRODUCT']),
                    float(self.metadata['CORNER_LL_PROJECTION_Y_PRODUCT']))
        east  = max(float(self.metadata['CORNER_UR_PROJECTION_X_PRODUCT']),
                    float(self.metadata['CORNER_LR_PROJECTION_X_PRODUCT']))
        west  = min(float(self.metadata['CORNER_UL_PROJECTION_X_PRODUCT']),
                    float(self.metadata['CORNER_LL_PROJECTION_X_PRODUCT']))
        
        xx,yy = [],[]
        for c in ['LL','UL','UR','LR']:
            xx.append(float(self.metadata['CORNER_{}_LON_PRODUCT'.format(c)]))
            yy.append(float(self.metadata['CORNER_{}_LAT_PRODUCT'.format(c)]))
            
        linRing =  list(zip(*[yy,xx]))
            
        self.shape = self._makeShape(linRing,shapeid=self.sensor)

        projStr = '+proj=utm +zone={} +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs'\
                    .format(self.metadata['UTM_ZONE'])

        bandList = ['QUALITY']
        for i in range(2,8):
            bandList.append(str(i))

        bandKey = 'FILE_NAME_BAND_{}'
        gainKey = 'REFLECTANCE_MULT_BAND_{}'
        biasKey = 'REFLECTANCE_ADD_BAND_{}'

        outKeys = ['BQA','B2','B3','B4','B5','B6','B7']

        bands = {}

        for i in range(len(bandList)):
            name = self.metadata[bandKey.format(bandList[i])]
            ds = gdal.Open('{0}{1}'.format(path,name),GA_ReadOnly)

            tmp = ds.ReadAsArray()
            tmp = np.ma.masked_where(tmp == 0,tmp)

            if i == 0:
                gt  = ds.GetGeoTransform()
                dims = tmp.shape

                bands[outKeys[i]] = tmp

                qa = self._extractBits(tmp,4,4)

            else:
                gain = float(self.metadata[gainKey.format(bandList[i])])
                bias = float(self.metadata[biasKey.format(bandList[i])])

                tmp = np.ma.masked_where(qa == 1, tmp)

                toa = (((gain * tmp) + bias) / np.cos(np.deg2rad(sZenith))) * 10000

                bands[outKeys[i]] = toa.astype(np.uint16)

            ds = None
            band = None 
            
        self.bands = bands
        self.bandNames = outKeys
        
        coords = {}
        coords['Lon'],coords['Lat'] = self._geoCoords(north,south,east,west,dims,projStr=projStr)
        self.coords = coords
        
        inProj = Proj(projStr)
        outProj = Proj(init='epsg:6933')
        
        newgt=list(gt)
        newgt[0],newgt[3]= transform(inProj,outProj,gt[0],gt[3])
        
        self.transform = Affine.from_gdal(*newgt)
        
        date = '{0} {1}{2}'.format(self.metadata['DATE_ACQUIRED'],self.metadata['SCENE_CENTER_TIME'][:-3], ' UTC')
        self.date = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')
        
        return
        
        
    def _parseMetadata(self):
        with open(self.src,'r') as f:
            data = f.read()

        split_metadata = data.split('\n')

        output = {}
        for x in split_metadata:
            if "=" in x:
                line = x.split("=")
                output[line[0].strip()] = line[1].strip()
                clean_output = {key: item.strip('"') for key, item in output.items()}

        return clean_output

    def plot(self,bandKey,latKey,lonKey,vmin=0,vmax=1,cmap='gray'):
        inProj = Proj(init='epsg:6933')
        outProj = Proj(init='epsg:4326')
        xx,yy = transform(inProj,outProj,self.coords[lonKey],self.coords[latKey])
        
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.pcolormesh(np.fliplr(xx), yy,np.fliplr(self.bands[bandKey])\
                       ,transform=ccrs.PlateCarree(),
                       vmin=vmin,vmax=vmax,cmap=cmap)

        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        
        return

            
    