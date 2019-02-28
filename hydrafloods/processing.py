import ee
import numpy as np
from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank

import rastersmith as rs

from . import geeutils, downscale


class eeProcess(object):
    def __init__(time_start,time_end,region):
        self.iniTime = time_start
        self.endTime = time_end
        self.region = region

    def export(self,image):



class Sentinel1(eeProcess):
    def __init__(self,time_start,time_end,region):
        eeProcess.__init__(self,time_start,time_end,region)
        return

    @staticmethod
    def exportFloodMap()

    def getFloodMap(gr,time_start,time_end,
                    canny_threshold=7,    # threshold for canny edge detection
                    canny_sigma=1,        # sigma value for gaussian filter
                    canny_lt=7,           # lower threshold for canny detection
                    smoothing=100,        # amount of smoothing in meters
                    connected_pixels=200, # maximum size of the neighborhood in pixels
                    edge_length=50,       # minimum length of edges from canny detection
                    smooth_edges=100,
                    ):

        geom = ee.Geometry.Rectangle([gr.west,gr.south,gr.east,gr.north])

        mapResult = geeutils.s1WaterMap(geom,time_start,time_end,canny_threshold,
                                        canny_sigma,canny_lt,smoothing,
                                        connected_pixels,edge_length,
                                        smooth_edges)

        return mapResult


class Atms(object):
    def __init__(self):
        return

    @staticmethod
    def maskClouds(ds,threshold=-20):
        rain = ds.sel(band='C16').astype(np.float) - ds.sel(band='C1').astype(np.float)

        cloudMask = rain > threshold

        return ds.raster.updateMask(cloudMask)

    @classmethod
    def getWaterFraction(cls,ds,cloudThresh=-20,constrain=True,maskClouds=True):

        if maskClouds:
            atmsNoClouds = cls.maskClouds(ds,threshold=cloudThresh)
        else:
            atmsNoClouds = ds.copy()

        dBtr = atmsNoClouds.sel(band='C4').astype(np.float) - atmsNoClouds.sel(band='C3').astype(np.float)
        dBtr.coords['band'] = 'dBtr'

        channels = xr.concat([atmsNoClouds.sel(band=['C3','C4','C16']).isel(time=0,z=0),
                              dBtr.isel(time=0,z=0)],dim='band')
        arr = channels.values
        arr[np.isnan(arr)] = -9999

        nClasses = 3
        nfindr = eea.NFINDR()
        U = nfindr.extract(arr, nClasses, maxit=100, normalize=True, ATGP_init=True)

        drop = np.argmin(list(map(lambda x:U[x,:].mean(),range(nClasses))))
        waterIdx = np.argmin(list(map(lambda x:np.delete(U,drop,axis=1)[x,:],range(nClasses-2))))

        if waterIdx == 0:
            bandList = ['water','land','mask']
        else:
            bandList = ['land','water','mask']

        nnls = amp.NNLS()
        amaps = nnls.map(arr, U, normalize=True)

        drop = np.argmin(list(map(lambda x:amaps[:,:,x].mean(),range(amaps.shape[2]))))

        unmixed = np.delete(amaps,drop,axis=2)

        unmixed[unmixed==0] = np.nan

        scaled = np.zeros_like(unmixed)
        for i in range(scaled.shape[2]):
            summed = unmixed[:,:,i]/unmixed.sum(axis=2)
            scaled[:,:,i] = (summed - np.nanmin(summed)) / (np.nanmax(summed) - np.nanmin(summed))

        scaled = scaled - 0.25
        scaled[scaled<0] = 0

        fWater = atmsNoClouds.sel(band=['C1','C2','mask']).copy()
        fWater[:,:,0,:2,0] = scaled[:,:,:]
        fWater.coords['band'] = bandList

        return fWater.raster.updateMask(atmsNoClouds.sel(band='mask'))



class Landsat8(object):
    def __init__():
        return


class Viirs(object):
    def __init__(self):
        return

    @staticmethod
    def getWaterMask(ds,transform=True):
        if transform:
            ds = 1 / (1 + (np.e ** ds))

        arr = ds.isel(time=0,z=0,band=0).values
        global_otsu = threshold_otsu(arr[~np.isnan(arr)])
        waterMask = ds >= global_otsu

        return waterMask.where(waterMask>0)


class Modis(object):
    def __init__():
        return


class Sentinel2(object):
    def __init__():
        return
