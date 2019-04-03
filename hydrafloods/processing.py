from __future__ import absolute_import
import os
import ee
from ee.ee_exception import EEException
from . import geeutils, downscale, fetch, preprocess, utils


class collectionDomain(object):
    def __init__(self,region,time_start,time_end):
        # TODO: add exceptions to check datatypes
        self.region = region # dtype = ee.Geometry
        self.iniTime = time_start
        self.endTime = time_end
        # self.iniEE = ee.Date(time_start) # dtype = str as YYYY-MM-DD
        # self.endEE = ee.Date(time_end) # dtype = ee.Date
        return



class Sentinel1(collectionDomain):
    def __init__(self,region,time_start,time_end):
        super(Sentinel1, self).__init__(region,time_start,time_end)

        self.collection = ee.ImageCollection('COPERNICUS/S1_GRD')\
                            .filterBounds(self.region)\
                            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
                            .filterDate(self.iniTime, self.endTime)\
                            .map(self._maskEdges)\
                            .select('VV')

        return


    def waterMap(self,target_date,
                canny_threshold=7,    # threshold for canny edge detection
                canny_sigma=1,        # sigma value for gaussian filter
                canny_lt=7,           # lower threshold for canny detection
                smoothing=100,        # amount of smoothing in meters
                connected_pixels=200, # maximum size of the neighborhood in pixels
                edge_length=50,       # minimum length of edges from canny detection
                smooth_edges=100,
                ):

        mapResult = geeutils.bootstrapOtsu(self.collection,target_date,canny_threshold,
                                           canny_sigma,canny_lt,smoothing,
                                           connected_pixels,edge_length,
                                           smooth_edges)

        return mapResult


    def _maskEdges(self,img):
        angles = img.select('angle')
        return img.updateMask(angles.lt(45).And(angles.gt(31)))



class Atms(collectionDomain):
    def __init__(self,region,time_start,time_end,collectionid=''):
        super(Atms, self).__init__(region,time_start,time_end)

        self.collection = ee.ImageCollection(collectionid)\
                            .filterBounds(self.region)\
                            .filterDate(self.iniTime,self.endTime)
        return


    def extract(self,date,region,outdir='./',creds=None,gridding_radius=50000):
        files = fetch.atms(date,region,outdir,creds)
        geotiffs = list(map(lambda x: preprocess.atms(x,gridding_radius), files))
        return geotiffs


    def load(self,files,gcsBucket='',eeAsset=''):
        if gcsBucket[-1] != '/':
            gcsBucket += '/'
        if eeAsset[-1] != '/':
            eeAsset += '/'

        for f in files:
            fName = os.path.basename(f)
            utils.push_to_gcs(f,gcsBucket)
            bso = gcsBucket + fName
            tstr = fName.split('.')[3]
            t = '{0}-{1}-{2}:{3}:00'.format(tstr[:4],tstr[4:6],tstr[6:11],tstr[11:])
            utils.push_to_gee(bso,eeAsset,properties={'time_start':t})

        return


    def waterMap(self,hand,permanent=None):
        inImage = self.collection.mean().divide(10000)
        mapResult = downscale.bathtub(inImage,hand,permanent)

        return mapResult



class Viirs(object):
    def __init__(self,region,time_start,time_end,collectionid=''):
        super(Viirs, self).__init__(region,time_start,time_end)

        self.collection = ee.ImageCollection(collectionid)\
                            .filterBounds(self.region)\
                            .filterDate(self.iniTime,self.endTime)
        return


    def extract(self,date,region,outdir='./',creds=None):

        return

    def load(self,files,gcsBucket='',eeAsset=''):

        return

    def downscale():

        return

    def waterMap():

        return


class Modis(object):
    def __init__():
        return
