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


    def waterMap(self,target_date,**kwargs):
        mapResult = geeutils.bootstrapOtsu(self.collection,target_date,**kwargs)

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


    def waterMap(self,hand,permanent=None,probablistic=False,nIters=100,elvStdDev=6,probTreshold=0.75):
        def _downscaleWrapper(iteration):
            i = ee.Number(iteration)
            handErr = hand.add(ee.Image.random(i).subtract(0.5).multiply(6)
                .multiply(elvStdDev))
            return downscale.bathtub(inImage,handErr,permanent)


        inImage = self.collection.mean().divide(10000)
        if probablistic:
            iters = ee.List.sequence(0,nIters-1)

            sims = ee.ImageCollection(iters.map(_downscaleWrapper))
            probs = sims.sum().divide(nIters).rename(['probability','error'])
            water = probs.select('probability').gt(probTreshold)
            mapResult = water.addBands(probs)

        else:
            mapResult = downscale.bathtub(inImage,hand,permanent)

        return mapResult



class Viirs(collectionDomain):
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


class Modis(collectionDomain):
    def __init__():
        super(Modis, self).__init__(region,time_start,time_end)
        return
