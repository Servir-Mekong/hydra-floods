from __future__ import absolute_import
import os
import ee
from ee.ee_exception import EEException
from . import geeutils, downscale, fetch, preprocess, utils


INITIME = ee.Date.fromYMD(1970,1,1)
BANDREMAP = ee.Dictionary({
  'landsat': ee.List(['B2','B3','B4','B5','B6','B7','time']),
  'viirs': ee.List(['M2','M4','M5','M7','M10','M11','time']),
  'sen2':  ee.List(['B2','B3','B4','B8','B11','B12','time']),
  'modis': ee.List(['sur_refl_b03','sur_refl_b04','sur_refl_b01','sur_refl_b02','sur_refl_b06','sur_refl_b07']),
  'new': ee.List(['blue','green','red','nir','swir1','swir2','time'])
})

class hfCollection(object):
    def __init__(self,region,time_start,time_end,collectionid=''):
        # TODO: add exceptions to check datatypes
        self.region = region # dtype = ee.Geometry
        self.iniTime = time_start
        self.endTime = time_end
        self.id = collectionid

        self.collection = ee.ImageCollection(self.id)\
                            .filterBounds(self.region)\
                            .filterDate(self.iniTime,self.endTime)

        return



class Sentinel1(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Sentinel1, self).__init__(*args,**kwargs)

        self.collection = self.collection\
                            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
                            .map(self._maskEdges)\
                            .select('VV')

        return


    def waterMap(self,target_date,**kwargs):
        mapResult = geeutils.bootstrapOtsu(self.collection,target_date,**kwargs)

        return mapResult


    def _maskEdges(self,img):
        angles = img.select('angle')
        return img.updateMask(angles.lt(45).And(angles.gt(31)))



class Atms(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Atms, self).__init__(*args,**kwargs)

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
            water = probs.select(['probability'],['water']).gt(probTreshold)
            mapResult = water.addBands(probs.multiply(10000).uint16())

        else:
            mapResult = downscale.bathtub(inImage,hand,permanent)

        return mapResult



class Viirs(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Viirs, self).__init__(*args,**kwargs)

        self.collection = self.collection.map(self._qaMask)\
            .select(BANDREMAP.get('viirs'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def _qaMask(self,img):
        viewing = img.select('SensorZenith').abs().multiply(0.01).lt(30)
        clouds = geeutils.extractBits(img.select('QF1'),2,3,'cloud_qa').lte(2)
        shadows = geeutils.extractBits(img.select('QF2'),3,3,'shadow_qa').eq(0)
        aerosols = geeutils.extractBits(img.select('QF2'),4,4,'aerosol_qa').eq(0)
        mask = clouds.And(shadows).And(aerosols)
        t = ee.Date(img.get('system:time_start'))
        nDays = t.difference(INITIME,'day')
        time = ee.Image(nDays).int16().rename('time')
        return img.updateMask(mask).addBands(time)


    def extract(self,date,region,outdir='./',creds=None):

        return

    def load(self,files,gcsBucket='',eeAsset=''):

        return

    def downscale(self,fineCollection,**kwargs):
        result = downscale.starfm(fineCollection,self.collection,**kwargs)
        self.downscaled = result
        return

    def waterMap(self,target_date,hand,probablistic=False,nIters=100,probTreshold=0.75,**kwargs):
        def _threholdWrapper(iteration):
            i = ee.Number(iteration)
            sim = geeutils.globalOtsu(self.downscaled,target_date,self.region,seed=i,**kwargs)
            return sim.updateMask(hand.lt(30))

        tDate = ee.Date(target_date)

        if probablistic:
            iters = ee.List.sequence(0,nIters-1)

            sims = ee.ImageCollection(iters.map(_threholdWrapper))
            probs = sims.sum().divide(nIters).rename(['probability'])
            water = probs.select(['probability']).gt(probTreshold).rename('water')
            mapResult = water.addBands(probs.multiply(10000).uint16())

        else:
            mapResult = geeutils.globalOtsu(self.downscaled,target_date,self.region,**kwargs)\
                .updateMask(hand.lt(30))


        return mapResult



class Modis(hfCollection):
    def __init__(self,*args,**kwargss):
        super(Modis, self).__init__(*args,**kwargs)
        return

    def _qaMask(img):

        return



class Landsat(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Landsat, self).__init__(*args,**kwargs)

        self.collection = self.collection.map(self._qaMask)\
            .select(BANDREMAP.get('landsat'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def _qaMask(self,img):
        qaCloud = geeutils.extractBits(img.select('pixel_qa'),5,5,'qa').neq(1)
        qaShadow = geeutils.extractBits(img.select('pixel_qa'),3,3,'qa').neq(1)
        mask = qaCloud.And(qaShadow)
        t = ee.Date(img.get('system:time_start'))
        nDays = t.difference(INITIME,'day')
        time = ee.Image(nDays).int16().rename('time')
        return img.updateMask(mask).addBands(time)




class Sentinel2(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Sentinel2, self).__init__(*args,**kwargs)

        return

    def _qaMask():
        return
