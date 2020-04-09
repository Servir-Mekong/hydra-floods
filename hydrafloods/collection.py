from __future__ import absolute_import
import os
import ee
import math
from pprint import pformat
from ee.ee_exception import EEException
from hydrafloods import geeutils, thresholding, downscale, fetch, preprocess, utils,filtering


INITIME = ee.Date('1971-01-01T00:00:00')
BANDREMAP = ee.Dictionary({
  'landsat': ee.List(['B2','B3','B4','B5','B6','B7','time']),
  'viirs': ee.List(['M2','M4','M5','M7','M10','M11','time']),
  'sen2':  ee.List(['B2','B3','B4','B8','B11','B12','time']),
  'modis': ee.List(['sur_refl_b03','sur_refl_b04','sur_refl_b01','sur_refl_b02','sur_refl_b06','sur_refl_b07','time']),
  'new': ee.List(['blue','green','red','nir','swir1','swir2','time'])
})

class hfCollection(object):
    def __init__(self,region,time_start,time_end,collectionid='',useQa=True):
        # TODO: add exceptions to check datatypes
        self.region = region # dtype = ee.Geometry
        self.iniTime = time_start
        self.endTime = time_end
        self.id = collectionid
        self.useQa = useQa

        self.collection = ee.ImageCollection(self.id)\
                            .filterBounds(self.region)\
                            .filterDate(self.iniTime,self.endTime)

        return

    def __repr__(self):
        objDict = {
            'name':self.__name__,
            'collectionid':self.id,
            'startTime':self.iniTime.strftime('%Y-%m-%d'),
            'endTime':self.endTime.strftime('%Y-%m-%d'),
            'region':self.region.coordinates().getInfo()
        }
        strRepr = pformat(objDict, depth=3)
        return f'HYDRAFloods Collection:\n{strRepr}'

    @property
    def __name__(self):
        return self.__class__.__name__

    def clip(self,img):
        return img.clip(self.region)



class Sentinel1(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Sentinel1, self).__init__(*args,**kwargs)

        if self.useQa:
            self.collection = self.collection.map(self._qa)

        self.collection = self.collection\
                            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))

        return

    def _qa(self,img):
        angles = img.select('angle')
        return img.updateMask(angles.lt(45).And(angles.gt(31)))


    def waterMap(self,target_date,**kwargs):
        mapResult = thresholding.bmaxOtsu(self.collection,target_date,self.region,**kwargs)

        return mapResult



class Atms(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Atms, self).__init__(*args,**kwargs)

        return

    def _qa(self,img):
        return


    def extract(self,date,region,credentials,outDir='./',gridding_radius=50000):
        files = fetch.atms(credentials,startTime=date,endTime=None,region=region,outDir=outDir)
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

        if self.useQa:
            self.collection = .map(self._qa)

        self.collection = self.collection\
            .select(BANDREMAP.get('viirs'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def _qa(self,img):
        cloudBit = int(math.pow(2,2))
        shadowBit = int(math.pow(2,3))
        snowBit = int(math.pow(2,5))

        viewing = img.select('SensorZenith').abs().multiply(0.01).lt(55)
        clouds = img.select('QF1').bitwiseAnd(shadowBit).eq(0)
        shadows = img.select('QF2').bitwiseAnd(shadowBit).eq(0)
        snows = img.select('QF2').bitwiseAnd(snowBit).eq(0)

        mask = clouds.And(shadows).And(snows).And(viewing)
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
            return sim.And(hand.lt(30))

        tDate = ee.Date(target_date)

        if probablistic:
            iters = ee.List.sequence(0,nIters-1)

            sims = ee.ImageCollection(iters.map(_threholdWrapper))
            probs = sims.sum().divide(nIters).rename(['probability'])
            water = probs.select(['probability']).gt(probTreshold).rename('water')
            mapResult = water.addBands(probs.multiply(10000).uint16())

        else:
            mapResult = geeutils.globalOtsu(self.downscaled,target_date,self.region,**kwargs)\
                .And(hand.lt(30))

        return mapResult



class Modis(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Modis, self).__init__(*args,**kwargs)

        if self.useQa:
            self.collection = self.collection.map(self._qa)

        self.collection = self.collection\
            .select(BANDREMAP.get('modis'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def _qa(self,img):
        cloudBit = int(math.pow(2,10))
        shadowBit = int(math.pow(2,2))
        snowBit = int(math.pow(2,15))
        viewing = img.select('SensorZenith').abs().multiply(0.01).lt(45)
        clouds = img.select('state_1km').bitwiseAnd(cloudBit).eq(0)
        shadows = img.select('state_1km').bitwiseAnd(shadowBit).eq(0)
        snows = img.select('state_1km').bitwiseAnd(snowBit).eq(0)
        mask = clouds.And(shadows).And(snows)#.And(viewing)
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
            return sim.And(hand.lt(30))

        tDate = ee.Date(target_date)

        if probablistic:
            iters = ee.List.sequence(0,nIters-1)

            sims = ee.ImageCollection(iters.map(_threholdWrapper))
            probs = sims.sum().divide(nIters).rename(['probability'])
            water = probs.select(['probability']).gt(probTreshold).rename('water')
            mapResult = water.addBands(probs.multiply(10000).uint16())

        else:
            mapResult = geeutils.globalOtsu(self.downscaled,target_date,self.region,**kwargs)\
                .And(hand.lt(30))

        return mapResult



class Landsat(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Landsat, self).__init__(*args,**kwargs)

        if self.useQa:
            self.collection = self.collection.map(self._qa)

        self.collection = self.collection\
            .select(BANDREMAP.get('landsat'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def waterMap(self,target_date,**kwargs):
        mapResult = thresholding.bmaxOtsu(self.collection,target_date,self.region,**kwargs)

        return mapResult


    def _qa(self,img):
        cloudBit = int(math.pow(2,5))
        shadowBit = int(math.pow(2,3))
        snowBit = int(math.pow(2,4))

        qaCloud = img.select('pixel_qa').bitwiseAnd(cloudBit).eq(0)
        qaShadow = img.select('pixel_qa').bitwiseAnd(shadowBit).eq(0)
        qaSnow = img.select('pixel_qa').bitwiseAnd(snowBit).eq(0)
        mask = qaCloud.And(qaShadow).And(qaSnow)
        t = ee.Date(img.get('system:time_start'))
        nDays = t.difference(INITIME,'day')
        time = ee.Image(nDays).int16().rename('time')
        return img.updateMask(mask).addBands(time).uint16()




class Sentinel2(hfCollection):
    def __init__(self,*args,**kwargs):
        super(Sentinel2, self).__init__(*args,**kwargs)

        if self.useQa:
            self.collection = self.collection.map(self._qa)

        self.collection = self.collection\
            .select(BANDREMAP.get('sen2'),BANDREMAP.get('new'))\
            .map(geeutils.addIndices)

        return

    def _qa(self,img):
        sclImg = img.select('SCL') # Scene Classification Map
        mask = sclImg.gte(4).And(sclImg.lte(6))
        t = ee.Date(img.get('system:time_start'))
        nDays = t.difference(INITIME,'day')
        time = ee.Image(nDays).int16().rename('time')
        return img.updateMask(mask).addBands(time).uint16()

    def _bandPassAdjustment(self,img):
        bands = ee.List(BANDREMAP.get('new'))
        # linear regression coefficients for adjustment
        gain = ee.Array([[0.9778], [1.0053], [0.9765], [0.9983], [0.9987], [1.003],[1.0]])
        bias = ee.Array([[-0.00411],[-0.00093],[0.00094],[-0.0001],[-0.0015],[-0.0012],[0.0]])
        # Make an Array Image, with a 1-D Array per pixel.
        arrayImage1D = img.select(bands).toArray()

        # Make an Array Image with a 2-D Array per pixel, 6x1.
        arrayImage2D = arrayImage1D.toArray(1)

        componentsImage = ee.Image(gain).multiply(arrayImage2D).add(ee.Image(bias))\
        .arrayProject([0])\
        .arrayFlatten([bands]).float()

        return componentsImage.uint16().set('system:time_start',img.get('system:time_start'))
