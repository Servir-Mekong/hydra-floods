from __future__ import absolute_import
import os
import ee
import fire
import glob
import yaml
import datetime
import warnings
import pandas as pd
import geopandas as gpd

from hydrafloods import utils
from hydrafloods.processing import *

class hydrafloods(object):
    def __init__(self,config=None):
        self.configuration = config
        return

    def _parse_config(self):
        if self.configuration:
            self.filePath = os.path.dirname(os.path.abspath(__file__))
            yamlFile = self.configuration

            with open(self.configuration,'r') as stream:
                try:
                    struct = yaml.load(stream,Loader=yaml.SafeLoader)
                except yaml.YAMLError as exc:
                    print(exc)

            conf = struct['configuration']
            self.conf = conf

            prcs = struct['process']
            self.prcs = prcs

            confKeys = list(conf.keys())
            prcsKeys = list(prcs.keys())


            parse top-level configuration key information
            if 'name' in confKeys:
                self.name = conf['name']
            else:
                raise AttributeError('provided yaml file does not have a name parameter in configuration')

            if 'region' in confKeys:
                shp = gpd.read_file(conf['region'])
                self.region = list(shp.bounds.values[0])

            elif 'country' in confKeys:
                world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
                country = world[world.name == conf['country']]
                if len(country) >= 1:
                    self.region = list(country.bounds.values[0])
                else:
                    raise ValueError('could not parse selected country from world shapefile')

            elif 'boundingbox' in confKeys:
                # from shapely import geometry
                # self.region = gpd.GeoDataFrame(pd.DataFrame({'id':[0],'geometry':[geometry.box(*conf['boundingbox'])]}))
                self.region = conf['boundingbox']

            else:
                raise AttributeError('provided yaml file does not have a specified region in configuration')

            if 'credentials' in confKeys:
                self.credentials = conf['credentials']
                self.earthdataLogin = self.credentials['earthdata'].values()
            else:
                self.credentials = None

            if 'stagingBucket' in confKeys:
                self.stagingBucket = conf['stagingBucket']
            else:
                self.stagingBucket = None

            if 'targetAsset' in confKeys:
                self.targetAsset = conf['targetAsset']
                # if ~self.targetAsset.endswith('/'):
                #     self.targetAsset += '/'
            else:
                self.targetAsset = None

            if 'workdir' in confKeys:
                self.workdir = conf['workdir']
            else:
                warning.warn("working dir for processing could not be parsed from provided yaml file, setting working dir to current dir",
                            UserWarning)
                self.workdir = './'

            # parse processing specific key information
            if 'hand' in prcsKeys:
                self.hand = prcs['hand']
            else:
                warning.warn("HAND assetID could not be parsed from provided yaml file, setting default global hand model",
                            UserWarning)
                self.hand = 'users/gena/GlobalHAND/30m/hand-5000'

            if 'atms' in prcsKeys:
                self.atmsParams = prcs['atms']

            if 'viirs' in prcsKeys:
                self.viirsParams = prcs['viirs']

            if 'sentinel1' in prcsKeys:
                self.s1Params = prcs['sentinel1']

        else:
            raise ValueError('yamlfile configuration argument must be specified')

        return


    def process(self,product, date,skipPreprocessing=False):
        try:
            ee.Initialize()
        except EEException as e:
            print(e)

        self._parse_config()

        if product in ['sentinel1','atms','viirs','modis']:
            dt = utils.decode_date(date)
            tomorrow = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')

            dateDir = os.path.join(self.workdir,dt.strftime('%Y%m%d'))
            prodDir = os.path.join(dateDir,product)

            geom = ee.Geometry.Rectangle(self.region)

            hand = ee.Image(self.hand)

            if product == 'atms':
                params = self.atmsParams
                paramKeys = list(params.keys())

                collId = self.atmsParams['waterFractionAsset']
                worker = Atms(geom,date,tomorrow,collectionid=collId)

                if skipPreprocessing == False:
                    if os.path.exists(dateDir) != True:
                        os.mkdir(dateDir)
                    if os.path.exists(prodDir) != True:
                        os.mkdir(prodDir)

                    geotiffs = worker.extract(dt,self.region,credentials=self.earthdataLogin,outDir=prodDir,gridding_radius=50000)
                    worker.load(geotiffs,self.stagingBucket,collId)

                if 'seed' in paramKeys:
                    permanentWater = ee.Image(self.atmsParams['seed'])
                else:
                    permanentWater = None

                if 'probablistic' in paramKeys:
                    runProbs = params['probablistic']
                else:
                    runProbs = False

                waterImage = worker.waterMap(hand,permanent=permanentWater,probablistic=runProbs)
                waterImage = waterImage.set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                assetTarget = self.targetAsset + '{0}_bathtub_{1}'.format(product,date.replace('-',''))

            elif (product == 'viirs') or (product == 'modis'):
                today = datetime.datetime.now()

                if (today - dt).days < 5:
                    avail = today - datetime.timedelta(5)
                    raise NotImplementedError('NRT processing for VIIRS or MODIS has not been implemented, please select a date prior to {}'.format(avail))
                else:
                    minDate = (dt - datetime.timedelta(45)).strftime('%Y-%m-%d')
                    maxDate = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')

                    if product == 'modis':
                        worker = Modis(geom,minDate,maxDate,collectionid='MODIS/006/MOD09GA')
                        params = self.viirsParams
                    else:
                        worker = Viirs(geom,minDate,maxDate,collectionid='NOAA/VIIRS/001/VNP09GA')
                        params = self.viirsParams

                    paramKeys = list(params.keys())

                    ls = Landsat(geom,minDate,maxDate,collectionid='LANDSAT/LC08/C01/T1_SR')
                    s2 = Sentinel2(geom,minDate,maxDate,collectionid='COPERNICUS/S2_SR')
                    highRes = ee.ImageCollection(ls.collection.merge(s2.collection))

                    worker.downscale(highRes,target_date=date,windowSize=33,A=0.5)

                    if 'probablistic' in paramKeys:
                        runProbs = params['probablistic']
                        params.pop('probablistic')
                    else:
                        runProbs = False
                        nIters=100
                        probTreshold=0.75

                    waterImage = worker.waterMap(date,hand,probablistic=runProbs,**params)
                    waterImage = waterImage\
                        .set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                    assetTarget = self.targetAsset + '{0}_downscaled_globalOtsu_{1}'.format(product,date.replace('-',''))

            elif product == 'sentinel1':
                previous = (dt + datetime.timedelta(0)).strftime('%Y-%m-%d')
                worker = Sentinel1(geom,previous,tomorrow,collectionid='COPERNICUS/S1_GRD')
                waterImage = worker.waterMap(date).And(hand.lt(30))
                waterImage = waterImage.rename('water')\
                    .set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                assetTarget = self.targetAsset + '{0}_bootstrapOtsu_{1}'.format(product,date.replace('-',''))

            else:
                raise NotImplementedError('select product is currently not implemented, please check back with later versions')

            description = '{0}_water_{1}'.format(product,date)
            geeutils.exportImage(waterImage,geom,assetTarget,description=description)

        else:
            raise NotImplementedError('select product is currently not implemented, please check back with later versions')

        return

    def run_tests(self):
        raise NotImplementedError('test functionality not implemented...please ')
        return

    @staticmethod
    def init_env():
        # authenticate earth engine
        cmd = "earthengine authenticate"
        os.system(cmd)

        # initialize gcloud environment
        cmd = "gcloud init"
        os.system(cmd)

        return


def main():
    fire.Fire(hydrafloods)
    return
