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

from . import utils
from .processing import *

ee.Initialize()

class hydrafloods(object):
    def __init__(self,configuration=None):
        if configuration:
            self.filePath = os.path.dirname(os.path.abspath(__file__))
            yamlFile = configuration

            with open(yamlFile,'r') as stream:
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


            # parse top-level configuration key information
            if 'name' in confKeys:
                self.name = conf['name']
            else:
                raise AttributeError('provided yaml file does not have a name parameter in configuration')

            if 'region' in confKeys:
                self.region = gpd.read_file(conf['region'])

            elif 'country' in confKeys:
                world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
                country = world[world.name == conf['country']]
                if len(country) >= 1:
                    self.region = country
                else:
                    raise ValueError('could not parse selected country from world shapefile')

            elif 'boundingbox' in confKeys:
                from shapely import geometry
                self.region = gpd.GeoDataFrame(pd.DataFrame({'id':[0],'geometry':[geometry.box(*conf['boundingbox'])]}))

            else:
                raise AttributeError('provided yaml file does not have a specified region in configuration')

            if 'credentials' in confKeys:
                self.credentials = conf['credentials']
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
        if product in ['sentinel1','atms','viirs']:
            dt = utils.decode_date(date)
            tomorrow = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')

            dateDir = os.path.join(self.workdir,dt.strftime('%Y%m%d'))
            prodDir = os.path.join(dateDir,product)

            geom = ee.Geometry.Rectangle(list(self.region.bounds.values[0]))

            hand = ee.Image(self.hand)

            if product == 'atms':
                if os.path.exists(dateDir) != True:
                    os.mkdir(dateDir)
                if os.path.exists(prodDir) != True:
                    os.mkdir(prodDir)

                collId = self.atmsParams['waterFractionAsset']
                worker = Atms(geom,date,tomorrow,collectionid=collId)
                params = self.atmsParams
                paramKeys = list(params.keys())

                if skipPreprocessing == False:
                    geotiffs = worker.extract(dt,self.region,outdir=prodDir,creds=self.credentials,gridding_radius=50000)
                    worker.load(geotiffs,self.stagingBucket,collId)

                if 'seed' in list(self.atmsParams.keys()):
                    permanentWater = ee.Image(self.atmsParams['seed'])
                else:
                    permanentWater = None

                waterImage = worker.waterMap(hand,permanent=permanentWater)
                mask = waterImage.select('water')
                waterImage = waterImage.updateMask(mask).set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                assetTarget = self.targetAsset + '{0}_bathtub_{1}'.format(product,date.replace('-',''))
                description= 'ATMS_WATER_'+date

                geeutils.exportImage(waterImage,geom,assetTarget,description=description)

            elif product == 'viirs':
                params = self.viirsParams
                paramKeys = list(params.keys())

            elif product == 'sentinel1':
                # tomorrow = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')
                # nextDay = (dt + datetime.timedelta(-12)).strftime('%Y-%m-%d')
                print(date,tomorrow)
                worker = Sentinel1(geom,date,tomorrow)
                waterImage = worker.waterMap(date,geom).And(hand.lt(30))
                waterImage = waterImage.updateMask(waterImage).rename('water')\
                    .set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                assetTarget = self.targetAsset + '{0}_logitTransform_{1}'.format(product,date.replace('-',''))
                description= 'SENTINEL1_WATER_'+date
                geeutils.exportImage(waterImage,geom,assetTarget,description=description)

            else:
                raise NotImplementedError('select product is currently not implemented, please check back with later versions')

        else:
            raise NotImplementedError('select product is currently not implemented, please check back with later versions')

        return

    def run_tests(self):
        raise NotImplementedError('test functionality not implemented...please ')
        return

def main():
    fire.Fire(hydrafloods)
    return
