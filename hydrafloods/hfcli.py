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
                    struct = yaml.load(stream)
                except yaml.YAMLError as exc:
                    print(exc)

            conf = struct['configuration']
            self.conf = conf
            ftch = struct['download']
            self.ftch = ftch
            prcs = struct['process']
            self.prcs = prcs

            confKeys = list(conf.keys())
            ftchKeys = list(ftch.keys())
            prcsKeys = list(prcs.keys())

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

            if 'credentials' in ftchKeys:
                self.credentials = ftch['credentials']
            else:
                self.credentials = None

            if 'outdir' in ftchKeys:
                self.ftchOut = ftch['outdir']
            else:
                warning.warn("output dir for downloading could not be parsed from provided yaml file, setting output dir to current dir",
                            UserWarning)
                self.ftchOut = './'
            if 'viirs' in ftchKeys:
                self.viirsFetch = ftch['viirs']
            if 'modis' in ftchKeys:
                self.modisFetch = ftch['modis']

            if 'atms' in prcsKeys:
                self.atmsParams = prcs['atms']
            if 'landsat' in prcsKeys:
                self.landsatParams = prcs['landsat']
            if 'viirs' in prcsKeys:
                self.viirsParams = prcs['viirs']

            if 'outdir' in prcsKeys:
                self.prcsOut = prcs['outdir']
            else:
                warning.warn("output dir for processing could not be parsed from provided yaml file, setting output dir to current dir",
                            UserWarning)
                self.prcsOut = './'

        else:
            raise ValueError('yamlfile configuration argument must be specified')

        return


    def process(self,product, date):
        if product in ['sentinel1','atms','viirs']:
            dt = utils.decode_date(date)

            dateDir = os.path.join(self.ftchOut,dt.strftime('%Y%m%d'))
            prodDir = os.path.join(dateDir,product)

            geom = ee.Geometry.Rectangle(list(self.region.bounds.values[0]))

            if product == 'atms':
                if os.path.exists(dateDir) != True:
                    os.mkdir(dateDir)
                if os.path.exists(prodDir) != True:
                    os.mkdir(prodDir)

                nextDay = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')

                collId = 'projects/servir-mekong/hydrafloods/atms_waterfraction'
                worker = Atms(geom,date,nextDay,collectionid=collId)
                params = self.atmsParams
                paramKeys = list(params.keys())

                geotiffs = worker.extract(dt,self.region,outdir=prodDir,creds=self.credentials,gridding_radius=50000)
                worker.load(geotiffs,'gs://servirmekong/hydrafloods/atms/',collId)

                permanentWater = ee.Image('projects/servir-mekong/yearly_primitives_smoothed/water/water2017').gt(50)
                hand = ee.Image('users/arjenhaag/SERVIR-Mekong/HAND_MERIT')

                waterImage = worker.waterMap(hand,permanent=permanentWater)
                mask = waterImage.select('water')
                waterImage = waterImage.updateMask(mask).set({'system:time_start':ee.Date(date).millis(),'sensor':product})
                assetTarget = 'projects/servir-mekong/hydrafloods/surface_water/' + '{0}_bathtub_{1}'.format(product,date.replace('-',''))

                geeutils.exportImage(waterImage,geom,assetTarget)

            elif product == 'viirs':
                params = self.viirsParams
                paramKeys = list(params.keys())

            elif product == 'sentinel1':
                tomorrow = (dt + datetime.timedelta(1)).strftime('%Y-%m-%d')
                nextDay = (dt + datetime.timedelta(-12)).strftime('%Y-%m-%d')
                worker = Sentinel1(geom,nextDay,tomorrow)
                waterImage = worker.waterMap(date)
                waterImage = waterImage.set({'system:time_start':ee.Date(date).millis(),'sensor':product}).rename('water')
                assetTarget = 'projects/servir-mekong/hydrafloods/surface_water/' + '{0}_bootStrapOtsu_{1}'.format(product,date.replace('-',''))
                geeutils.exportImage(waterImage,geom,assetTarget)

            else:
                raise NotImplementedError('select product is currently not implemented, please check back with later versions')

        else:
            raise NotImplementedError('select product is currently not implemented, please check back with later versions')

        return

def main():
    fire.Fire(hydrafloods)
    return
