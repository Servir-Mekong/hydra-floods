# coding: utf-8 MODIS DNNS approach

import ee
import numpy as np
import math
import matplotlib.pyplot as pl

def perm_water_mask():
    return ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24").select(['water_mask'], ['b1'])

def DEM():
    return ee.Image('CGIAR/SRTM90_V4')
#ALOS DEM
#ee.Image('JAXA/ALOS/AW3D30_V1_1') 
#SRTM DEM
#ee.Image('USGS/SRTMGL1_003')


def GEE_classifier(img, name_classifier, arguments={}):
    arguments = {
                'training_image'    : perm_water_mask(),
                'training_band'     : "b1",
                'training_region'   : img
               }
    c_arguments = {
                   'subsampling'       : 0.1, 
                   'max_classification': 2,
                   'classifier_mode' : 'classification',
                   'name_classifier'   : name_classifier
                  }
    arguments.update(c_arguments)
    classes1 = img.trainClassifier(**arguments)  
    classes2 = classify(classes1).select(['classification'], ['b1'])
    return classes2;


def dnns(img):
    
    kernel_scale = 40
    VIIRS_ave_scale = 500 
    
    k1 = ee.Kernel.square(kernel_scale, 'pixels', False)
    
    composite = img['b1'].addBands(img['b2']).addBands(img['b6'])
    
    classes = GEE_classifier(img, 'Pegasos', {'classifier_mode' : 'probability'})
    water = classes.gte(0.9)
    land  = classes.lte(0.1)
    mix = water.Not().And(land.Not())
    ave_water = water.mask(water).multiply(composite)
    ave_water_img = ee.Image([ave_water['b1'], ave_water['b2'], ave_water['b6']])
    
    N_nmin_water = 1
    N_water = water.convolve(k1)
    
    water_rf = water.multiply(composite).convolve(k1).multiply(N_water.gte(N_nmin_water)).divide(N_water)
    water_rf = water_rf.add(ave_water_img.multiply(water_rf.Not()))

    ave_land = composite.mask(land)
    
    ave_land_img = ee.Image([ave_land['b1'], ave_land['b2'], ave_land['b6']])
    
    # Constraints 
    R1 = img['b1'].divide(img['b6'])
    R2 = img['b2'].divide(img['b6'])
    R3 = R1.subtract(water_rf.select('b1').divide(img['b6']) )
    R4 = R2.subtract(water_rf.select('b2').divide(img['b6']) )
    
    NR1 = R1.neighborhoodToBands(k1) 
    NR2 = R2.neighborhoodToBands(k1)
    NI1 = img['b1'].neighborhoodToBands(k1)
    NI2 = img['b2'].neighborhoodToBands(k1)
    NI3 = img['b6'].neighborhoodToBands(k1)
    
    
    M1 = (NR1.gt(R3)).And(NR1.lt(R1))
    M2 = (NR2.gt(R4)).And(NR2.lt(R2))
    nLP = M1.And(M2)
    
    NnLP = nLP.reduce(ee.Reducer.sum())
    ave_nI1 = NI1.multiply(nLP).reduce(ee.Reducer.sum()).divide(numnLP)
    ave_nI2 = NI2.multiply(nLP).reduce(ee.Reducer.sum()).divide(numnLP)
    ave_nI3 = NI3.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)

    N_nmin_land = 1
    ave_land = ave_nI1.addBands(ave_nI2).addBands(ave_nI3)
    ave_land = ave_land.multiply(NnLP.gte(N_nmin_land)).add(ave_land_img.multiply(NnLP.lt(N_nmin_land)) )

    ave_landI3 = ave_land.select('b6')
    f_water = (ave_landI3.subtract(img['b6'])).divide(ave_landI3.subtract(water_rf.select('b6'))).clamp(0, 1)
    f_water = f_water.add(water).subtract(land).clamp(0, 1)
    
    return f_water


def DEM_downscale(img, f_water):
    
    MODIS_Pixel_meters = 500
    DEM_d = img.DEM()
    
    water_present = f_water.gt(0.0)

    h_min = DEM_d.mask(f_water).focal_min(MODIS_Pixel_meters, 'square', 'meters')
    h_max = DEM_d.mask(f_water).focal_max(MODIS_Pixel_meters, 'square', 'meters')
    
    
    water_high = h_min.add(h_max.subtract(h_min).multiply(f_water))
    water_high = water_high.multiply(f_water.lt(1.0)).multiply(f_water.gt(0.0)) 
    return f_water.eq(1.0).select(['elevation'], ['b1'])

