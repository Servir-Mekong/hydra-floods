# coding: utf-8 VIIRS DNNS

import ee
import numpy as np
import math
import matplotlib.pyplot as pl


def perm_water_mask():
    return ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24").select(["water_mask"], ["b1"])


def DEM():
    return ee.Image("CGIAR/SRTM90_V4")


# ALOS DEM
# ee.Image('JAXA/ALOS/AW3D30_V1_1')
# SRTM DEM
# ee.Image('USGS/SRTMGL1_003')


# def GEE_classifier(img, name_classifier, arguments={}):
#     arguments = {
#         "training_image": perm_water_mask(),
#         "training_band": "b1",
#         "training_region": img,
#     }
#     c_arguments = {
#         "subsampling": 0.1,
#         "max_classification": 2,
#         "classifier_mode": "classification",
#         "name_classifier": name_classifier,
#     }
#     arguments.update(c_arguments)
#     classes1 = img.trainClassifier(**arguments)
#     classes2 = classify(classes1).select(["classification"], ["b1"])
#     return classes2


def dnns(img):

    kernel_scale = 40
    VIIRS_ave_scale = 500

    k1 = ee.Kernel.square(kernel_scale, "pixels", False)

    composite = img["I1"].addBands(img["I2"]).addBands(img["I3"])

    # classes = GEE_classifier(img, "Pegasos", {"classifier_mode": "probability"})
    classes = ee.Image([0.5, 0.5])
    water = classes.gte(0.9)
    land = classes.lte(0.1)
    mix = water.Not().And(land.Not())
    ave_water = water.mask(water).multiply(composite)
    ave_water_img = ee.Image([ave_water["I1"], ave_water["I2"], ave_water["I3"]])

    N_nmin_water = 1
    N_water = water.convolve(k1)

    water_rf = (
        water.multiply(composite)
        .convolve(k1)
        .multiply(N_water.gte(N_nmin_water))
        .divide(N_water)
    )
    water_rf = water_rf.add(ave_water_img.multiply(water_rf.Not()))

    ave_land = composite.mask(land)

    ave_land_img = ee.Image([ave_land["I1"], ave_land["I2"], ave_land["I3"]])

    # Constraints
    R1 = img["I1"].divide(img["I3"])
    R2 = img["I2"].divide(img["I3"])
    R3 = R1.subtract(water_rf.select("I1").divide(img["I3"]))
    R4 = R2.subtract(water_rf.select("I2").divide(img["I3"]))

    NR1 = R1.neighborhoodToBands(k1)
    NR2 = R2.neighborhoodToBands(k1)
    NI1 = img["I1"].neighborhoodToBands(k1)
    NI2 = img["I2"].neighborhoodToBands(k1)
    NI3 = img["I3"].neighborhoodToBands(k1)

    M1 = (NR1.gt(R3)).And(NR1.lt(R1))
    M2 = (NR2.gt(R4)).And(NR2.lt(R2))
    nLP = M1.And(M2)

    NnLP = nLP.reduce(ee.Reducer.sum())
    ave_nI1 = NI1.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)
    ave_nI2 = NI2.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)
    ave_nI3 = NI3.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)

    N_nmin_land = 1
    ave_land = ave_nI1.addBands(ave_nI2).addBands(ave_nI3)
    ave_land = ave_land.multiply(NnLP.gte(N_nmin_land)).add(
        ave_land_img.multiply(NnLP.lt(N_nmin_land))
    )

    ave_landI3 = ave_land.select("I3")
    f_water = (
        (ave_landI3.subtract(img["I3"]))
        .divide(ave_landI3.subtract(water_rf.select("I3")))
        .clamp(0, 1)
    )
    f_water = f_water.add(water).subtract(land).clamp(0, 1)

    return f_water


def DEM_downscale(img, f_water):

    VIIRS_Pixel_meters = 500
    DEM_d = img.DEM()

    water_present = f_water.gt(0.0)

    h_min = DEM_d.mask(f_water).focal_min(VIIRS_Pixel_meters, "square", "meters")
    h_max = DEM_d.mask(f_water).focal_max(VIIRS_Pixel_meters, "square", "meters")

    water_high = h_min.add(h_max.subtract(h_min).multiply(f_water))
    water_high = water_high.multiply(f_water.lt(1.0)).multiply(f_water.gt(0.0))
    return f_water.eq(1.0).select(["elevation"], ["I1"])
