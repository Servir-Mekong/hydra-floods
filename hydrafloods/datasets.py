import ee
import numpy as np
import xarray as xr

import rastersmith as rs

from . import geeutils, downscale

class Sentinel1(object):
    def __init__(self,gr,
                 canny_threshold=7,    # threshold for canny edge detection
                 canny_sigma=1,        # sigma value for gaussian filter
                 canny_lt=7,           # lower threshold for canny detection
                 smoothing=100,        # amount of smoothing in meters
                 connected_pixels=200, # maximum size of the neighborhood in pixels
                 edge_length=50,       # minimum length of edges from canny detection
                 smooth_edges=100):


        self.gr = gr
        self.canny_threshold = canny_threshold
        self.canny_sigma = canny_sigma
        self.canny_lt = canny_lt
        self.smoothing = smoothing
        self.connected_pixels = connected_pixels
        self.edge_length = edge_length
        self.smooth_edges = smooth_edges

        return

    def getFloodMap(time_start,time_end):

        geom = ee.Geometry.Rectangle([gr.west,gr.south,gr.east,gr.north])

        mapResult = geeutils.s1WaterMap(geom,time_start,time_end,self.canny_threshold,
                                        self.canny_sigma,self.canny_lt,self.smoothing,
                                        self.connected_pixels,self.edge_length,
                                        self.smooth_edges)

        return mapResult


class Atms(object):
    def __init__():

        return


class Landsat8(object):
    def __init__():
        return


class Viirs(object):
    def __init__():
        return


class Sentinel2(object):
    def __init__():
        return
