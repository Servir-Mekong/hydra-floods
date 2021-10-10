import ee
from functools import partial

from hydrafloods import filtering


def fwdet(water_img, dem, band=None, outlier_test=False, force_projection=False):
    """Implementation of the Flood Water Depth Estimation Tool
    Used to calculate water depth given a water map and elevation data
    Original paper: https://doi.org/10.5194/nhess-19-2053-2019
    Earth Engine paper: https://doi.org/10.1109/LGRS.2020.3031190

    args:
        water_img (ee.Image): earth engine image object representing the water extent. Assumes values of 1 equal water.
        dem (ee.Image): earth engine image object representing elevation
        band (str | None): band name to use for algorithm, if set to `None` will use first band in image. default = None
        outlier_test (bool): flag used to find and filter outlier dem values around water boundaries. default = False
        force_projection (bool): flag used to force the output image projection to that of the input dem. default = False

    returns:
        ee.Image: image of water depth in meters
    """

    proj = dem.projection()

    # select band to use for algorithm
    if band is None:
        water_img = water_img.select([0]).selfMask()

    else:
        water_img = water_img.select(ee.String(band)).selfMask()

    if outlier_test:
        filled = filtering.modified_median_zscore(dem)

        expand = water_img.focal_max(
            kernel=ee.Kernel.square(radius=res, units="meters")
        )

        dem_mask = dem.updateMask(water_image.gt(0))

        boundary = dem_mask.add(expand)

        filled = hf.modified_median_zscore(boundary, filled)

    else:
        filled = dem

    # cumulative cost model

    mod = filled.updateMask(water_img.mask().eq(0))
    source = mod.mask()

    val = 10000
    push = 5000

    cost0 = ee.Image(val).where(source, 0).cumulativeCost(source, push)
    cost1 = ee.Image(val).where(source, 1).cumulativeCost(source, push)
    cost2 = mod.unmask(val).cumulativeCost(source, push)

    costFill = cost2.subtract(cost0).divide(cost1.subtract(cost0))

    costSurface = mod.unmask(0).add(costFill)

    # Kernel for low-pass filter
    boxcar = ee.Kernel.square(radius=3, units="pixels", normalize=True)

    # Floodwater depth calculation and smoothing using a low-pass filter
    depths = (
        costSurface.subtract(filled).convolve(boxcar).updateMask(water_img.mask())
    ).rename("depth")

    # force min values to be 0
    depths = depths.max(ee.Image.constant(0))

    # force the output project to be that of the input dem if force_projection is true
    if force_projection:
        depths = depths.reproject(proj)

    return depths


def fwdet_experimental(
    water_img,
    dem,
    iter=10,
    band=None,
    pixel_edge=0,
    boundary_definition="mean",
    smooth_depths=True,
    force_projection=False,
):
    """Experimantal changes to the implementation of the Flood Water Depth Estimation Tool
    Used to calculate water depth given a water map and elevation data

    Difference from original FwDET implementation include ...

    args:
        water_img (ee.Image): earth engine image object representing the water extent.
        dem (ee.Image): earth engine image object representing elevation
        iter (int|ee.Number, optional): keyword to set number of iterations to fill up every pixel.
            The algorithm will calculate minimum iterations needed and use the max between the two.
            default = 10
        band (str | None,optional): band name to use for algorithm, if set to `None` will use first band in image. default = None
        pixel_edge (int): value that represents the boundary/edge of water. Note: must be client-side int object. default = 0
        boundary_definition (str): method for defining the elevation at water edges. Options are 'mean' or 'nearest'
            Note: 'nearest' is the FwDET original method, 'mean' is an updated method. default=mean
        smooth_depths (bool): flag to control if the resulting water depths should be smoothed or not, uses a 3x3 pixel mean. default=True
        force_projection (bool): flag used to force the output image projection to that of the input dem. default = False


    returns:
        ee.Image: image of water depth in meters
    """
    proj = dem.projection()

    res = water_img.projection().nominalScale().multiply(1.5)

    # select band to use for algorithm
    if band is None:
        water_img = water_img.select([0])

    else:
        water_img = water_img.select(ee.String(band))

    # detect water edge
    watermap_edge = (
        water_img.selfMask().unmask(-999).focal_min(res, "square", "meters").eq(-999)
    )
    watermap_edge = watermap_edge.updateMask(water_img.unmask(0))
    watermap_edge = watermap_edge.selfMask()
    # watermap_edge = watermap_edge.reproject(water_img.projection());

    # detect watermap extent outer edge (not to be used in algorithm)
    if pixel_edge > 0:
        outer_edge = (
            water_img.mask()
            .unmask(0)
            .gte(0)
            .unmask(0, False)
            .gt(0)
            .focal_min(res.multiply(ee.Number(pixel_edge).add(1)), "square", "meters")
            .eq(0)
        )
        outer_edge = outer_edge.updateMask(outer_edge)
        # outer_edge = outer_edge.reproject(water_img.projection())
        # update watermap edge
        watermap_edge = watermap_edge.updateMask(outer_edge.unmask(0).eq(0))
        # watermap_edge = watermap_edge.reproject(water_img.projection());
    else:
        outer_edge = water_img.mask().unmask(0).lt(0)

    DEM = dem  # being super lazy here by setting a variable vs renaming things...

    # get elevation at edge
    # DEM_watered_edge = watermap_edge.multiply(DEM); # simply using boundary cell elevations
    DEM_masked = DEM.updateMask(
        water_img.unmask(0).eq(0).add(watermap_edge.unmask(0)).eq(1)
    )
    DEM_watered_edge = watermap_edge.multiply(
        DEM_masked.focal_mean(res, "square", "meters")
    )
    DEM_watered_edge = DEM_watered_edge  # .reproject(water_img.projection())
    # calculate (approximation of) minimum needed iterations to fill up every pixel
    req_iter = (
        water_img.eq(0)
        .add(watermap_edge)
        .fastDistanceTransform(50, "pixels", "squared_euclidean")
        .sqrt()
        .round()
        .updateMask(water_img)
        .updateMask(outer_edge.unmask(0).eq(0))
        .reproject(water_img.projection())
    )
    req_iter = req_iter.reduceRegion(
        reducer=ee.Reducer.max(),
        geometry=water_img.geometry().bounds(1e3),
        scale=res,
        maxPixels=1e12,
    ).get("distance")
    # set iterations to be used
    iter = ee.Number(req_iter).max(iter)
    # construct list to iterate over (list contents unused, just for iterations)
    iter_list = ee.List.sequence(0, iter)

    # use look up which algorithm to use for defining the boundary elevation values
    boundary_algos = {"mean": _mbce, "nearest": _nbce}
    f_boundary = partial(
        boundary_algos[boundary_definition], mask=water_img, resolution=res
    )

    # apply boundary condition algo, starting with edge values and filling it up further with each iteration
    DEM_watered_fill = ee.Image(iter_list.iterate(f_boundary, DEM_watered_edge))
    DEM_watered = DEM_watered_fill.reproject(water_img.projection())
    # derive water depths
    depths = DEM_watered.subtract(DEM).max(0)

    # smooth water depths
    # Cohen et al. (2018) use a 3x3 pixel mean (but FwDET-GEE seems to use a boxcar kernel?)
    if smooth_depths:
        depths = depths.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.square(res, units="meters"),
            optimization="boxcar",  # optimization for mean.
        )

    # force the output project to be that of the input dem if force_projection is true
    if force_projection:
        depths = depths.reproject(proj)

    return depths.rename("depth").set("depth_iter", iter)


# Mean Boundary Cell Elevation (MBCE), adapted from NBCE above
# MBCE function
def _mbce(i, img, mask=None, resolution=30):
    """Mean Boundary Cell Elevation (MBCE), unpublished A. Haag (2021)"""
    # obtain mean boundary cell elevation
    new_img = ee.Image(img).focal_mean(resolution, "square", "meters")
    # make sure the original / previous step values are kept intact (comment out for more smoothing)
    # new_img = new_img.where(ee.Image(img).gt(-999), ee.Image(img));
    # make sure it does not expand to other side of floodmap (where there is no water)
    return new_img.updateMask(mask.unmask(0))  # .reproject(floodmap.projection())


# Nearest Boundary Cell Elevation (NBCE) algorithm
def _nbce(i, img, mask=None, resolution=30):
    """Nearest Boundary Cell Elevation (NBCE), Cohen et al. (2018)"""
    # obtain nearest boundary cell elevations using different kernels
    new_img_0 = ee.Image(img).focal_min(resolution, "square", "meters")
    new_img_1 = ee.Image(img).focal_min(resolution, "plus", "meters")
    new_img_2 = ee.Image(img).focal_min(resolution, "cross", "meters")
    # make sure the original / previous step values are kept intact
    new_img_0 = new_img_0.where(ee.Image(img).gt(-999), ee.Image(img))
    new_img_1 = new_img_1.where(ee.Image(img).gt(-999), ee.Image(img))
    new_img_2 = new_img_2.where(ee.Image(img).gt(-999), ee.Image(img))
    # add the NBCE in right order (minimum distance at the end, overwrites diagonal)
    new_img = new_img_0.where(new_img_2.gt(-999), new_img_2)
    new_img = new_img.where(new_img_1.gt(-999), new_img_1)
    # make sure it does not expand to other side of floodmap (where there is no water)
    return new_img.updateMask(mask.unmask(0))


def downscale_wfraction(fraction, dem):
    """Algorithm to downscale water fraction maps using DEM data
    Converts fractional water data to higher resolution binary water data
    Paper: https://doi.org/10.1016/j.rse.2013.03.015

    args:
        fraction (ee.Image): earth engine image object representing the fractional water extent within each pixel.
        dem (ee.Image): earth engine image object representing elevation

    returns:
        ee.Image: downscaled binary image of water
    """

    resolution = fraction.projection().nominalScale()

    water_present = fraction.gt(0.0)

    h_min = dem.updateMask(water_present).focal_min(resolution, "square", "meters")
    h_max = dem.updateMask(water_present).focal_max(resolution, "square", "meters")

    water_high = h_min.add(h_max.subtract(h_min).multiply(fraction))

    return dem.gte(h_min).And(dem.lt(water_high))
