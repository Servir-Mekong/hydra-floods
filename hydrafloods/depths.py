import ee
from functools import partial


def fwdet(
    water_img,
    dem,
    iter=10,
    band="water",
    pixel_edge=0,
    boundary_definition="mean",
    smooth_depths=True,
):
    proj = dem.projection()

    res = water_img.projection().nominalScale().multiply(1.5)

    # do we need proejction here
    water_img = water_img.select(band)  # .reproject(proj);

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
            .focal_min(res.multiply(ee.Number(pix_edge).add(1)), "square", "meters")
            .eq(0)
        )
        outer_edge = outer_edge.updateMask(outer_edge)
        outer_edge = outer_edge.reproject(water_img.projection())
        # update watermap edge
        watermap_edge = watermap_edge.updateMask(outer_edge.unmask(0).eq(0))
        # watermap_edge = watermap_edge.reproject(water_img.projection());
    else:
        outer_edge = water_img.mask().unmask(0).lt(0)

    DEM = dem

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

    # apply MBCE, starting with edge values and filling it up further with each iteration
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

    return depths.rename("depth").set("depth_iter", iter)


# Mean Boundary Cell Elevation (MBCE), adapted from NBCE above
# MBCE function
def _mbce(i, img, mask=None, resolution=30):
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
