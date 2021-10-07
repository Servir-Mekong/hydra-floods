import ee


def fwdet(
    water_img, dem, iter=10, band="water", pixel_edge=0, boundary_definition="mb"
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
    DEM_watered_edge = DEM_watered_edge.reproject(water_img.projection())
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
    # MBCE function
    def MBCE_helper(i, img):
        # obtain mean boundary cell elevation
        new_img = ee.Image(img).focal_mean(res, "square", "meters")
        # make sure the original / previous step values are kept intact (comment out for more smoothing)
        # new_img = new_img.where(ee.Image(img).gt(-999), ee.Image(img));
        # make sure it does not expand to other side of watermap (where there is no water)
        return new_img.updateMask(water_img.unmask(0)).reproject(water_img.projection())

    # apply MBCE, starting with edge values and filling it up further with each iteration
    DEM_watered_fill = ee.Image(iter_list.iterate(MBCE_helper, DEM_watered_edge))
    DEM_watered_MBCE = DEM_watered_fill.reproject(water_img.projection())
    # derive water depths
    depths_MBCE = DEM_watered_MBCE.subtract(DEM).max(0)
    # smooth water depths
    # Cohen et al. (2018) use a 3x3 pixel mean (but FwDET-GEE seems to use a boxcar kernel?)
    # depths_smoothed_MBCE = depths_MBCE.focal_mean(res.multiply(3), 'square', 'meters');
    # depths_smoothed_MBCE = depths_smoothed_MBCE.updateMask(water_img.unmask(0));
    # return ee.Dictionary({
    #   'MBCE': depths_MBCE,
    #   'MBCE_smoothed': depths_smoothed_MBCE
    # });
    # return depths_MBCE;
    # return depths_smoothed_MBCE;
    # return water_img.addBands(depths_MBCE.rename(["depth_MBCE"])).set(
    #     "depth_iter", iter
    # )
    return depths_MBCE.rename("depth").set("depth_iter", iter)
    # , 'depth_method', 'MBCE');


# Mean Boundary Cell Elevation (MBCE), adapted from NBCE above
def mbce(floodmap, DEM_flooded_edge, res, iter):
    # construct list to iterate over (list contents unused, just for iterations)
    iter_list = ee.List.sequence(0, iter)
    # MBCE function
    def MBCE(i, img):
        # obtain mean boundary cell elevation
        new_img = ee.Image(img).focal_mean(res, "square", "meters")
        # make sure the original / previous step values are kept intact (comment out for more smoothing)
        # new_img = new_img.where(ee.Image(img).gt(-999), ee.Image(img));
        # make sure it does not expand to other side of floodmap (where there is no water)
        return new_img.updateMask(floodmap.unmask(0)).reproject(floodmap.projection())

    # apply MBCE, starting with edge values and filling it up further with each iteration
    DEM_flooded_fill = ee.Image(iter_list.iterate(MBCE, DEM_flooded_edge))
    return DEM_flooded_fill.reproject(floodmap.projection())


# Nearest Boundary Cell Elevation (NBCE) algorithm
def nbce(floodmap, DEM_flooded_edge, res, iter):
    # construct list to iterate over (list contents unused, just for iterations)
    iter_list = ee.List.sequence(0, iter)
    # NBCE function
    def NBCE(i, img):
        # obtain nearest boundary cell elevations using different kernels
        new_img_0 = ee.Image(img).focal_min(res, "square", "meters")
        new_img_1 = ee.Image(img).focal_min(res, "plus", "meters")
        new_img_2 = ee.Image(img).focal_min(res, "cross", "meters")
        # make sure the original / previous step values are kept intact
        new_img_0 = new_img_0.where(ee.Image(img).gt(-999), ee.Image(img))
        new_img_1 = new_img_1.where(ee.Image(img).gt(-999), ee.Image(img))
        new_img_2 = new_img_2.where(ee.Image(img).gt(-999), ee.Image(img))
        # add the NBCE in right order (minimum distance at the end, overwrites diagonal)
        new_img = new_img_0.where(new_img_2.gt(-999), new_img_2)
        new_img = new_img.where(new_img_1.gt(-999), new_img_1)
        # make sure it does not expand to other side of floodmap (where there is no water)
        return new_img.updateMask(floodmap.unmask(0))

    # apply NBCE, starting with edge values and filling it up further with each iteration
    DEM_flooded_fill = ee.Image(iter_list.iterate(NBCE, DEM_flooded_edge))
    return DEM_flooded_fill
