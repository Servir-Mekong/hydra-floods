from __future__ import print_function, division
import ee
from hydrafloods import decorators, ml, thresholding


def starfm(
    coarseCollection, fineCollection=None, targetDate="1970-01-01", windowSize=33, A=0.5
):
    @decorators.keep_attrs
    def apply_starfm(img):
        t = ee.Date(img.get("system:time_start"))

        Tijk = img.select("time").subtract(base.select("time_first"))

        Sijk = img.subtract(slv).abs().convolve(Dijk).rename(bandList)

        Cijk = Sijk.multiply(Tijk).convolve(Dijk)
        Wijk = one.divide(Cijk).divide(one.divide(Cijk).convolve(Dijk))

        nBands = ee.Number(base.bandNames().length())
        expected = ee.Number(len(bandList))

        outImg = ee.Algorithms.If(
            nBands.neq(expected),
            None,
            base.divide(Wijk.convolve(Dijk))
            .add(Sijk)
            .set("system:time_start", t)
            .rename(bandList),
        )

        return ee.Image(outImg)

    target = ee.Date(targetDate)

    if type(fineCollection) is not ee.imagecollection.ImageCollection:
        try:
            fineCollection = getattr(fineCollection, "collection")
        except Exception as e:
            raise TypeError(
                "keyword fineCollection needs to be either of type ee.ImageCollection "
                "or hf.hfCollection"
            )

    bandList = ee.Image(coarseCollection.first()).bandNames().getInfo()

    one = ee.Image.constant(1)
    centerPos = ee.Number((windowSize - 1) / 2)
    pos = ee.Array(ee.List.sequence(0, windowSize - 1))
    xPos = pos.repeat(1, windowSize)
    yPos = pos.repeat(1, windowSize).transpose()

    dijk = (
        ee.Array(centerPos)
        .subtract(xPos)
        .pow(2)
        .add(ee.Array(centerPos).subtract(yPos).pow(2))
        .sqrt()
    )

    weightMax = dijk.get([0, 0]).add(1)
    # flip weights to have heighest in middle
    # normalize values for convolutions
    dijk_norm = dijk.subtract(weightMax).abs().divide(weightMax.subtract(1))

    dW = ee.Array(1).add(dijk.divide(ee.Array(A)))

    square = ee.Kernel.square(windowSize, "pixels")
    Dijk = ee.Kernel.fixed(
        windowSize, windowSize, dW.toList(), centerPos, centerPos, True
    )

    base = (
        fineCollection.filterDate(target.advance(-2, "month"), target)
        .sort("system:time_start", False)
        .reduce(ee.Reducer.firstNonNull())
    )

    slv = (
        coarseCollection.filterDate(target.advance(-2, "month"), target)
        .sort("system:time_start", False)
        .reduce(ee.Reducer.firstNonNull())
    )

    result = coarseCollection.filterDate(
        target.advance(-5, "day"), target.advance(5, "day")
    ).map(apply_starfm, True)

    return result


def bathtub(wfrac, hand, permanent=None):
    """
    Function to fit hi-resolution HAND model to water fraction estimate
    args:
        wFrac (ee.Image): water fraction image, values must be 0-1
        hand (ee.Image): height above nearest drainage (HAND) image, units in meters
        permanent (ee.Image): permanent water image to seed HAND filling
    """

    def fillGrids(d):
        d = ee.Number(d)
        dMap = (
            hand.lte(d)
            .Or(permWater)
            .reduceResolution(
                reducer=ee.Reducer.mean(), bestEffort=True, maxPixels=1024
            )
            .reproject(crs=proj)
        )

        diff = wfrac.subtract(dMap).abs()
        return diff

    def minimizeDepth(d):
        d = ee.Number(d)
        dMap = hand.lte(d).Or(permWater)
        rd = dMap.reduceResolution(
            reducer=ee.Reducer.mean(), bestEffort=True, maxPixels=1024
        ).reproject(crs=proj)

        diff = wfrac.subtract(rd).abs()

        out = nil.where(diff.eq(minDiff), dMap)
        err = nil.where(diff.eq(minDiff), diff)
        return out.rename("water").addBands(err.rename("error"))

    nil = ee.Image(0)
    if permanent:
        permWater = permanent
    else:
        permWater = ee.Image(0)
    proj = wfrac.projection()

    depths = ee.List.sequence(0, 20)

    depthMaps = ee.ImageCollection(depths.map(fillGrids))
    minDiff = depthMaps.min()

    waterMap = ee.ImageCollection(depths.map(minimizeDepth))
    final = (
        waterMap.select("water").max().addBands(waterMap.select("error").min().float())
    )

    return final


@decorators.keep_attrs
def pca_fusion(
    img,
    img_stats=None,
    match_stats=None,
    eigen_vecs=None,
    scaling_dict=None,
    img_band=None,
    match_band=None,
    inverse_relationship=True,
    reverse_pca_bands=False,
):
    pca_bands = [img_band, match_band]

    if reverse_pca_bands:
        pca_bands = pca_bands[::-1]

    img_z = (
        img.select(img_band)
        .subtract(img_stats.select(img_band + "_mean"))
        .divide(img_stats.select(img_band + "_stdDev"))
    )

    if inverse_relationship:
        img_z = img_z.multiply(-1)

    match_synth = (
        match_stats.select(match_band + "_mean")
        .add(img_z.multiply((match_stats.select(match_band + "_stdDev"))))
        .rename(match_band)
    )

    if reverse_pca_bands:
        input_img = match_synth.addBands(img)
    else:
        input_img = img.addBands(match_synth)

    if scaling_dict is not None:
        input_img = ml.standard_image_scaling(input_img, scaling_dict, pca_bands)

    pcs = ml.apply_image_pca(input_img, eigen_vecs, pca_bands)

    df_index = pcs.select([0]).clamp(-3.75, 3.75).rename("df_index")

    return df_index

@decorators.keep_attrs
def dnns(img, region=None):
    """Dynamic Nearest Neighbor Search algorithm for estimating 
    
    """
    pimg = img.select(["red","nir","swir1"])

    kernel_size = 40  #pixels?
    scale = pimg.projection().nominalScale()

    kernel = ee.Kernel.square(kernel_size, "pixels", False)
    kernel_norm = ee.Kernel.square(kernel_size, 'pixels', True)

    water = img.select("mndwi").gt(0.15)#hf.bmax_otsu(img, band="mndwi", initial_threshold=0.1, region=region).Not()
    land = img.select("mndwi").lt(-0.15)

    mix = water.Not().And(land.Not())
    ave_water = water.multiply(pimg).reduceRegion(ee.Reducer.mean(),region,scale,maxPixels=1e4,bestEffort=True).toImage()

    N_nmin_water = 1
    N_water = water.convolve(kernel)

    water_rf = (
        water.multiply(pimg)
        .convolve(kernel)
        .multiply(N_water.gte(N_nmin_water))
        .divide(N_water)
    )
    water_rf = water_rf.add(ave_water.multiply(water_rf.Not()))

    ave_land = pimg.multiply(land).reduceRegion(ee.Reducer.mean(),region,scale,maxPixels=1e4,bestEffort=True).toImage()

    R1 = pimg.select("red").divide(pimg.select("swir1"))
    R2 = pimg.select("nir").divide(pimg.select("swir1"))
    R3 = R1.subtract(water_rf.select("red").divide(pimg.select("swir1")))
    R4 = R2.subtract(water_rf.select("red").divide(pimg.select("swir1")))

    NR1 = R1.neighborhoodToBands(kernel)
    NR2 = R2.neighborhoodToBands(kernel)
    NI1 = img.select("red").neighborhoodToBands(kernel)
    NI2 = img.select("nir").neighborhoodToBands(kernel)
    NI3 = img.select("swir1").neighborhoodToBands(kernel)

    M1 = (NR1.gt(R3)).And(NR1.lt(R1))
    M2 = (NR2.gt(R4)).And(NR2.lt(R2))
    nLP = M1.And(M2)

    NnLP = nLP.reduce(ee.Reducer.sum())
    ave_nI1 = NI1.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)
    ave_nI2 = NI2.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)
    ave_nI3 = NI3.multiply(nLP).reduce(ee.Reducer.sum()).divide(NnLP)

    N_nmin_land = 1
    ave_pureland = ee.Image.cat([ave_nI1,ave_nI2,ave_nI3])
    ave_pureland = ave_pureland.multiply(NnLP.gte(N_nmin_land)).add(
        ave_land.multiply(NnLP.lt(N_nmin_land))
    )

    ave_landI3 = ave_land.select([2])
    f_water_all = (
        (ave_pureland.subtract(pimg))
        .divide(ave_pureland.subtract(water_rf))
        .clamp(0, 1)
    )
    f_water = f_water_all.add(water).subtract(land).clamp(0, 1)

    return f_water.reduce(ee.Reducer.mean()).rename("water_fraction").toFloat()
