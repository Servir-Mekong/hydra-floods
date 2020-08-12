from __future__ import print_function, division
import ee
from hydrafloods import decorators


def starfm(
    coarseCollection, fineCollection=None, targetDate="1970-01-01", windowSize=33, A=0.5
):
    @decorators.carry_metadata
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
