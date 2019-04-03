from __future__ import print_function,division
import ee


def starfm(fineCollection,coarseCollection,target_date='1970-01-01',windowSize=33,A=0.5):
    def apply_starfm(img):
        t = ee.Date(img.get('system:time_start'))

        base = fineCollection.filterDate(t.advance(-2,'month'),t)\
                .sort('system:time_start',False).reduce(ee.Reducer.firstNonNull())

        slv =  coarseCollection.filterDate(t.advance(-2,'month'),t)\
                .sort('system:time_start',False).reduce(ee.Reducer.firstNonNull())

        Tijk = img.select('time').subtract(base.select('time_first'))

        Sijk = img.subtract(slv).abs().reduceNeighborhood(ee.Reducer.sum(),square,'kernel',True,'boxcar')\
                .rename(bandRemap.get('new'))

        Cijk = Sijk.multiply(Tijk).convolve(Dijk)
        Wijk = one.divide(Cijk).divide(one.divide(Cijk)\
                  .reduceNeighborhood(ee.Reducer.sum(),Dijk))

        nBands = ee.Number(base.bandNames().length())
        expected = ee.Number(ee.List(bandRemap.get('new')).length())

        outImg = ee.Algorithms.If(nBands.neq(expected), None,
                    base.divide(Wijk.reduceNeighborhood(ee.Reducer.sum(),Dijk)).add(Sijk)
                    .set('system:time_start',t)
                    .rename(bandRemap.get('new'))
                    )

        return outImg

    iniTime = ee.Date.fromYMD(1970,1,1)
    target = ee.Date(target_date)

    bandRemap = ee.Dictionary({
      'landsat': ee.List(['B2','B3','B4','B5','B6','B7','time']),
      'viirs': ee.List(['M2','M4','M5','M7','M10','M11','time']),
      'new': ee.List(['blue','green','red','nir','swir1','swir2','time'])
    });

    one = ee.Image.constant(1)
    centerPos = ee.Number((windowSize-1)/2)
    pos = ee.Array(ee.List.sequence(0,windowSize-1))
    xPos = pos.repeat(1,windowSize)
    yPos = pos.repeat(1,windowSize).transpose()

    dijk = ee.Array(centerPos).subtract(xPos).pow(2).add(
           ee.Array(centerPos).subtract(yPos).pow(2)).sqrt()

    dW = ee.Array(1).add(dijk.divide(ee.Array(A)))

    square = ee.Kernel.square(windowSize,'pixels')
    Dijk = ee.Kernel.fixed(windowSize,windowSize,dW.toList(),centerPos,centerPos,True)

    result = coarseCollection.filterDate(target.advance(-3,'month'),target.advance(3,'month'))\
                .map(apply_starfm,True)

    final = ee.Image(result.filterDate(target,target.advance(1,'day')).first())

    return final

def bathtub(wfrac,hand,permanent=None):
    '''
    Function to fit hi-resolution HAND model to water fraction estimate
    args:
        wFrac (ee.Image): water fraction image, values must be 0-1
        hand (ee.Image): height above nearest drainage (HAND) image, units in meters
        permanent (ee.Image): permanent water image to seed HAND filling
    '''
    def fillGrids(d):
        d = ee.Number(d)
        dMap = hand.lte(d).Or(permWater).reduceResolution(
         reducer=ee.Reducer.mean(),
         bestEffort=True,
         maxPixels = 1024).reproject(crs=proj)

        diff = wfrac.subtract(dMap).abs()
        return diff

    def minimizeDepth(d):
        d = ee.Number(d)
        dMap = hand.lte(d).Or(permWater)
        rd = dMap.reduceResolution(
          reducer = ee.Reducer.mean(),
          bestEffort = True,
          maxPixels = 1024
        ).reproject(crs=proj)

        diff = wfrac.subtract(rd).abs()

        out = nil.where(diff.eq(minDiff),dMap)
        err = nil.where(diff.eq(minDiff),diff)
        return out.rename('water').addBands(err.rename('error'))

    nil = ee.Image(0)
    if permanent:
        permWater = permanent
    else:
        permWater = ee.Image(0)
    proj = wfrac.projection()

    depths = ee.List.sequence(0,20)

    depthMaps = ee.ImageCollection(depths.map(fillGrids))
    minDiff = depthMaps.min()

    waterMap = ee.ImageCollection(depths.map(minimizeDepth))
    final = waterMap.select('water').max().addBands(waterMap.select('error').min())

    return final
