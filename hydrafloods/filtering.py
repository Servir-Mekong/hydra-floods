import ee

def powerToDb(img):
  return ee.Image(10).multiply(img.log10())

def dbToPower(img):
  return ee.Image(10).pow(img.divide(10))


def leeSigma(img):
    meanReducer = ee.Reducer.mean(),
    varianceReducer = ee.Reducer.variance(),
    sumReducer = ee.Reducer.sum()
    bandNames = img.bandNames()

    img = dbToPower(img)

    kernelWeights = ee.List.repeat(ee.List.repeat(1,9),9)
    kernel = ee.Kernel.fixed(9,9,kernelWeights,5,5)

    targetWeights = ee.List.repeat(ee.List.repeat(1,3),3)
    targetKernel = ee.Kernel.fixed(3,3,targetWeights,1,1)

    Tk = ee.Image(7)
    sigma = '0.9'
    looks = '4'

    # Lookup table for range and eta values for intensity
    sigmaLookup = ee.Dictionary({
    1: ee.Dictionary({
      0.5: ee.Dictionary({
        'A1': 0.436,
        'A2': 1.92,
        'η': 0.4057
      }),
      0.6: ee.Dictionary({
        'A1': 0.343,
        'A2': 2.21,
        'η': 0.4954
      }),
      0.7: ee.Dictionary({
        'A1': 0.254,
        'A2': 2.582,
        'η': 0.5911
      }),
      0.8: ee.Dictionary({
        'A1': 0.168,
        'A2': 3.094,
        'η': 0.6966
      }),
      0.9: ee.Dictionary({
        'A1': 0.084,
        'A2': 3.941,
        'η': 0.8191
      }),
      0.95: ee.Dictionary({
        'A1': 0.043,
        'A2': 4.840,
        'η': 0.8599
      })
    }),
    2: ee.Dictionary({
      0.5: ee.Dictionary({
        'A1': 0.582,
        'A2': 1.584,
        'η': 0.2763
      }),
      0.6: ee.Dictionary({
        'A1': 0.501,
        'A2': 1.755,
        'η': 0.3388
      }),
      0.7: ee.Dictionary({
        'A1': 0.418,
        'A2': 1.972,
        'η': 0.4062
      }),
      0.8: ee.Dictionary({
        'A1': 0.327,
        'A2': 2.260,
        'η': 0.4819
      }),
      0.9: ee.Dictionary({
        'A1': 0.221,
        'A2': 2.744,
        'η': 0.5699
      }),
      0.95: ee.Dictionary({
        'A1': 0.152,
        'A2': 3.206,
        'η': 0.6254
      }),
    }),
    3: ee.Dictionary({
      0.5: ee.Dictionary({
        'A1': 0.652,
        'A2': 1.458,
        'η': 0.2222
      }),
      0.6: ee.Dictionary({
        'A1': 0.580,
        'A2': 1.586,
        'η': 0.2736
      }),
      0.7: ee.Dictionary({
        'A1': 0.505,
        'A2': 1.751,
        'η': 0.3280
      }),
      0.8: ee.Dictionary({
        'A1': 0.419,
        'A2': 1.865,
        'η': 0.3892
      }),
      0.9: ee.Dictionary({
        'A1': 0.313,
        'A2': 2.320,
        'η': 0.4624
      }),
      0.95: ee.Dictionary({
        'A1': 0.238,
        'A2': 2.656,
        'η': 0.5084
      }),
    }),
    4: ee.Dictionary({
      0.5: ee.Dictionary({
        'A1': 0.694,
        'A2': 1.385,
        'η': 0.1921
      }),
      0.6: ee.Dictionary({
        'A1': 0.630,
        'A2': 1.495,
        'η': 0.2348
      }),
      0.7: ee.Dictionary({
        'A1': 0.560,
        'A2': 1.627,
        'η': 0.2825
      }),
      0.8: ee.Dictionary({
        'A1': 0.480,
        'A2': 1.804,
        'η': 0.3354
      }),
      0.9: ee.Dictionary({
        'A1': 0.378,
        'A2': 2.094,
        'η': 0.3991
      }),
      0.95: ee.Dictionary({
        'A1': 0.302,
        'A2': 2.360,
        'η': 0.4391
      }),
    })
    })

    # extract data from lookup
    looksDict = ee.Dictionary(sigmaLookup.get(ee.String(looks)))
    sigmaImage = ee.Dictionary(looksDict.get(ee.String(sigma))).toImage()
    a1 = sigmaImage.select('A1')
    a2 = sigmaImage.select('A2')
    aRange = a2.subtract(a1)
    eta = sigmaImage.select('η').pow(2)

    # MMSE estimator
    mmseMask = img.gte(a1).or(img.lte(a2))
    mmseIn = img.updateMask(mmseMask)
    oneImg = ee.Image(1)
    z = mmseIn.reduceNeighborhood(meanReducer,kernel,null,true)
    varz = mmseIn.reduceNeighborhood(varianceReducer,kernel)
    varx = (varz.subtract(z.abs().pow(2).multiply(eta))).divide(oneImg.add(eta))
    b = varx.divide(varz)
    mmse = oneImg.subtract(b).multiply(z.abs()).add(b.multiply(mmseIn))

    # workflow
    z99 = ee.Dictionary(img.reduceRegion({
    reducer: ee.Reducer.percentile([99],null,255,0.001,1e6),
    geometry: img.geometry(),
    scale:10,
    bestEffort:true
    })).toImage()

    overThresh = img.gte(z99)

    K = overThresh.reduceNeighborhood(sumReducer,targetKernel,null,true)

    retainPixel = K.gte(Tk)
    xHat = powerToDb(img.updateMask(retainPixel).unmask(mmse))

    return ee.Image(xHat).rename(bandNames).copyProperties(img)


# The RL speckle filter
def refinedLee(image):
    def applyFilter(b):
        img = image.select([b])

        # img must be in natural units, i.e. not in dB!
        # Set up 3x3 kernels
        weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
        kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false)

        mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
        variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

        # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
        sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

        sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false)

        # Calculate mean and variance for the sampled windows and store as 9 bands
        sample_mean = mean3.neighborhoodToBands(sample_kernel)
        sample_= variance3.neighborhoodToBands(sample_kernel)

        # Determine the 4 gradients for the sampled windows
        gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
        gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
        gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
        gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

        # And find the maximum gradient amongst gradient bands
        max_gradient = gradients.reduce(ee.Reducer.max())

        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(max_gradient)

        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask)

        # Determine the 8 directions
        directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
        directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
        directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
        directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).not().multiply(5))
        directions = directions.addBands(directions.select(1).not().multiply(6))
        directions = directions.addBands(directions.select(2).not().multiply(7))
        directions = directions.addBands(directions.select(3).not().multiply(8))

        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask)

        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum())

        #pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000']
        #Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false)

        sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

        # Calculate localNoiseVariance
        sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

        # Set up the 7*7 kernels for directional statistics
        rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

        diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],
          [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]])

        rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false)
        diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false)

        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
        dir_= img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
        dir_= dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

        # and add the bands for rotated kernels
        for i in range(1,4):
          dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
          dir_= dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
          dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
          dir_= dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))


        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
        dir_mean = dir_mean.reduce(ee.Reducer.sum())
        dir_= dir_var.reduce(ee.Reducer.sum())

        # A finally generate the filtered value
        varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))

        b = varX.divide(dir_var)

        # return multi-band image band from array
        return dir_mean.add(b.multiply(img.subtract(dir_mean)))\
          .arrayProject([0])\
          .arrayFlatten([['sum']])\
          .float()

    bandNames = image.bandNames()
    image = dbToPower(image)

    result = ee.ImageCollection(bandNames.map(applyFilter)).toBands().rename(bandNames)
    return powerToDb(ee.Image(result))


def gammaMap(img):

    ksize = 7
    enl = 5
    bandNames = img.bandNames()

    # Convert image from dB to natural values
    nat_img = dbToPower(img)

    # Square kernel, ksize should be odd (typically 3, 5 or 7)
    weights = ee.List.repeat(ee.List.repeat(1,ksize),ksize)

    # ~~(ksize/2) does integer division in JavaScript
    kernel = ee.Kernel.fixed(ksize,ksize, weights, ~~(ksize/2), ~~(ksize/2), false)

    # Get mean and variance
    mean = nat_img.reduceNeighborhood(ee.Reducer.mean(), kernel)
    variance = nat_img.reduceNeighborhood(ee.Reducer.variance(), kernel)

    # "Pure speckle" threshold
    ci = variance.sqrt().divide(mean)  # square root of inverse of enl

    # If ci <= cu, the kernel lies in a "pure speckle" area -> return simple mean
    cu = 1.0/Math.sqrt(enl)

    # If cu < ci < cmax the kernel lies in the low textured speckle area -> return the filtered value
    cmax = Math.sqrt(2.0) * cu

    alpha = ee.Image(1.0 + cu*cu).divide(ci.multiply(ci).subtract(cu*cu))
    b = alpha.subtract(enl + 1.0)
    d = mean.multiply(mean).multiply(b).multiply(b).add(alpha.multiply(mean).multiply(nat_img).multiply(4.0*enl))
    f = b.multiply(mean).add(d.sqrt()).divide(alpha.multiply(2.0))

    caster = ee.Dictionary.fromLists(bandNames,ee.List.repeat('float',3))
    img1 = powerToDb(mean.updateMask(ci.lte(cu))).rename(bandNames).cast(caster)
    img2 = powerToDb(f.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))).rename(bandNames).cast(caster)
    img3 = img.updateMask(ci.gte(cmax)).rename(bandNames).cast(caster)

    # If ci > cmax do not filter at all (i.e. we don't do anything, other then masking)
    result = ee.ImageCollection([img1,img2,img3])
    .reduce(ee.Reducer.firstNonNull()).rename(bandNames)

    # Compose a 3 band image with the mean filtered "pure speckle", the "low textured" filtered and the unfiltered portions
    return result
