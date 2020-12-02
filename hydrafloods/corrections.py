import ee
import math
from hydrafloods import geeutils, decorators


@decorators.carry_metadata
def slope_correction(image, elevation, model="volume", buffer=0, scale=1000):
    """This function applies the slope correction on a Sentinel-1 image.
    Function based on https:# doi.org/10.3390/rs12111867.
    Adapted from https:# github.com/ESA-PhiLab/radiometric-slope-correction/blob/master/notebooks/1%20-%20Generate%20Data.ipynb
       
    args:
        image (ee.Image): Sentinel-1 to perform correction on
        elevation (ee.Image): Input DEM to calculate slope corrections from
        model (str, optional): physical reference model to be applied. Options are 'volume' or 'surface'.
            default = volume
        buffer (int, optional): buffer in meters for layover/shadow mask. If zero then no buffer will be applied. default = 0
        scale (int, optional): reduction scale to process satellite heading compared to ground. Increasing will reduce
            chance of OOM errors but reduce local scale correction accuracy. default = 1000
        
    returns:
        ee.Image: slope corrected SAR imagery with look and local incidence angle bands

    raises:
        NotImplementedError: when keyword model is not of 'volume' or 'surface'
    """

    def _volumetric_model_SCF(theta_iRad, alpha_rRad):
        """Closure funnction for calculation of volumetric model SCF
        
        args:
            theta_iRad (ee.Image): incidence angle in radians
            alpha_rRad (ee.Image): slope steepness in range
        
        returns:
            ee.Image
        """

        # model
        nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan()
        denominator = (ninetyRad.subtract(theta_iRad)).tan()
        return nominator.divide(denominator)

    def _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad):
        """Closure funnction for calculation of direct model SCF
        
        args:
            theta_iRad (ee.Image): incidence angle in radians
            alpha_rRad (ee.Image): slope steepness in range
            alpha_azRad (ee.Image): slope steepness in azimuth
        
        returns:
            ee.Image
        """

        # model
        nominator = (ninetyRad.subtract(theta_iRad)).cos()
        denominator = alpha_azRad.cos().multiply(
            (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos()
        )

        return nominator.divide(denominator)

    def _erode(image, distance):
        """Closure function to buffer raster values

        args:
            image (ee.Image): image that should be buffered
            distance (int): distance of buffer in meters
        
        returns: 
            ee.Image
      """

        d = (
            image.Not()
            .unmask(1)
            .fastDistanceTransform(10)
            .sqrt()
            .multiply(ee.Image.pixelArea().sqrt())
        )

        return image.updateMask(d.gt(distance))

    def _masking(alpha_rRad, theta_iRad, buffer):
        """Closure function for masking of layover and shadow
        
        args:
            alpha_rRad (ee.Image): slope steepness in range
            theta_iRad (ee.Image): incidence angle in radians
            buffer (int): buffer in meters
        
        returns: 
            ee.Image
        """
        # layover, where slope > radar viewing angle
        layover = alpha_rRad.lt(theta_iRad).rename("layover")

        # shadow
        shadow = alpha_rRad.gt(
            ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))
        ).rename("shadow")

        # add buffer to layover and shadow
        if buffer > 0:
            layover = _erode(layover, buffer)
            shadow = _erode(shadow, buffer)

        # combine layover and shadow
        no_data_mask = layover.And(shadow).rename("no_data_mask")

        return no_data_mask

    # get the image geometry and projection
    geom = image.geometry(scale)
    proj = image.select(1).projection()

    # image to convert angle to radians
    to_radians = ee.Image.constant((math.pi / 180))
    # create a 90 degree image in radians
    ninetyRad = ee.Image.constant(90).multiply(to_radians)

    # calculate the look direction
    heading = (
        ee.Terrain.aspect(image.select("angle"))
        .reduceRegion(ee.Reducer.mean(), geom, scale)
        .get("aspect")
    )

    # Sigma0 to Power of input image
    sigma0Pow = geeutils.db_to_power(image)

    # the numbering follows the article chapters
    # 2.1.1 Radar geometry
    theta_iRad = image.select("angle").multiply(to_radians)
    phi_iRad = ee.Image.constant(heading).multiply(to_radians)

    # 2.1.2 Terrain geometry
    alpha_sRad = (
        ee.Terrain.slope(elevation)
        .select("slope")
        .multiply(to_radians)
        .setDefaultProjection(proj)
    )

    phi_sRad = (
        ee.Terrain.aspect(elevation)
        .select("aspect")
        .multiply(to_radians)
        .setDefaultProjection(proj)
    )

    # 2.1.3 Model geometry
    # reduce to 3 angle
    phi_rRad = phi_iRad.subtract(phi_sRad)

    # slope steepness in range (eq. 2)
    alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

    # slope steepness in azimuth (eq 3)
    alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

    # local incidence angle (eq. 4)
    theta_liaRad = (
        alpha_azRad.cos().multiply((theta_iRad.subtract(alpha_rRad)).cos())
    ).acos()
    theta_liaDeg = theta_liaRad.multiply(180 / math.pi)

    # 2.2
    # Gamma_nought
    gamma0 = sigma0Pow.divide(theta_iRad.cos())
    gamma0dB = geeutils.db_to_power(gamma0).select(
        ["VV", "VH"], ["VV_gamma0", "VH_gamma0"]
    )

    if model == "volume":
        scf = _volumetric_model_SCF(theta_iRad, alpha_rRad)

    elif model == "surface":
        scf = _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)

    else:
        raise NotImplementedError(
            f"Defined model, {model}, has not been implemented. Options are 'volume' or 'surface'"
        )

    # apply model for Gamm0_f
    gamma0_flat = gamma0.divide(scf)
    gamma0_flatDB = geeutils.power_to_db(gamma0_flat).select(["VV", "VH"])

    # calculate layover and shadow mask
    masks = _masking(alpha_rRad, theta_iRad, buffer)

    return (
        gamma0_flatDB.updateMask(masks)
        .addBands(image.select("angle"))
        .addBands(theta_liaDeg.rename("local_inc_angle"))
    )


@decorators.carry_metadata
def illumination_correction(image, elevation, model="rotation", scale=90, sensor="LC8"):
    """This function applies a terrain correction to optical imagery based on solar and viewing geometry
     
    args:
        image (ee.Image): Optical image to perform correction on
        elevation (ee.Image): Input DEM to calculate illumination corrections from
        model (str, optional): correction model to be applied. Options are 'cosine', 'c', 'scsc', or 'rotation'
            default = rotation
        scale (int, optional): reduction scale to process satellite heading compared to ground. Increasing will reduce
            chance of OOM errors but reduce local scale correction accuracy. default = 90
        sensor (str, optional): name of sensor to correct. options are 'LC8' or 'S2' (lower case also accepted).
            default = LC8
        
    returns:
        ee.Image: illumination corrected optical imagery

    raises:
        NotImplementedError: when keyword sensor is not of 'LC8' or 'S2'
        NotImplementedError: when keyword model is not of 'cosine', 'c', 'scsc', or 'rotation'
    """

    def _get_band_coeffs(band_name):
        """Closure function to find illumination correction fit across the different bands

        args:
            band_name (str | ee.String): band name to find correction coefficients for
        """

        #  Create the image to apply the linear regression.The first band
        # is the cosi and the second band is the response variable, the reflectance (the bands).
        # L (y) = a + b*cosi(x); a = intercept, b = slope
        #  Dependent: Reflectance
        y = image.select([band_name])
        #  create an image with the three variables by concatenating them
        reg_image = ee.Image.cat([cosi, one, y])
        #  specify the linear regression reducer
        lr_reducer = ee.Reducer.linearRegression(numX=2, numY=1)
        #  fit the model
        fit = reg_image.reduceRegion(
            reducer=lr_reducer, geometry=image.geometry(), scale=scale, maxPixels=1e10
        )

        #  Get the coefficients as a nested list, cast it to an array, and get
        #  just the selected column
        slope = ee.Array(fit.get("coefficients")).get([0, 0])
        intercept = ee.Array(fit.get("coefficients")).get([1, 0])

        return ee.List([slope, intercept])

    if sensor.lower() == "lc8":
        sz_property = "SOLAR_ZENITH_ANGLE"
        sa_property = "SOLAR_AZIMUTH_ANGLE"
    elif sensor.lower() == "s2":
        sz_property = "MEAN_SOLAR_ZENITH_ANGLE"
        sa_property = "MEAN_SOLAR_AZIMUTH_ANGLE"
    else:
        raise NotImplementedError(
            f"Selected sensor, {sensor}, is not available. Options are 'LC8' or 'S2' (lower case also accepted)"
        )

    # value convert angle to radians
    to_radians = ee.Number((math.pi / 180))
    # constand image of 1
    one = ee.Image.constant(1).rename("one")

    # calculate terrain info from elevation data
    terrain = ee.Algorithms.Terrain(elevation)
    # Extract slope in radians for each pixel in the image
    p = terrain.select(["slope"]).multiply(to_radians)
    # Extract aspect in radians for each pixel in the image
    o = terrain.select(["aspect"]).multiply(to_radians)
    # Extract solar zenith angle from the image
    z = ee.Image.constant(ee.Number(image.get(sz_property)).multiply(to_radians))
    # Extract solar azimuth from the image
    az = ee.Image.constant(ee.Number(image.get(sa_property)).multiply(to_radians))

    cosao = (o.subtract(az)).cos()
    # cos(ϕa−ϕo)
    # Calculate the cosine of the local solar incidence for every pixel in the image in radians (cosi=cosp*cosz+sinp*sinz*cos(ϕa−ϕo)
    cosi = image.expression(
        "((cosp * cosz) + (sinp * sinz * cosao))",
        {
            "cosp": p.cos(),
            "cosz": z.cos(),
            "sinp": p.sin(),
            "sinz": z.sin(),
            "cosao": cosao,
        },
    )

    if model == "cosine":
        # if cosine model correction, return early as we don't need to do extra processing
        return image.expression(
            "((image * cosz) / cosi) ", {"image": image, "cosz": z.cos(), "cosi": cosi}
        )

    bnames = image.bandNames()
    ab = ee.Array(bnames.map(_get_band_coeffs))

    # get the coefficients as images
    a = ee.Image(ee.Array(ab.slice(1, 0, 1))).arrayProject([0]).arrayFlatten([bnames])
    b = ee.Image(ee.Array(ab.slice(1, 1, 2))).arrayProject([0]).arrayFlatten([bnames])
    C = b.divide(a)

    if model == "c":
        newimage = image.expression(
            "((image * (cosz + C)) / (cosi + C))",
            {"image": image, "cosz": z.cos(), "cosi": cosi, "C": C},
        )

    elif model == "scsc":
        newimage = image.expression(
            "((image * ((cosp * cosz) + C))/(cosi + C))",
            {"image": image, "cosp": p.cos(), "cosz": z.cos(), "cosi": cosi, "C": C},
        )

    elif model == "rotation":
        # Apply the empirical rotation model
        newimage = image.expression(
            "image - a * (cosi - cosz)",
            {"image": image, "cosz": z.cos(), "cosi": cosi, "a": a},
        )

    else:
        raise NotImplementedError(
            f"Defined model, {model}, has not been implemented. Options are 'cosine', 'c', 'scsc', or 'rotation'"
        )

    return newimage
