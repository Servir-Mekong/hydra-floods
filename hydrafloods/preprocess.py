import os
import numpy as np
import xarray as xr
from scipy import ndimage, interpolate
from osgeo import gdal, osr
from pyproj import Proj, transform
from pyresample import bilinear, geometry, utils


def viirs(infile):
    tree = "//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/"
    field = "SurfReflect_{0}{1}_1"
    base = 'HDF5:"{0}":{1}{2}'

    outEpsg = 3857
    nd = -999

    # sinuWkt = '''PROJCS["Sphere_Sinusoidal",GEOGCS["GCS_Sphere",DATUM["D_Sphere",SPHEROID["Sphere",6371000.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]'''

    proj = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    outProj = Proj(proj)
    inProj = Proj(init="epsg:4326")

    srs = osr.SpatialReference()
    srs.ImportFromProj4(outProj.definition_string())

    m = [i for i in range(12) if i not in [0, 6, 9]]
    i = [i for i in range(1, 4)]
    bands = [m, i]
    flatBands = [item for sublist in bands for item in sublist]

    res = ["1km", "500m"]
    mode = ["M", "I"]

    band = gdal.Open(base.format(infile, tree.format("1km"), field.format("QF", 1)))
    metadata = band.GetMetadata()
    # print(metadata)
    cloudQA = _extractBits(band.ReadAsArray(), 2, 3)
    hiresCloudQA = ndimage.zoom(cloudQA, 2, order=0)
    band = None

    band = gdal.Open(base.format(infile, tree.format("1km"), field.format("QF", 2)))
    shadowQA = _extractBits(band.ReadAsArray(), 3, 3)
    hiresShadowQA = ndimage.zoom(shadowQA, 2, order=0)

    qa = ~(hiresCloudQA > 0) & (hiresShadowQA < 1)

    ringLatitude = metadata["GRingLatitude"].split(" ")[:-1]
    ringLongitude = metadata["GRingLongitude"].split(" ")[:-1]

    ll, ul, ur, lr = [
        transform(inProj, outProj, float(ringLongitude[i]), float(ringLatitude[i]))
        for i in range(len(ringLatitude))
    ]
    iniX, iniY = ul[0], ul[1]
    print(ll, ul, ur, lr)

    bandNames = [
        "{0}{1}".format(mode[i], bands[i][j])
        for i in range(len(res))
        for j in range(len(bands[i]))
    ]

    subdata = [
        [infile, res[i], mode[i], bands[i][j]]
        for i in range(len(res))
        for j in range(len(bands[i]))
    ]

    data = np.zeros([qa.shape[0], qa.shape[1], len(subdata) + 1])

    for i, s in enumerate(subdata):
        infile, r, m, b = s

        subdataset = base.format(infile, tree.format(r), field.format(m, b))

        band = gdal.Open(subdataset)
        if m == "M":
            data = ndimage.zoom(band.ReadAsArray(), 2, order=0)

        else:
            sd = np.array(band.ReadAsArray())
        band = None

        sd[np.where(sd < 0)] = nd

        data[:, :, i] = sd

    data[:, :, -1] = qa

    res = float(metadata["CharacteristicBinSize500M"])
    yskew = math.sin(_getSkew(ur, ul)) * res
    xskew = math.cos(_getSkew(ul, ll)) * res
    print(xskew, yskew)
    yDim, xDim = qa.shape
    gt = (10007554.677, res, 0, 2223901.039333, 0, -res)
    print(gt)

    name, _ = os.path.splitext(infile)
    out_name = name + ".TIF"

    write_geotiff(out_name, data, gt, outEpsg, no_data=nd)

    return out_name


def atms(
    infile, gridding_radius=25000,
):
    ds = xr.open_dataset(infile)

    outEpsg = 3857
    nd = -999

    outProj = Proj(init="epsg:{0}".format(outEpsg))
    inProj = Proj(init="epsg:4326")

    lons, lats = ds.lon.values, ds.lat.values

    xx, yy = transform(inProj, outProj, lons, lats)
    minx, miny = transform(inProj, outProj, -180, -86)
    maxx, maxy = transform(inProj, outProj, 180, 86)
    res = 16000

    eastings = np.arange(round(minx), round(maxx), res)
    northings = np.arange(round(miny), round(maxy), res)

    ee = eastings[np.where((eastings > xx.min()) & (eastings < xx.max()))]
    nn = northings[np.where((northings > yy.min()) & (northings < yy.max()))]

    swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
    area_def = geometry.AreaDefinition(
        "mercator",
        "WGS 84 / Pseudo-Mercator - Projected",
        "mercator",
        {
            "x_0": "0.0",
            "y_0": "0.0",
            "lat_ts": "0.00",
            "lon_0": "0.00",
            "proj": "merc",
            "k": "1.0",
            "datum": "WGS84",
            "ellps": "WGS84",
            "a": "6378137",
            "b": "6378137",
        },
        ee.size,
        nn.size,
        [ee.min(), nn.min(), ee.max(), nn.max()],
    )

    data = None

    # TODO: dynamically estimate sigama based on beam footprints
    eps = 0.1

    result = bilinear.resample_bilinear(
        ds.land_frac.where(ds["sat_zen"] < 50).values,
        swath_def,
        area_def,
        radius=gridding_radius,
        neighbours=32,
        nprocs=1,
        fill_value=nd,
        reduce_data=True,
        segments=None,
        epsilon=eps,
    )

    result[np.where(result >= 0)] = np.abs(result[np.where(result >= 0)] - 1) * 10000

    name, _ = os.path.splitext(infile)
    out_name = name + "_waterfrac.TIF"

    gt = (ee.min(), res, 0, nn.max(), 0, -res)

    write_geotiff(out_name, result, gt, outEpsg, no_data=nd)

    return out_name


def write_geotiff(out_name, data, gt, epsg, no_data=None):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)

    if len(data.shape) == 3:
        yDim, xDim = data.shape[:2]
        nBands = data.shape[2]

    elif len(data.shape) == 2:
        yDim, xDim = data.shape
        nBands = 1
    else:
        raise

    driver = gdal.GetDriverByName("GTiff")

    outDs = driver.Create(out_name, xDim, yDim, nBands, gdal.GDT_Int16)
    outDs.SetGeoTransform(gt)
    # set to something that is not user defined
    outDs.SetProjection(srs.ExportToWkt())

    for b in range(nBands):
        band = outDs.GetRasterBand(b + 1)
        if no_data:
            band.Setno_dataValue(-999)
        if nBands > 1:
            band.WriteArray(data[:, :, b])
        else:
            band.WriteArray(data)

        band = None

    outDs.FlushCache()

    return
