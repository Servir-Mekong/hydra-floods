from __future__ import print_function,division
import numpy as np
from scipy import ndimage
from PIL import Image, ImageDraw


def find_nearest(xx,yy,xval,yval):
    xidx = (np.abs(xx-xval)).argmin()
    yidx = (np.abs(yy-yval)).argmin()

    ridx = yidx / xx.shape[1]
    cidx = xidx % xx.shape[1]

    return [ridx,cidx]

def geoverts_2_imgverts(poly,yy,xx):
    """
    FUNCTION:  verts = geoverts_2_imgverts(xx,yy,poly)
    ARGUMENTS: xx - x-coordinate grid for output raster
               yy - y-coordinate grid for output raster
               poly - swath polygon vertices
    KEYWORDS:  N/A
    RETURNS:   vertsOut - the polygon vertices calculated from the swath
    NOTES:     Takes polygon geographic coordinates and converts them into
               image pixel coordinates for creating a mask image.
    """

    vertsOut = []

    # Iterate through each polygon vertex
    for i in range(len(poly)):
        # Extract x-y coordinate from vertex
        xval = poly[i][0]
        yval = poly[i][1]

        # Find closest image index to the x-y coordinate
        xidx = (np.abs(xx-xval)).argmin()
        yidx = (np.abs(yy-yval)).argmin()

        # Convert the 1-d index to 2-d
        ridx = yidx / xx.shape[1]
        cidx = xidx % xx.shape[1]

        vertsOut.append((int(cidx),int(ridx)))

    return vertsOut

def hist_match(source, template):
    """
    Adjust the pixel values of a grayscale image such that its histogram
    matches that of a target image

    Arguments:
    -----------
        source: np.ndarray
            Image to transform; the histogram is computed over the flattened
            array
        template: np.ndarray
            Template image; can have different dimensions to source
    Returns:
    -----------
        matched: np.ndarray
            The transformed output image
    """

    oldshape = source.shape
    source = source.ravel()
    template = template.ravel()

    # get the set of unique pixel values and their corresponding indices and
    # counts
    s_values, bin_idx, s_counts = np.unique(source, return_inverse=True,
                                            return_counts=True)
    t_values, t_counts = np.unique(template, return_counts=True)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

    return interp_t_values[bin_idx].reshape(oldshape)


def clip(toBeClipped, toClipTo, targetCoordslonKey='ILon',latKey='ILat'):
    if len(targetCoords) != 2:
        raise Error('Wrong input for targetCoords argument, expected list of [lat, lon] keys')

    templateLon = toBeClipped[target]
    templateLat = toBeClipped[latKey]

    mask = np.zeros_like(templateLon)

    minx,maxx = toClipTo['Lon'].min(),toClipTo['Lon'].max()
    miny,maxy = toClipTo['Lat'].min(),toClipTo['Lat'].max()

    mask = np.ma.masked_where(templateLon < minx, mask)
    mask = np.ma.masked_where(templateLon > maxx, mask)
    mask = np.ma.masked_where(templateLat < miny, mask)
    mask = np.ma.masked_where(templateLat > maxy, mask)

    return

def mask(linearRing,yCoords,xCoords):

    verts = geoverts_2_imgverts(linearRing,yCoords,xCoords)

    # Create a mask from the convex hull vertices
    img = Image.new('L', (yCoords.shape[1], yCoords.shape[0]), 0)
    ImageDraw.Draw(img).polygon(verts, outline=1, fill=1)
    mask = np.array(img)

    return mask
