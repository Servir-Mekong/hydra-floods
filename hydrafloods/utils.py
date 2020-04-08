from __future__ import print_function,division
import datetime
import ee
import os
import sys
import subprocess
import numpy as np



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


def push_to_gcs(file,bucketPath):
    if os.path.exists(file):
        cmd = "gsutil cp {0} {1}".format(file,bucketPath)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = proc.communicate()
    else:
        raise ValueError('file "{0} does not exist'.format(file))
    return

def push_to_gee(bucketObj,assetCollection,properties=None,deleteBucketObj=True):
    name = os.path.basename(bucketObj).replace('.','_')
    asset = assetCollection + name

    pStr = ''
    for i in properties:
         pStr += '--{0} {1} '.format(i,properties[i])

    binPath = os.path.dirname(sys.executable)
    cmd = "{0}/earthengine upload image --asset_id={1} {2} {3}".format(binPath,asset,pStr,bucketObj)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = proc.communicate()
    if properties:
        pStr = ''

        running = True
        while running == True:
            tasks = ee.batch.Task.list()
            if 'COMPLETED' in str(tasks[0]):
                running = False
            elif 'FAILED' in str(tasks[0]):
                print('EE upload process failed for image {}, check Earth Engine for error'.format(bucketObj))
                sys.exit(1)

    if deleteBucketObj:
        cmd = "gsutil rm {0}".format(bucketObj)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = proc.communicate()

    return


def decode_date(string):
  """Decodes a date from a command line argument, returning msec since epoch".
  Args:
    string: See AssetSetCommand class comment for the allowable
      date formats.
  Returns:
    long, ms since epoch
  Raises:
    argparse.ArgumentTypeError: if string does not conform to a legal
      date format.
  """

  try:
    return int(string)
  except ValueError:
    date_formats = ['%Y%m%d',
                    '%Y-%m-%d',
                    '%Y-%m-%dT%H:%M:%S',
                    '%Y-%m-%dT%H:%M:%S.%f']
    for date_format in date_formats:
      try:
        dt = datetime.datetime.strptime(string, date_format)
        return dt
      except ValueError:
        continue
  raise argparse.ArgumentTypeError(
      'Invalid value for property of type "date": "%s".' % string)

def parse_atms_time(infile):

    return

def parse_viirs_time(infile):

    return
