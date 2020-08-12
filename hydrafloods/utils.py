from __future__ import print_function, division
import datetime
import ee
import os
import sys
import subprocess
import numpy as np
import gcsfs


def list_gcs_objs(bucket_path, pattern=None, output_url=False):
    fs = gcsfs.GCSFileSystem()
    if pattern is not None:
        bucket_path = (
            bucket_path + "/" if not bucket_path.endswith("/") else bucket_path
        )
        blobs = fs.glob(f"{bucket_path}{pattern}")
    else:
        blobs = fs.ls(bucket_path)

    base = "https://storage.cloud.google.com/{0}" if output_url else "gs://{0}"

    return [base.format(blob) for blob in blobs]


def push_to_gcs(file, bucketPath):
    if os.path.exists(file):
        cmd = "gsutil cp {0} {1}".format(file, bucketPath)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        out, err = proc.communicate()
    else:
        raise ValueError('file "{0} does not exist'.format(file))
    return


def push_to_gee(bucket_obj, asset_collection, properties=None, delete_bucket_obj=False):
    name = os.path.basename(bucketObj).replace(".", "_")
    asset = assetCollection + name

    pStr = ""
    for i in properties:
        pStr += "--{0} {1} ".format(i, properties[i])

    binPath = os.path.dirname(sys.executable)
    cmd = "{0}/earthengine upload image --asset_id={1} {2} {3}".format(
        binPath, asset, pStr, bucketObj
    )
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    out, err = proc.communicate()
    if properties:
        pStr = ""

        running = True
        while running == True:
            tasks = ee.batch.Task.list()
            if "COMPLETED" in str(tasks[0]):
                running = False
            elif "FAILED" in str(tasks[0]):
                print(
                    "EE upload process failed for image {}, check Earth Engine for error".format(
                        bucketObj
                    )
                )
                sys.exit(1)

    if deleteBucketObj:
        cmd = "gsutil rm {0}".format(bucketObj)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
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
        date_formats = [
            "%Y%m%d",
            "%Y-%m-%d",
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%dT%H:%M:%S.%f",
        ]
        for date_format in date_formats:
            try:
                dt = datetime.datetime.strptime(string, date_format)
                return dt
            except ValueError:
                continue
    raise argparse.ArgumentTypeError(
        'Invalid value for property of type "date": "%s".' % string
    )

