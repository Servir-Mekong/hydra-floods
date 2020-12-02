from __future__ import print_function, division
import datetime
import ee
import os
import sys
import subprocess
import gcsfs


def list_gcs_objs(bucket_path, pattern=None, output_url=False, project=None):
    """Function to list objects in Google Cloud Storage Bucket

    args:
        bucket_path (str): Google Cloud Storage bucket name
        pattern (str | None, optional): regex pattern to search in bucket. 
            Can seach folders by adding folder names (i.e. pattern = 'subfolder/*.txt).
            If None then will not use search pattern. default = None
        output_url (bool, optional): boolean switch to output google cloud storage http url
            or google cloud storage object uri. If false will output gcs uri. default = False
        project (str | None): Cloud project name to use when initiation file spec. If None then
            use default gcloud config. default = None

    returns:
        list[str]: List of objects in bucket that match pattern
    """
    fs = gcsfs.GCSFileSystem(project=project)
    if pattern is not None:
        bucket_path = (
            bucket_path + "/" if not bucket_path.endswith("/") else bucket_path
        )
        blobs = fs.glob(f"{bucket_path}{pattern}")
    else:
        blobs = fs.ls(bucket_path)

    base = "https://storage.cloud.google.com/{0}" if output_url else "gs://{0}"

    return [base.format(blob) for blob in blobs]


def push_to_gcs(file, bucket_path):
    """Helper function to copy local files to Google Cloud Storage
    Thinly wraps `gsutil cp` command line

    args:
        file (str): file path to push to GCS
        bucket_path (str): path on GCS to copy file to

    """
    if os.path.exists(file):
        cmd = "gsutil cp {0} {1}".format(file, bucket_path)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        out, err = proc.communicate()
    else:
        raise ValueError('file "{0} does not exist'.format(file))
    return


def push_to_ee(bucket_obj, asset_collection, properties=None, delete_bucket_obj=False):
    """Helper function to begin ingest process for imagery on GCS to GEE
    Thinly wraps `earthengine upload image`

    args:
        bucket_obj (str): GCS bucket object to ingest into GEE. Expects that object has mime type of image/tiff
        asset_collection (str): Earth Engine asset collection to push object to
        properties (list[str], optional): list of properties to set when ingesting files. If None then no properties
            will be set. default = None
        delete_bucket_obj (bool, optional): boolean switch to delete GCS object once ingested into EE. If set to False
            then file will remain on GCS. default = False
    """
    name = os.path.basename(bucket_obj).replace(".", "_")
    asset = asset_collection + name

    pStr = ""
    for i in properties:
        pStr += "--{0} {1} ".format(i, properties[i])

    binPath = os.path.dirname(sys.executable)
    cmd = "{0}/earthengine upload image --asset_id={1} {2} {3}".format(
        binPath, asset, pStr, bucket_obj
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
                        bucket_obj
                    )
                )
                sys.exit(1)

    if delete_bucket_obj:
        cmd = "gsutil rm {0}".format(bucket_obj)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        out, err = proc.communicate()

    return


def decode_date(date):
    """Decodes a date from a command line argument, returning msec since epoch".
    
    args:
        date (str): date value in a format that can be parsed into datetime object
        
    returns:
        datetime.datetime: decoded datetime value

    raises:
        TypeError: if string does not conform to a legal date format.
    """

    date_formats = [
        "%Y%m%d",
        "%Y-%m-%d",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%dT%H:%M:%S.%f",
    ]
    for date_format in date_formats:
        try:
            dt = datetime.datetime.strptime(date, date_format)
            return dt
        except ValueError:
            continue
    raise TypeError(f"Invalid value for property of type 'date': '{date}'.")

