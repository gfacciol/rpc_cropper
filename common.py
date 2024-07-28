# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os
import sys
import logging
import datetime
import warnings
import subprocess
import numpy as np
import rasterio
from scipy import ndimage
from typing import Optional
import tempfile

logger = logging.getLogger()

# silent rasterio NotGeoreferencedWarning
warnings.filterwarnings("ignore",
                        category=rasterio.errors.NotGeoreferencedWarning)

# add the bin folder to system path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
bin_dir = os.path.join(parent_dir, 'bin')
os.environ['PATH'] = bin_dir + os.pathsep + os.environ['PATH']


def remove(target):
    try:
        os.remove(target)
    except OSError:
        pass

def tmpfile(ext=''):
    """
    Creates a temporary file in the /tmp directory.

    Args:
        ext: desired file extension

    Returns:
        absolute path to the created file

    The path of the created file is added to the garbage list to allow cleaning
    at the end of the pipeline.
    """
    fd, out = tempfile.mkstemp(suffix = ext, prefix = 's2p_', dir = '.')
    os.close(fd)           # http://www.logilab.org/blogentry/17873
    return out


def run(cmd, env=os.environ, timeout=None, shell=False):
    """
    Runs a shell command, and print it before running.

    Arguments:
        cmd: list of a command and its arguments, or as a fallback,
            a string to be passed to a shell that will be split into a list.
        env (optional, default value is os.environ): dictionary containing the
            environment variables
        timeout (optional, int): time in seconds after which the function will
            raise an error if the command hasn't returned

        TODO: remove the temporary `shell` argument once all commands use shell=False
        shell (bool): run the command in a subshell. Defaults to False.

    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """
    logging.info("RUN: %s", cmd)
    t = datetime.datetime.now()
    if not isinstance(cmd, list) and not shell:
        cmd = cmd.split()
    subprocess.run(cmd, shell=shell, stdout=sys.stdout, stderr=sys.stderr,
                   env=env, timeout=timeout, check=True)
    logging.info("execution time: %s", datetime.datetime.now() - t)


def matrix_translation(x, y):
    t = np.eye(3)
    t[0, 2] = x
    t[1, 2] = y
    return t


def rio_read_as_array_with_nans(im):
    """
    Read an image replacing gdal nodata value with np.nan

    Args:
        im: path to the input image file

    Returns:
        array: raster as numpy array
    """
    with rasterio.open(im, 'r') as src:
        array = src.read()
        nodata_values = src.nodatavals

    for band, nodata in zip(array, nodata_values):
        if nodata is not None:
            band[band == nodata] = np.nan

    return array.squeeze()



def image_size(im):
    """
    Reads the width and height of a geotiff

    Args:
        im: path to the input image file
    Returns:
        a tuple of size 2, giving width and height
    """
    try:
        with rasterio.open(im, 'r') as f:
            return (f.width, f.height)
    except :
        print("cannot state the size of the image %s" % str(im))
        sys.exit()


def rasterio_write(path, array, profile={}, tags={}):
    """
    Write a numpy array in a tiff or png file with rasterio.

    Args:
        path (str): path to the output tiff/png file
        array (numpy array): 2D or 3D array containing the image to write.
        profile (dict): rasterio profile (ie dictionary of metadata)
        tags (dict): dictionary with additional geotiff tags
    """
    # determine the driver based on the file extension
    extension = os.path.splitext(path)[1].lower()
    if extension in ['.tif', '.tiff']:
        driver = 'GTiff'
    elif extension in ['.png']:
        driver = 'png'
    else:
        raise NotImplementedError('format {} not supported'.format(extension))

    # read image size and number of bands
    array = np.atleast_3d(array)
    height, width, nbands = array.shape

    # define image metadata dict
    profile.update(driver=driver, count=nbands, width=width, height=height,
                   dtype=array.dtype)

    # write to file
    with rasterio.Env():
        with rasterio.open(path, 'w', **profile) as dst:
            dst.write(np.transpose(array, (2, 0, 1)))
            dst.update_tags(**tags)


def bounding_box2D(pts):
    """
    bounding box for the points pts
    """
    dim = len(pts[0])  # should be 2
    bb_min = [min([t[i] for t in pts]) for i in range(dim)]
    bb_max = [max([t[i] for t in pts]) for i in range(dim)]
    return bb_min[0], bb_min[1], bb_max[0] - bb_min[0], bb_max[1] - bb_min[1]


def maximum_filter_ignore_nan(array, *args, **kwargs):
    nans = np.isnan(array)
    replaced = np.where(nans, -np.inf, array)
    replaced = ndimage.maximum_filter(replaced, *args, **kwargs)
    return np.where(np.isinf(replaced), np.nan, replaced)

def minimum_filter_ignore_nan(array, *args, **kwargs):
    nans = np.isnan(array)
    replaced = np.where(nans, +np.inf, array)
    replaced = ndimage.minimum_filter(replaced, *args, **kwargs)
    return np.where(np.isinf(replaced), np.nan, replaced)

def cargarse_basura(inputf, outputf):
    se=5
    im = rio_read_as_array_with_nans(inputf)

    tmp1 = minimum_filter_ignore_nan(im,   size=se)
    tmp1 = maximum_filter_ignore_nan(tmp1, size=se)

    tmp2 = maximum_filter_ignore_nan(im,   size=se)
    tmp2 = minimum_filter_ignore_nan(tmp2, size=se)

    # put to nan if dilation minus erosion is larger than 5
    tmpM = np.where( np.abs(tmp1 - tmp2) > 5, np.nan, im)

    # remove small connected components
    rasterio_write(outputf, tmpM)
    run('remove_small_cc %s %s %d %d' % (outputf, outputf, 200, 5))
    #tmpM = specklefilter(tmpM,200,5)
    #rasterio_write(outputf, tmpM)


_t0 = datetime.datetime.now()
_t1: Optional[datetime.datetime] = None

def print_elapsed_time(since_first_call: bool = False) -> None:
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    global _t1
    t2 = datetime.datetime.now()
    if since_first_call:
        logging.info("Total elapsed time: %s", t2 - _t0)
    else:
        if _t1 is not None:
            logging.info("Elapsed time: %s", t2 - _t1)
        else:
            logging.info("Elapsed time: %s", t2 - _t0)
    _t1 = t2


def reset_elapsed_time() -> None:
    global _t1
    _t1 = datetime.datetime.now()


def linear_stretching_and_quantization_8bit(img, p=1):
    """
    Simple 8-bit quantization with linear stretching.

    Args:
        img (np.array): image to quantize
        p (float): percentage of the darkest and brightest pixels to saturate,
            from 0 to 100.

    Returns:
        numpy array with the quantized uint8 image
    """
    a, b = np.nanpercentile(img, (p, 100 - p))
    img = np.round(255 * (np.clip(img, a, b) - a) / (b - a))
    img = np.nan_to_num(img, nan=0)
    return img.astype(np.uint8)
