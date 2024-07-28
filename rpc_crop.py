#!/usr/bin/env python3

# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>
# Copyright (C) 2013, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>

import os
import sys

import rpcm
from rpcm.rpc_file_readers import read_rpc_file
import rpc_utils
import srtm4 as srtm
import common
import rasterio



def procedure1(poly, a11, a22, b1, b2):
   """
   given the polynomial of degree three of three variables poly(x,y,z)
   computes the polynomial coefficients implementing
   the variable change x = a11 x' + b1
                       y = a22 y' + b2
   here the variable z is height of the RPC model
   VERIFIED!
   """
   newpoly = [None] * 20
   newpoly[0]    =  poly[0] + b2*b2*b2* poly[11] + b1*b1* b2* poly[12] + b1* b2*b2* poly[14] + b1*b1*b1* poly[15] + b2* poly[1] + b1* poly[2] + b1* b2* poly[4] + b2*b2* poly[7] + b1*b1* poly[8]
   newpoly[1]    =  (3* a22* b2*b2* poly[11] + a22* b1*b1 *poly[12] + 2* a22* b1* b2* poly[14] + a22* poly[1] + a22* b1* poly[4] + 2* a22* b2* poly[7])
   newpoly[2]    =  (2* a11* b1* b2* poly[12] + a11* b2*b2* poly[14] + 3* a11* b1*b1* poly[15] + a11* poly[2] + a11* b2* poly[4] + 2* a11* b1* poly[8])
   newpoly[3]    =  (b1* b2* poly[10] + b2*b2* poly[17] + b1*b1* poly[18] + poly[3] + b2* poly[5] + b1* poly[6])
   newpoly[4]    =  (2* a11* a22* b1* poly[12] + 2* a11* a22* b2* poly[14] + a11* a22* poly[4])
   newpoly[5]    =  (a22* b1* poly[10] + 2* a22* b2* poly[17] + a22* poly[5])
   newpoly[6]    =  (a11* b2* poly[10] + 2* a11* b1* poly[18] + a11* poly[6])
   newpoly[7]    =  (3* a22*a22* b2* poly[11] + a22*a22* b1* poly[14] + a22*a22* poly[7])
   newpoly[8]    =  (a11*a11* b2* poly[12] + 3 *a11*a11* b1* poly[15] + a11*a11* poly[8])
   newpoly[9]    =  (b2* poly[13] + b1* poly[16] + poly[9])
   newpoly[10]   =  a11* a22* poly[10]
   newpoly[11]   =  a22*a22*a22* poly[11]
   newpoly[12]   =  a11*a11* a22* poly[12]
   newpoly[13]   =  a22* poly[13]
   newpoly[14]   =  a11* a22*a22* poly[14]
   newpoly[15]   =  a11*a11*a11* poly[15]
   newpoly[16]   =  a11* poly[16]
   newpoly[17]   =  a22*a22* poly[17]
   newpoly[18]   =  a11*a11* poly[18]
   newpoly[19]   =  poly[19]
   return newpoly


def poly_variable_change_in(polyNum, polyDen, a11, a22, b1, b2):
   """
   given the RPC polynomials polyNum(x,y,z)/polyDen(x,y,z)
   computes the polynomial coefficients implementing
   the variable change x = a11 x' + b1
                       y = a22 y' + b2
   VERIFIED!
   """
   #print a11,a22,b1,b2
   newNum = procedure1(polyNum, a11, a22, b1, b2)
   newDen = procedure1(polyDen, a11, a22, b1, b2)
   return newNum, newDen


def poly_variable_change_out(polyNum, polyDen, a11, b1):
   """
   given the RPC polynomials polyNum(x,y,z)/polyDen(x,y,z)
   computes the polynomial coefficients implementing
   the operation   a11*(polyNum(x,y,z)/polyDen(x,y,z)) + b1
   VERIFIED!
   """
   import numpy as np
   newNum = list(float(a11) * np.array(polyNum) + float(b1) * np.array(polyDen))
   newDen = polyDen
   return newNum, newDen



def rpc_apply_crop_to_rpc_model(rpc, x0, y0, w, h):
   import copy
   rpcout = copy.deepcopy(rpc)

   rpcout.col_offset -= x0
   rpcout.row_offset -= y0

#   ## compute the scale and shift parameter for the normalized RPC
#   a11 = (float(w)/2) / (rpc.col_scale)
#   a22 = (float(h)/2) / (rpc.row_scale)
##   a11 = 1.0
##   a22 = 1.0
#   b1  = float(x0)/float(rpc.col_scale)
#   b2  = float(y0)/float(rpc.row_scale)
#   ## apply the transformation to the direct polynomials, BEWARE!! The order of the variables is reversed
#   rpcout.directLonNum, rpcout.directLonDen = poly_variable_change_in(rpc.lon_num, rpc.lon_den, a22,a11,b2,b1)
#   rpcout.directLatNum, rpcout.directLatDen = poly_variable_change_in(rpc.lat_num, rpc.lat_den, a22,a11,b2,b1)
#
#
#   # scale the RPC domain (I'm not sure its [-1,1]^2)
#   #   # TODO correct RPC so that the validity domain is still the square [-1,1]^2
#   rpcout.col_scale= float(w)/2
#   rpcout.row_scale= float(h)/2
##   rpcout.col_scale= rpc.col_scale  ## keep it unchanged (it also works)
##   rpcout.row_scale= rpc.row_scale
#

#   ## compute the scale and shift parameters for the normalized RPC
#   b1 = float(x0)/float(rpcout.col_scale)
#   b2 = float(y0)/float(rpcout.row_scale)
#   a11 = float(rpc.col_scale)/float(rpcout.col_scale)
#   a22 = float(rpc.row_scale)/float(rpcout.row_scale)
##   a11 = 1.0
##   a22 = 1.0
#   ## apply the transform to the inverse polynomials
#   rpcout.inverseColNum, rpcout.inverseColDen  =  poly_variable_change_out(rpcout.col_num, rpcout.col_den, a11, -b1)
#   rpcout.inverseLinNum, rpcout.inverseLinDen  =  poly_variable_change_out(rpcout.row_num, rpcout.row_den, a22, -b2)

   return rpcout


def test_me():
    import numpy as np
    import rpcm
    from rpcm.rpc_file_readers import read_rpc_file
    import rpc_crop
    r1 = rpcm.RPCModel(read_rpc_file('RPC_PHR1A_P_201309231105136_SEN_756965101-001.XML'))
    r2 = rpc_crop.rpc_apply_crop_to_rpc_model(r1, 10000,20000,2000,2000)
    
    #print "direct estimate error:"
    geo1 = np.array(r1.localization(11000,20100,10, return_normalized=False))
    geo2 = np.array(r2.localization(1000,100,10, return_normalized=False))
    print (geo1 - geo2)
    
    #print "inverse estimate error:"
    pix1 = np.array(r1.projection(geo1[0], geo1[1], geo1[2]))
    pix2 = np.array(r2.projection(geo1[0], geo1[1], geo1[2]))
    print (pix1 - pix2)
    
    r2.write_to_file('cropped.xml')


def crop_rpc_and_image(out_dir, img, rpc, rpc_ref, x, y, w, h):
    """
    Crops an image and its rpc. The ROI may be defined on another image.

    Args:
        out_dir: path to the output directory. The cropped image and rpc files
            will be written there.
        img: path to the input image
        rpc: path to the input rpc
        rpc_ref: path to the rpc file of the reference image
        x, y, w, h: 4 integers defining a rectangular ROI in the reference
            image
    """
    # get the rpc of the first image
    r = rpcm.RPCModel(read_rpc_file(rpc))

    # recompute the roi if the input image is not the reference image
    if rpc_ref is not rpc:
        r_ref = rpcm.RPCModel(read_rpc_file(rpc_ref))
        cfg = {'exogenous_dem':None , 'use_srtm':True, 'rpc_alt_range_scale_factor': 1.0, 'exogenous_dem_geoid_mode': False}
        x, y, w, h = rpc_utils.corresponding_roi(cfg, r_ref, r, x, y, w, h)

    # output filenames
    crop_rpc_and_image.counter += 1
    s = "_%02d" % crop_rpc_and_image.counter
    out_img_file = os.path.join(out_dir, "img%s.tif" % s)
    out_rpc_file = os.path.join(out_dir, "rpc%s.xml" % s)
    out_prv_file = os.path.join(out_dir, "prv%s.png" % s)

    # do the crop
    out_r = rpc_apply_crop_to_rpc_model(r, x, y, w, h)
    out_r.write_to_file(out_rpc_file)

    # generate a Geotiff
    out_tags = out_r.to_geotiff_dict()
    out_profile = {}
    with rasterio.open(img, "r") as f:
        out_img_data = f.read(window=((y, y + h), (x, x + w))).squeeze()
        out_profile = f.profile

    common.rasterio_write(out_img_file, out_img_data, profile=out_profile)

    # add RPC attributes
    with rasterio.open(out_img_file, "r+") as f:
        f.update_tags(ns="RPC", **out_tags)

    # save preview
    common.rasterio_write(out_prv_file, common.linear_stretching_and_quantization_8bit(out_img_data))


def main():

    # verify input
    if len(sys.argv) in [8, 10, 12]:
       out_dir = sys.argv[1]
       img1 = sys.argv[2]
       rpc1 = sys.argv[3]
       if len(sys.argv) == 8:
           x = float(sys.argv[4])
           y = float(sys.argv[5])
           w = float(sys.argv[6])
           h = float(sys.argv[7])
       elif len(sys.argv) == 10:
           img2 = sys.argv[4]
           rpc2 = sys.argv[5]
           x = float(sys.argv[6])
           y = float(sys.argv[7])
           w = float(sys.argv[8])
           h = float(sys.argv[9])
       else:
           img2 = sys.argv[4]
           rpc2 = sys.argv[5]
           img3 = sys.argv[6]
           rpc3 = sys.argv[7]
           x = float(sys.argv[8])
           y = float(sys.argv[9])
           w = float(sys.argv[10])
           h = float(sys.argv[11])
    else:
       print ("Tool to crop an image and its RPC.")
       print ("Incorrect syntax, use:")
       print ("  >  %s out_dir img1 rpc1 [img2 rpc2 [img3 rpc3]] x y w h" % sys.argv[0])
       print ("   generates files in the out_dir directory:  img_01.tif  rpc_01.tif [img_02.tif rpc_02.tif ... ] ")
       print ("   the image files are geotiffs and contain the rpc metadata so can be visualized in qgis")
       exit(1)

    try:
       os.stat(img1)
       os.stat(rpc1)

    except OSError:
       exit(1)

    # create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # do the crops
    crop_rpc_and_image.counter = 0 # used to name the output files
    crop_rpc_and_image(out_dir, img1, rpc1, rpc1, x, y, w, h)
    if 'img2' in locals() and 'rpc2' in locals():
        crop_rpc_and_image(out_dir, img2, rpc2, rpc1, x, y, w, h)
    if 'img3' in locals() and 'rpc3' in locals():
        crop_rpc_and_image(out_dir, img3, rpc3, rpc1, x, y, w, h)


if __name__ == '__main__': main()
