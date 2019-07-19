#!/usr/bin/env python
"""
Skim data from a file.
"""
import os,sys
import logging
import fitsio
from collections import OrderedDict as odict
import numpy as np

from matplotlib import mlab

import utils

COLUMNS = ['COADD_OBJECT_ID', 'RA', 'DEC',
           'ALPHAWIN_J2000','DELTAWIN_J2000',
           'SOF_PSF_MAG_G','SOF_PSF_MAG_R','SOF_PSF_MAG_I',
           'SOF_PSF_MAG_ERR_G','SOF_PSF_MAG_ERR_R','SOF_PSF_MAG_ERR_I',
           'DELTA_MAG_V3_G','DELTA_MAG_V3_R','DELTA_MAG_V3_I',
           'DELTA_MAG_CHROM_G','DELTA_MAG_CHROM_R','DELTA_MAG_CHROM_I',
           'A_SED_SFD98_G','A_SED_SFD98_R','A_SED_SFD98_I',
           'SOF_CM_MAG_CORRECTED_G','SOF_CM_MAG_CORRECTED_R',
           'SOF_CM_MAG_ERR_G','SOF_CM_MAG_ERR_R',
           'WAVG_MAG_PSF_G','WAVG_MAG_PSF_R','WAVG_MAG_PSF_I',
           'WAVG_MAGERR_PSF_G','WAVG_MAGERR_PSF_R','WAVG_MAGERR_PSF_I',
           'WAVG_SPREAD_MODEL_I','WAVG_SPREADERR_MODEL_I',
           'SOF_CM_T','SOF_CM_T_ERR',
           'EXTENDED_CLASS_MASH_SOF',
           'FLAGS_GOLD','SOF_FLAGS','FLAGS_BADREGIONS',
           'FLAGS_FOOTPRINT','FLAGS_FOREGROUND',
           ]

# full output
#OUTPUT = COLUMNS + ['SOF_PSF_MAG_CORRECTED_G','SOF_PSF_MAG_CORRECTED_R',
#                    'WAVG_MAG_PSF_CORRECTED_G','WAVG_MAG_PSF_CORRECTED_R']

# small output 
OUTPUT = ['COADD_OBJECT_ID', 'RA', 'DEC',
          'SOF_PSF_MAG_CORRECTED_G','SOF_PSF_MAG_CORRECTED_R','SOF_PSF_MAG_CORRECTED_I',
          'SOF_PSF_MAG_ERR_G','SOF_PSF_MAG_ERR_R','SOF_PSF_MAG_ERR_I',
          'A_SED_SFD98_G','A_SED_SFD98_R','A_SED_SFD98_I',
          'WAVG_MAG_PSF_G','WAVG_MAG_PSF_R','WAVG_MAG_PSF_I',
          #'WAVG_MAG_PSF_CORRECTED_G','WAVG_MAG_PSF_CORRECTED_R','WAVG_MAG_PSF_CORRECTED_I',
          'WAVG_MAGERR_PSF_G','WAVG_MAGERR_PSF_R','WAVG_MAGERR_PSF_I',
          'WAVG_SPREAD_MODEL_I','WAVG_SPREADERR_MODEL_I',
          'SOF_CM_T','SOF_CM_T_ERR',
          'FLAGS_GOLD','EXTENDED_CLASS_MASH_SOF',
          ]

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',
                        help="input FITS file")
    parser.add_argument('outfile',
                        help="output FITS file")
    parser.add_argument('--class',choices=['star','galaxy','all'],dest='cls',
                        help='object class to skim')
    parser.add_argument('--ext',nargs=2,default=[0,4],type=int,
                        help='object class to skim')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output files')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='verbosity')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)

    if os.path.exists(args.outfile) and not args.force:
        logging.info("Found %s; skipping..."%args.outfile)
        sys.exit(0)
    else:
        outdir = os.path.dirname(args.outfile)
        if outdir: utils.mkdir(outdir)

    logging.info("Reading %s..."%args.infile)
    data,header=fitsio.read(args.infile,header=True,columns=COLUMNS)

    extclass = data['EXTENDED_CLASS_MASH_SOF']

    # Do the selection
    logging.info("Skimming objects with: %i <= ext <= %i"%(args.ext[0],args.ext[1]))
    sel = ((extclass >= args.ext[0]) & (extclass <= args.ext[1]))

    # Only objects with measurements (SOF sentinel = 37.5)
    sel &= ( (data['SOF_PSF_MAG_G'] < 35) | (data['WAVG_MAG_PSF_G'] < 35) )
    sel &= ( (data['SOF_PSF_MAG_R'] < 35) | (data['WAVG_MAG_PSF_R'] < 35) )
    #sel &= ( (data['FLAG_FOREGROUND'] & 16) == 0 )
    #sel &= (data['FLAG_FOOTPRINT'] == 1)
    sel &= ( (data['FLAGS_GOLD'] & 0b111100) == 0 )
    data = data[sel]

        
    psf_sfd_g = data['SOF_PSF_MAG_G'] + data['DELTA_MAG_V3_G'] \
        + data['DELTA_MAG_CHROM_G'] - data['A_SED_SFD98_G']
    psf_sfd_r = data['SOF_PSF_MAG_R'] + data['DELTA_MAG_V3_R'] \
        + data['DELTA_MAG_CHROM_R'] - data['A_SED_SFD98_R']
    psf_sfd_i = data['SOF_PSF_MAG_I'] + data['DELTA_MAG_V3_I'] \
        + data['DELTA_MAG_CHROM_I'] - data['A_SED_SFD98_I']

    wavg_sfd_g = data['WAVG_MAG_PSF_G'] + data['DELTA_MAG_V3_G'] \
        + data['DELTA_MAG_CHROM_G'] - data['A_SED_SFD98_G']
    wavg_sfd_r = data['WAVG_MAG_PSF_R'] + data['DELTA_MAG_V3_R'] \
        + data['DELTA_MAG_CHROM_R'] - data['A_SED_SFD98_R']
    wavg_sfd_i = data['WAVG_MAG_PSF_I'] + data['DELTA_MAG_V3_I'] \
        + data['DELTA_MAG_CHROM_I'] - data['A_SED_SFD98_I']

    new = odict([
            ('SOF_PSF_MAG_CORRECTED_G',psf_sfd_g),
            ('SOF_PSF_MAG_CORRECTED_R',psf_sfd_r),
            ('SOF_PSF_MAG_CORRECTED_I',psf_sfd_i),
            #('WAVG_MAG_PSF_CORRECTED_G',wavg_sfd_g),
            #('WAVG_MAG_PSF_CORRECTED_R',wavg_sfd_r),
            #('WAVG_MAG_PSF_CORRECTED_I',wavg_sfd_i),
            ])
    data = mlab.rec_append_fields(data,new.keys(),new.values())
    
    #drop = []
    #if drop:
    #    data = mlab.rec_drop_fields(data,drop)

    if len(data):
        logging.info("Writing %s..."%args.outfile)
        fitsio.write(args.outfile,data[OUTPUT],clobber=True)
    else:
        logging.warn("No data passing cuts; skipping...")
