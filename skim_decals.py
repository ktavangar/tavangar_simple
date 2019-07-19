#!/usr/bin/env python
"""
Skim data from a file.
"""
import os,sys
import logging
import fitsio
from collections import OrderedDict as odict
import numpy as np
import numpy.lib.recfunctions as recfn

from matplotlib import mlab

import utils

COLUMNS = ['BRICKID', 'OBJID', 'RA', 'DEC',
           'FLUX_G','FLUX_R','FLUX_Z',
           'FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z',
           'SHAPEEXP_R','DCHISQ','TYPE',
           'ANYMASK_G','ANYMASK_R','ANYMASK_Z',
           'FRACFLUX_G','FRACFLUX_R','FRACFLUX_Z',
           'NOBS_G','NOBS_R','NOBS_Z','EBV',
           'PARALLAX','PMRA','PMDEC','PMRA_IVAR','PMDEC_IVAR',
           'MW_TRANSMISSION_G','MW_TRANSMISSION_R','MW_TRANSMISSION_Z',
           ]

OUTPUT = ['UID'] + COLUMNS + \
         ['EXTENDED_CLASS',
          'MAG_G','MAG_R','MAG_Z',
          'MAG_ERR_G','MAG_ERR_R','MAG_ERR_Z',
          'MAG_SFD_G','MAG_SFD_R','MAG_SFD_Z',
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
    parser.add_argument('--ext',nargs=2,default=[0,1],type=int,
                        help='object class to skim (0=stars; 1=galaxies)')
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

    # Do the stellar selection
    # Switch so stars = 0; galaxies = 1
    extclass = (data['TYPE'] != 'PSF ')
    # Could use SHAPEEXP_R and/or DCHISQ...

    logging.info("Skimming objects with: %i <= ext <= %i"%(args.ext[0],args.ext[1]))
    sel = ((extclass >= args.ext[0]) & (extclass <= args.ext[1]))

    # Only objects with measurements

    # Selection objects passing quality cuts
    sel &= ((data['FLUX_G'] > 0) & (data['FLUX_R'] > 0) )
    sel &= (data['ANYMASK_G']==0) & (data['ANYMASK_R']==0) 
    sel &= (data['FRACFLUX_G'] < 0.05) & (data['FRACFLUX_R'] < 0.05) 

    # Select objects with >10sigma detection in r-band
    #sel &= (data['RFPSFMAGERR'] < 0.1)

    extclass = extclass[sel]
    data = data[sel]

    # Calculate magnitudes
    mag_g = -2.5*np.log10(data['FLUX_G']) + 22.5
    magerr_g = 2.5/np.log(10) * data['FLUX_IVAR_G']**-0.5 / data['FLUX_G']

    mag_r = -2.5*np.log10(data['FLUX_R']) + 22.5
    magerr_r = 2.5/np.log(10) * data['FLUX_IVAR_R']**-0.5 / data['FLUX_R']

    mag_z = -2.5*np.log10(data['FLUX_Z']) + 22.5
    magerr_z = 2.5/np.log(10) * data['FLUX_IVAR_Z']**-0.5 / data['FLUX_Z']

    # Extinction
    mag_sfd_g = -2.5*np.log10(data['FLUX_G']*data['MW_TRANSMISSION_G']) + 22.5
    mag_sfd_r = -2.5*np.log10(data['FLUX_R']*data['MW_TRANSMISSION_R']) + 22.5
    mag_sfd_z = -2.5*np.log10(data['FLUX_Z']*data['MW_TRANSMISSION_Z']) + 22.5
   
    # Unique ID
    ZEROSTR = '%06d%05d'
    uid = np.array([ZEROSTR%(data['BRICKID'][i],data['OBJID'][i]) for i in xrange(len(data))]
    new = odict([
            ('UID', uid),
            ('EXTENDED_CLASS', extclass),
            ('MAG_G', mag_g),
            ('MAG_R', mag_r),
            ('MAG_Z', mag_z),

            ('MAG_ERR_G', magerr_g),
            ('MAG_ERR_R', magerr_r),
            ('MAG_ERR_Z', magerr_z),

            ('MAG_SFD_G', mag_sfd_g),
            ('MAG_SFD_R', mag_sfd_r),
            ('MAG_SFD_Z', mag_sfd_z),
            ])
    data = recfn.rec_append_fields(data,new.keys(),new.values())
    
    if len(data):
        logging.info("Writing %s..."%args.outfile)
        fitsio.write(args.outfile,data[OUTPUT],clobber=True)
    else:
        logging.warn("No data passing cuts; skipping...")
