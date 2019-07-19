import healpy as hp
import numpy as np
import fitsio
import sys
import os
import mask_bits

def make_mask_table(workdir,footprint_filename,mask_filename,fracdet_filename,nside,identifier):

    print 'Making mask table'
    print 'Reading footprint'
    #footprint = hp.read_map(footprint_filename, nest=True) 
    footprint = fitsio.read(footprint_filename,ext=1)['I'].ravel()
    print 'Reading mask'
    #foremask = hp.read_map(foremask_filename, nest=True) 
    mask = fitsio.read(mask_filename,ext=1)['I'].ravel()
    print 'Reading detection fraction'
#    fracdet = hp.read_map(fracdet_filename, nest=True)
    fracdet = fitsio.read(fracdet_filename)
    fracdet_exp = np.zeros(hp.nside2npix(4096))
    for c,i in enumerate(fracdet):
        if c % 1000000 == 0:
            print c
        fracdet_exp[fracdet[c][0]] = fracdet[c][1]
##    if os.path.isfile(fracdet_full_filename):
##        fracdet = fitsio.read(fracdet_full_filename,ext=1)['I'].ravel()
##    else:
#    fracdet = np.full(hp.nside2npix(4096),hp.UNSEEN)
#    print 'Filling fracdet full-sky map'
#    for p in range(len(fracdet_partial)):
#        fracdet[fracdet_partial[p]['PIXEL']] = fracdet_partial[p]['SIGNAL']
###    fracdet_full_filename = workdir+'y3a2_griz_fracdet_full.fits'
###    hp.write_map(fracdet_full_filename,fracdet,nest=True,coord='C',dtype=np.int32)

    pixarea = hp.nside2pixarea(nside,degrees=True)
    griz = (footprint >= 1)
    if identifier == 'bad':
        bt = mask_bits.BAD_BITS
    elif identifier == 'foreground':
        bt = mask_bits.FOREGROUND_BITS
    else:
        print 'Identifier',identifier,'not found when making map'
        sys.exit()

    for b in bt.keys():
        grizarea = (~(mask >= b) * griz).sum() * pixarea
        fracarea = (~(mask >= b) * griz * fracdet_exp).sum() * pixarea
        #maskarea = (((mask & b) == b) * griz * fracdet_exp).sum() * pixarea
        maskarea = 0.0
        print '| %4s | %8.2f | %8.2f | %8.2f |'%(b,maskarea,fracarea,grizarea)      
    
