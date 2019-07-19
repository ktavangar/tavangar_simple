#!/usr/bin/env python
'''
This script will create a footprint file based on the fracdet and exptime maps 
Author: Nacho Sevilla (based on Eli Rykoff's code)
'''
import healpy as hp
import numpy as np
import sys
import os
from optparse import OptionParser

def make_footprint_file(syspath,minexp,minfrac,nside,toggle_Y, exptime):
    '''
    Create footprint file based on the exposure time and detection fraction
    '''
    footprint_accum = np.zeros(hp.nside2npix(nside),dtype=np.int32) 
    footprint = np.zeros(hp.nside2npix(nside),dtype=np.int32) 

    if toggle_Y == True:
        bands = ['g','r','i','z','Y']
        detfrac_filename = '%s/decals_grizY_o.%s_t.32768_coverfoot_EQU.fits.gz' % (syspath,nside)
    else:
        bands = ['g','r'] # could also have 'i' and 'z' here
        detfrac_filename = 'decals_dr7_pseudo_fracdet.fits' # from pseudo_fracdet code
        #detfrac_filename = '%s/decals_griz_o.%s_t.32768_coverfoot_EQU.fits.gz' % (syspath,nside)

    nbands = len(bands)
    try:
        detfrac = hp.read_map(detfrac_filename,nest=True)
    except IOError:
        print detfrac_filename,'not found'
        sys.exit()
    print 'Reading',detfrac_filename

    if exptime:
    #    for band in bands: # used if exposure time is important
    #        print 'Processing band',band
    #        exptime_filename = '%s/decals_%s_o.%s_t.32768_N_IMAGES_EQU.fits.gz' % (syspath,band,nside)
    #        try:
    #            exptime = hp.read_map(exptime_filename,nest=True)
    #        except IOError:
    #            print exptime_filename,'not found'        
    #        print 'Reading',exptime_filename
    #       ## for individual band coverage
    #        #detfrac_filename = '%s/y3a1_%s_o.%s_t.32768_frac_EQU.fits.gz' % (syspath,band,nside)
    #        #print 'Reading',detfrac_filename
    #        #detfrac = hp.read_map(detfrac_filename,nest=True)
    #        mask, = np.where((detfrac>minfrac) & (exptime>=minexp))
    #        footprint_accum[mask] = footprint_accum[mask] + 1
    #        mask, = np.where(footprint_accum==nbands) # must have passed the criteria for each filter
    else:
        mask, = np.where((detfrac>minfrac) # set to 0.5 by default
        footprint_accum[mask] = footprint_accum[mask] + 1

    footprint[mask] = 1
    write_footprint(footprint,syspath,minexp,toggle_Y)

def write_footprint(hpfoot,syspath,minexp,toggle_Y):
    '''
    Write footprint to file
    '''
    hpfoot_float = np.zeros(hpfoot.size,dtype=np.float32) + hp.UNSEEN
    ok,=np.where(hpfoot>0)
    hpfoot_float[ok] = hpfoot[ok].astype(np.float32)

    if toggle_Y == True:
        band_combination = 'grizY'
    else:
        band_combination = 'gr'        

    outfilename = '%s/decals_footprint_%s.fits' % (syspath,band_combination)
    #outfilename = '%s/decals_footprint_%s_%dexp.fits' % (syspath,band_combination,minexp)
    hp.write_map(outfilename,hpfoot_float,coord='C',nest=True)

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='/masking')
    parser.add_option("-Y","--grizY",action="store_true",dest="toggle_Y",help="Add Y band to requirement for footprint",default=False)
    parser.add_option("-t","--exp",action="store_true",dest="exptime",help="Use exposure time as a parameter",default=False)
    parser.add_option("-e","--minimum_exposures",type="int",dest="minexp",help="Minimum number of exposures",default=1)
    parser.add_option("-f","--minimum_detection_fraction",type="float",dest="minfrac",help="Minimum detection fraction",default=0.5)
    parser.add_option("-n","--nside",type="int",dest="nside",help="Healpix nside",default=1024)

    # Parse command line
    (options, args) = parser.parse_args()
    
    make_footprint_file(options.syspath,options.minexp,options.minfrac,options.nside,options.toggle_Y)

    print 'DONE'

if __name__ == '__main__':
    sys.exit(main())
    
