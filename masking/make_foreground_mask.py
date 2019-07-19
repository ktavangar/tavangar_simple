#!/usr/bin/env python
'''
This script will create a foreground objects map file using the Y3 Gold footprint 
Author: Nacho Sevilla (based Eli Rykoff's and Alex Drlica-Wagner's code)
'''
import healpy as hp
import numpy as np
#from astropy.io import fits
import fitsio
import sys
import os
from optparse import OptionParser
import copy
import re
import pandas
from mask_table import make_mask_table
import mask_bits

def make_foremask_file(workdir,footprint_filename,foremask_filename,nside,bits):
    
    ''' 
    Inputs:
    workdir - working directory
    footprint_filename - footprint file created with make_footprint.py
    foremask_filename - output file for the mask
    nside - 1024 for decals
    bits - from make_bits.py
    '''

    print 'Reading footprint in',footprint_filename
    
    #footprint = hp.read_map(footprint_filename,nest=True) ## this is too slow
    if footprint_filename.endswith('.fits') or footprint_filename.endswith('.fit') or footprint_filename.endswith('.fits.gz'):
        footprint = fitsio.read(footprint_filename,ext=1)['I'].ravel()
    else:
        print 'Filename extension not .fits, .fit or fits.gz Exiting...'
        sys.exit()
    foremask = np.zeros(hp.nside2npix(nside),dtype=np.int32)    
    pixarea = hp.nside2pixarea(nside,degrees=True)

    for b,d in bits.items():
        mask_name = d['name']
        if mask_name == 'Total area':
            continue # don't do the following lines and proceed to the next iteration
        print 'Analyzing',mask_name
        #if re.search('lmc',mask_name, re.IGNORECASE):
            #hpix,=np.where(footprint >= 0)
            #ra,dec=hp.pix2ang(nside,hpix,lonlat=True,nest=True)
            #bad,=np.where((ra > 60) & (ra < 100) & (dec > -70) & (dec < -58)) # this is the LMC region, N/A for decals
            #foremask[hpix[bad]] = foremask[hpix[bad]] | b
            #print mask_name,'masked',len(hpix[bad])*pixarea,'square degrees'
    else:
            badpix = compute_radmask_badpix(workdir,footprint,d)
            foremask[badpix] = foremask[badpix] | b
            print mask_name,'masked',len(badpix)*pixarea,'square degrees'

    print "Writing map to",foremask_filename
    hp.write_map(foremask_filename,foremask,nest=True,coord='C',dtype=np.int32)

def compute_radmask_badpix(workdir,footprint,mask):
    
    # this is for all foreground objects that are not in the LMC
    
    nside=hp.npix2nside(footprint.size)
    #nside = 1024

    #let's degrade the footprint mask to make radii on a coarser resolution, for speed
    #therefore the centers of some objects may lie outside of the hires footprint

    nside_ref=32
    footprint_low = hp.ud_grade(footprint,nside_ref,False,'NEST','NEST') #degrades the footprint

    # we now filter out the objects really outside the coarse footprint, and create radii for objects without them
    if mask['mag'] == None: 
        #hdu = fits.open(workdir+mask['filename'])
        #cat = hdu[1].data
        cat = fitsio.read(workdir+mask['filename'],ext=1) # need to have all these files
        pix_ext = hp.ang2pix(nside_ref,cat['ra'],cat['dec'],lonlat=True,nest=True)
        gd,=np.where(footprint_low[pix_ext] >= 0)
        cat=cat[gd]
    else:
        #hdu = fits.open(workdir+mask['filename'])
        #cat = hdu[1].data
        cat = fitsio.read(workdir+mask['filename'],ext=1)

        #if re.search('leda',mask['name'], re.IGNORECASE):
        #    print np.where(cat['RADIUS'] != 0)

        masknan = np.isnan(cat['ra'])
        cat = cat[~masknan]
        pix_ext = hp.ang2pix(nside_ref,cat['ra'],cat['dec'],lonlat=True,nest=True)
        gd,=np.where(footprint_low[pix_ext] >= 0)
        cat=cat[gd] # take away points that are negative

        m=(mask['r2']-mask['r1'])/(mask['m2']-mask['m1'])
        maskrad=m*cat[mask['mag']] - m*mask['m1'] + mask['r1']
        hi,=np.where(maskrad > mask['maxrad'])
        lo,=np.where(maskrad < mask['minrad'])
        maskrad[hi] = mask['maxrad']
        maskrad[lo] = mask['minrad']
        
        norad, = np.where(cat['rad'] == 0)
        cat['rad'][norad] = maskrad[norad]

    vec=hp.ang2vec(cat['ra'],cat['dec'],lonlat=True)
    map_temp = footprint.copy()

    print 'Defining foreground mask pixels...'
    for i in xrange(cat['rad'].size): #for each object in coarse footprint
        #check all pixels inside the avoidance radii + some cushion in high resolution
        pixint=hp.query_disc(nside,vec[i,:],(cat['rad'][i]+mask['cushion'])*np.pi/180.,inclusive=False,fact=8,nest=True)
        map_temp[pixint] = -1*(i+1)

    #badpix,=np.where((map_temp < 0) & (footprint >= 0))
    badpix,=np.where((map_temp < 0) & (footprint > 0)) #pixels marked bad should in the radius AND in the footprint

    return badpix

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='/masking')
    parser.add_option("--table", action="store_true", dest="toggle_table", default=False, help="Printout coverage table")
    parser.add_option("-n","--nside",type="int",dest="nside",help="Healpix nside",default=4096)
    parser.add_option("--footprint_filename",dest="footprint_filename",help="Footprint filename (input, fits format)",default='decals_footprint_gr.fits.gz')
    parser.add_option("--foremask_filename",dest="foremask_filename",help="Foreground mask filename (output, fits format)",default='decals_foreground_mask.fits.gz')

    # Parse command line
    (options, args) = parser.parse_args()

    workdir = options.syspath
    footprint_filename = workdir + options.footprint_filename
    foremask_filename = workdir + options.foremask_filename
    fracdet_filename = workdir + 'decals_dr7_pseudo_fracdet.fits'
    nside = options.nside
    mask_bits.init()
    bits = copy.deepcopy(mask_bits.FOREGROUND_BITS)

    if options.toggle_table:
        make_mask_table(workdir,footprint_filename,foremask_filename,fracdet_filename,nside,'foreground')
    else:
        make_foremask_file(workdir,footprint_filename,foremask_filename,nside,bits)

if __name__ == '__main__':
    sys.exit(main())
