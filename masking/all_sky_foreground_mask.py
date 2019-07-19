#!/usr/bin/env python
'''
This script will create a foreground objects map file for the entire sky
Author: Kiyan Tavangar (based Eli Rykoff's and Alex Drlica-Wagner's code)
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
#from mask_table import make_mask_table
import mask_bits

nside = 4096
footprint_all_sky = np.zeros(hp.nside2npix(nside),dtype=np.int32)
footprint_all_sky = footprint_all_sky + 1 # footprint is the entire sky

def make_foremask_file(workdir,footprint,milky_way_file,foremask_filename,nside,bits):
    
    ''' 
    Inputs:
    workdir - working directory
    foremask_filename - output file for the mask
    nside - 4096
    bits - from make_bits.py
    '''

    foremask = np.zeros(hp.nside2npix(nside),dtype=np.int32)    
    hpix, = np.where(footprint_all_sky >= 0)
    print(hpix)
    print(len(hpix))
    pixarea = hp.nside2pixarea(nside,degrees=True)

    for b,d in bits.items():
        mask_name = d['name']
        if mask_name == 'Total area':
            continue # don't do the following lines and proceed to the next iteration
        print 'Analyzing',mask_name
        if re.search('milky',mask_name, re.IGNORECASE):
            milky_way = hp.read_map(milky_way_file, nest=True)
            bad,=np.where(milky_way > 0.2)
            foremask[hpix[bad]] = foremask[hpix[bad]] | b
            print mask_name,' masked ',len(hpix[bad])*pixarea,'square degrees'
        elif re.search('the',mask_name, re.IGNORECASE):
            ra = [80.89, 13.19] # ra for center of LMC and SMC respectively
            dec = [-69.76, -72.83] # dec for center of LMC and SMC respectively
            vec=hp.ang2vec(ra,dec,lonlat=True)
            print(vec)
            #bad,=np.where((ra > 60) & (ra < 100) & (dec > -70) & (dec < -58)) # this is the LMC regio
            #foremask[hpix[bad]] = foremask[hpix[bad]] | b
            pixint_lmc = hp.query_disc(nside,vec[0,:],15*np.pi/180.,inclusive=False,fact=8,nest=True) # create a radius of 15 degrees
            pixint_smc = hp.query_disc(nside,vec[1,:],10*np.pi/180.,inclusive=False,fact=8,nest=True) # create a radius of 10 degrees
            map_temp = footprint.copy()
            map_temp[pixint_lmc] = -1 # make all these pixels negative
            map_temp[pixint_smc] = -1
            badpix,=np.where(map_temp < 0)

            foremask[badpix] = foremask[badpix] | b

            print mask_name,' masked ',len(badpix)*pixarea,'square degrees'
        else:
            badpix = compute_radmask_badpix(workdir,footprint,d)
            foremask[badpix] = foremask[badpix] | b
            print mask_name,' masked ',len(badpix)*pixarea,'square degrees'

    print "Writing map to ",foremask_filename
    hp.write_map(foremask_filename,foremask,nest=True,coord='C',dtype=np.int32)
    
def compute_radmask_badpix(workdir,footprint,mask):
    
    # this is for all foreground objects that are not in the LMC
    
    #nside=hp.npix2nside(footprint.size)
    nside = 4096
    filename = workdir + mask['filename']

    #let's degrade the footprint mask to make radii on a coarser resolution, for speed
    #therefore the centers of some objects may lie outside of the hires footprint

    #nside_ref=32
    #footprint_low = hp.ud_grade(footprint,nside_ref,False,'NEST','NEST') #degrades the footprint

    # we now filter out the objects really outside the coarse footprint, and create radii for objects without them
    if mask['mag'] == None:
        cat = fitsio.read(workdir+mask['filename'],ext=1) # need to have all these files
    else:
        #hdu = fits.open(workdir+mask['filename'])
        #cat = hdu[1].data
        cat = fitsio.read(workdir+mask['filename'],ext=1)
        if re.search('leda',mask['name'], re.IGNORECASE):
            print np.where(cat['radius'] != 0)

        masknan = np.isnan(cat['ra'])
        cat = cat[~masknan]
        #pix_ext = hp.ang2pix(nside_ref,cat['ra'],cat['dec'],lonlat=True,nest=True)
        #gd,=np.where(footprint_low[pix_ext] >= 0)
        #cat=cat[gd] # take away points that are negative

        m=(mask['r2']-mask['r1'])/(mask['m2']-mask['m1'])
        maskrad=m*cat[mask['mag']] - m*mask['m1'] + mask['r1']
        hi,=np.where(maskrad > mask['maxrad'])
        lo,=np.where(maskrad < mask['minrad'])
        maskrad[hi] = mask['maxrad']
        maskrad[lo] = mask['minrad']
        
        norad, = np.where(cat['radius'] == 0)
        cat['radius'][norad] = maskrad[norad]

    vec=hp.ang2vec(cat['ra'],cat['dec'],lonlat=True)
    map_temp = footprint.copy()

    print 'Defining foreground mask pixels...'
    for i in xrange(cat['radius'].size): #for each object in coarse footprint
        #check all pixels inside the avoidance radii + some cushion in high resolution
        pixint=hp.query_disc(nside,vec[i,:],(cat['radius'][i]+mask['cushion'])*np.pi/180.,inclusive=False,fact=8,nest=True)
        map_temp[pixint] = -1*(i+1) # make all these pixels negative

    badpix,=np.where(map_temp < 0)

    return badpix

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='/home/s1/tavangar/data/all_sky/')
    parser.add_option("--table", action="store_true", dest="toggle_table", default=False, help="Printout coverage table")
    parser.add_option("-n","--nside",type="int",dest="nside",help="Healpix nside",default=4096)
    parser.add_option("--milky_way_file",dest="milky_way_file",help="Milky Way Filename (input, fits format)",default='ebv_sfd98_fullres_nside_4096_ring_equatorial.fits')
    #parser.add_option("--footprint_filename",dest="footprint_filename",help="Footprint filename (input, fits format)",default='decals_footprint_gr.fits.gz')
    parser.add_option("--foremask_filename",dest="foremask_filename",help="Foreground mask filename (output, fits format)",default='all_sky_foreground_mask.fits')

    # Parse command line
    (options, args) = parser.parse_args()

    workdir = options.syspath
    #footprint_filename = workdir + options.footprint_filename
    milky_way_file = workdir + options.milky_way_file
    foremask_filename = workdir + options.foremask_filename
    nside = options.nside
    mask_bits.init()
    bits = copy.deepcopy(mask_bits.FOREGROUND_BITS)

    if options.toggle_table:
        make_mask_table(workdir,footprint,foremask_filename,fracdet_filename,nside,'foreground')
    else:
        make_foremask_file(workdir,footprint_all_sky,milky_way_file,foremask_filename,nside,bits)

if __name__ == '__main__':
    sys.exit(main())
