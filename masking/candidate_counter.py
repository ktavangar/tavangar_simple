'''
Author: Kiyan Tavangar
Purpose: Using a given candidate list, create a graph of the candidates that would be excluded under certain types of masking
'''


import healpy as hp
import numpy as np
from astropy.io import fits
import fitsio
import sys
import os
from optparse import OptionParser
import copy
import matplotlib.pyplot as plt
import re
import pandas
import mask_bits

def count_candidates(workdir,candidate_list,foreground_mask,nside):
    if candidate_list.endswith('.fits') or candidate_list.endswith('.fit') or candidate_list.endswith('.fits.gz'):
        candidates = fitsio.read(candidate_list)
    else:
        print 'Filename extension not .fits, .fit or fits.gz Exiting...'
        sys.exit()

    foreground_mask = fitsio.read(foreground_mask)
    foreground_mask = foreground_mask['T'] #This creates one long array of the values associated with each pixel
    foreground_mask = np.concatenate(foreground_mask)

    sig = candidates['SIG']
    ra = candidates['RA']
    dec = candidates['DEC']
    mod = candidates['MODULUS']
    r = candidates['r']
    N_OBS = candidates['N_OBS']
    N_OBS_HALF = candidates['N_OBS_HALF']
    N_MODEL = candidates['N_MODEL']
    MC_SOURCE_ID = candidates['MC_SOURCE_ID']
    pix = hp.ang2pix(nside,ra,dec,lonlat=True, nest=True)
   
    cut_milky_way = np.where(foreground_mask[pix] & 0b00010000, False, True)
    cut_lmc = np.where(foreground_mask[pix] & 0b00001000, False, True)
    cut_leda = np.where(foreground_mask[pix] & 0b00000010, False, True)
    cut_faint_stars = np.where(foreground_mask[pix] & 0b00000001, False, True)
    cut_bright_stars = np.where(foreground_mask[pix] & 0b00000100, False, True)
    cut_yale = np.where(foreground_mask[pix] & 0b00100000, False, True)

    bins = np.arange(0., 40., 0.5) # step=0.5, sigs max out at 37.5
    axis_label = 'SIG'
    
    plt.ylim(0.5, 40000)
    plt.yscale('log')
    plt.hist(sig, bins=bins, color='red', histtype='step', cumulative=-1, linewidth=1.2, label='All')
    plt.hist(sig[cut_milky_way], bins=bins, color='blue', histtype='step', cumulative=-1, linewidth=1.2, label= r'Cut Milky Way')
    plt.hist(sig[cut_milky_way & cut_lmc], bins=bins, color='green', histtype='step', cumulative=-1, linewidth=1.2, label='cut lmc')
    plt.hist(sig[cut_milky_way & cut_lmc & cut_leda], bins=bins, color='darkorange', histtype='step', cumulative=-1, linewidth=1.2, label='cut leda galaxies')
    plt.hist(sig[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale], bins=bins, color='purple', histtype='step', cumulative=-1, linewidth=1.2, label='cut stars')
    plt.ion()
    plt.show()

    #starting process of writing the reduced file
    sig = sig[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    ra = ra[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    dec = dec[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    mod = mod[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    r = r[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    N_OBS = N_OBS[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    N_OBS_HALF = N_OBS_HALF[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    N_MODEL = N_MODEL[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]
    MC_SOURCE_ID = MC_SOURCE_ID[cut_milky_way & cut_lmc & cut_leda & cut_faint_stars & cut_bright_stars & cut_yale]

    c0 = fits.Column(name='SIG',          format='E', array=sig)
    c1 = fits.Column(name='RA',        format='E', array=ra)
    c2 = fits.Column(name='DEC',        format='E', array=dec)
    c3 = fits.Column(name='MODULUS',      format='E', array=mod)
    c4 = fits.Column(name='r',            format='E', array=r)
    c5 = fits.Column(name='N_OBS',        format='E', array=N_OBS)
    c6 = fits.Column(name='N_OBS_HALF',   format='E', array=N_OBS_HALF)
    c7 = fits.Column(name='N_MODEL',      format='E', array=N_MODEL)
    c8 = fits.Column(name='MC_SOURCE_ID', format='E', array=MC_SOURCE_ID)
    
    t = fits.BinTableHDU.from_columns([c0, c1, c2, c3, c4, c5, c6, c7, c8])
    t.writeto('reduced_candidate_list0.fits', overwrite=True)

    # Diagnostic output
    data = fitsio.read('reduced_candidate_list0.fits')
    print("{} hotspots found.").format(len(data))
    cut_0 = (data['SIG'] > 5.5)
    print("{} hotspots found with SIG > 5.5.").format(len(data[cut_0]))
    cut_1 = (data['SIG'] > 10)
    print("{} hotspots found with SIG > 10.").format(len(data[cut_1]))
    cut_2 = (data['SIG'] > 15)
    print("{} hotspots found with SIG > 15.").format(len(data[cut_2]))
    cut_3 = (data['SIG'] > 20)
    print("{} hotspots found with SIG > 20.").format(len(data[cut_3]))
    cut_4 = (data['SIG'] > 25)
    print("{} hotspots found with SIG > 25.").format(len(data[cut_4]))
    cut_5 = (data['SIG'] >= 37.5)
    print("{} hotspots found with SIG >= 37.55").format(len(data[cut_5]))


def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='/home/s1/tavangar/tav_simple/')
    parser.add_option("--candidate_list",dest="candidate_list",help="Candidate List (input, fits format)",default='candidate_list.fits')
    parser.add_option("--foreground_mask",dest="foreground_mask",help="Foreground Mask (input, fits format)",default='all_sky_foreground_mask0.fits')
    parser.add_option("-n","--nside",type="int",dest="nside",help="Healpix nside",default=4096)

    # Parse command line	
    (options, args) = parser.parse_args()	
    workdir = options.syspath
    candidate_list = workdir + options.candidate_list
    foreground_mask = '/home/s1/tavangar/data/all_sky/' + options.foreground_mask
    nside = options.nside

    count_candidates(workdir,candidate_list,foreground_mask,nside)

if __name__ == '__main__':
    sys.exit(main())
