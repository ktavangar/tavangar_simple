import numpy as np
import fitsio
import glob
import healpy as hp

import os
from sklearn.mixture import GMM
#import seaborn as sns
import matplotlib.pyplot as plt
import pylab
import sys

import utils

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


#fnames_healpix = glob.glob('/data/des40.b/data/decals/dr7/healpix/*.fits')
fnames_skim = glob.glob('/data/des40.b/data/decals/dr7/skim/*.fits')
nside=1024
#hpxmap_healpix = np.zeros(hp.nside2npix(nside))
#hpxmap_skim = np.zeros(hp.nside2npix(nside))
hpxmap_stars = np.zeros(hp.nside2npix(nside))

mag_cut = float(sys.argv[1])
str_mag_cut = int(10 * mag_cut)

if os.path.exists('/home/s1/tavangar/data/spatial_plots/decals_dr7/diff_cutoffs_data/hpxstar_data_{}.npy'.format(str_mag_cut)):
    print('Loading data')
    data = np.load('/home/s1/tavangar/data/spatial_plots/decals_dr7/diff_cutoffs_data/hpxstar_data_{}.npy'.format(str_mag_cut))
    hp.mollzoom(np.log10(data), title = 'STARS SPATIAL MAP WITH {} MAG_G CUT'.format(str_mag_cut))
    plt.ion()
    plt.show()
else:
    print('saving data')
    for f in fnames_skim:
        data = fitsio.read(f, columns=['RA','DEC', 'EXTENDED_CLASS', 'MAG_G'])
        to_del_class = np.where(data['EXTENDED_CLASS'] == 1) # filter out all objects taht are not stars
        data = np.delete(data, to_del_class)
    
        to_del_mag = np.where(data['MAG_G'] > mag_cut) # filter out all objects above a certain magnitude
        d = np.delete(data, to_del_mag)
                  
        pixels_stars = hp.ang2pix(nside,np.radians(-d['DEC']+90.),np.radians(d['RA']))
        pix_stars,cts_stars = np.unique(pixels_stars,return_counts=True)
        hpxmap_stars[pix_stars] = cts_stars
    np.save('/home/s1/tavangar/data/spatial_plots/decals_dr7/diff_cutoffs_data/hpxstar_data_{}.npy'.format(str_mag_cut), hpxmap_stars)
