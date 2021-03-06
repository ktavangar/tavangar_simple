# Set the backend first!
import matplotlib
matplotlib.use('Agg')

import fitsio as fits
import numpy as np
import pylab as plt
import glob
import yaml
import ugali
import ugali.candidate.associate

import healpy as hp

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    candidate_list = cfg[survey]['candidate_list']
    fracdet = cfg[survey]['fracdet']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

candidate_list = cfg[survey]['candidate_list']

data = fits.read(candidate_list)
data = data[data['SIG'] > 15]

fracdet = fits.read('{}_pseudo_fracdet.fits'.format(survey))
#fracdet = fits.read(fracdet)
fracdetrav = fracdet['I'].ravel()
#fracdetrav = fracdet['SIGNAL'].ravel()

plt.figure()

lonmin = -180
lonmax = 180
#lonmin = 200-180
#lonmax = 260-180
latmin = min(data[basis_2])
#latmax = latmin+20
latmax = max(data[basis_2])

hp.cartview(fracdetrav, lonra=[lonmin, lonmax], latra=[latmin, latmax], return_projected_map=True, cmap='binary')
hp.projscatter(data[basis_1], data[basis_2], lonlat=True, edgecolor='none', s=0.5, c='red')

plt.xlabel(basis_1)
plt.ylabel(basis_2)
plt.savefig("{}_fracdet_check.png".format(survey), bbox_inches='tight', dpi=1000)
plt.close()
