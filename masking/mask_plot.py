import numpy as np
import healpy as hp
import fitsio as fits
import matplotlib.pyplot as plt

foreground_mask = fits.read('/home/s1/tavangar/data/all_sky/all_sky_foreground_mask.fits')
foreground_mask = np.concatenate(foreground_mask['T'])


hp.mollzoom(np.log10(foreground_mask+1), nest=True) #log scale (approximate), had to add 1 so there wasn't an error
#hp.mollzoom(foreground_mask, nest=True) # linear scale
plt.ion()
plt.show()
