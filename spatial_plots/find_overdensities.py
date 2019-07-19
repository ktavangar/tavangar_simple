import numpy as np
import healpy as hp
import sys
import os
import query_image

mag_cut = float(sys.argv[1])
print(mag_cut)
str_mag_cut = int(10 * mag_cut)

print('Loading data')
data0 = np.load('/home/s1/tavangar/diff_cutoffs_data/hpxstar_data_{}.npy'.format(str_mag_cut))


# determine which pixels will be flagged as over dense
# Find the mean of the pixels that are not zero
nside = 1024
data = np.where(data0 != 0)
ave = np.mean(data0[data])
std = np.std(data0[data])
print(ave)
print(std)

# different cuts need different thresholds
if mag_cut == 22.5:
    threshold = round(10**2.1, 0)
elif mag_cut == 22.0:
    threshold = round(10**2.05, 0)
elif mag_cut == 21.5:
    threshold = round(10**2.0, 0)
elif mag_cut == 21.0:
    threshold = round(10**1.95, 0)
elif mag_cut ==20.5:
    threshold = round(10**1.9, 0)

mask = np.where(data0 >= threshold)

theta, phi = hp.pix2ang(nside, mask)
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)

count = len(ra[0])
print(count)

print('Creating coordinates')
coord = np.zeros((count, 1, 3))
for i in range(count):
    pix = hp.ang2pix(nside, theta[0,i], phi[0,i])
    value = data0[pix]
    coord[i] = [ra[0,i], dec[0,i], value]
print(coord)
print(len(coord))
print('Saving Coordinates')
np.save('coord_data_thrsh={}_cut={}.npy'.format(threshold, str_mag_cut), coord)



save_dir = 'save_dirs/save_dir_{}'.format(str_mag_cut)
plots = os.listdir(save_dir)

for i in range(count):
    image_file = 'image_{:0.2f}_{:0.2f}.png'.format(ra[0,i], dec[0,i])
    if image_file not in plots:
        print('Image not found; retrieving image for ({}, {})'.format(ra[0,i], dec[0,i]))
image_file = 'image_{:0.2f}_{:0.2f}.png'.format(ra[0,i], dec[0,i])
query_image.retrieve_image(image_file, ra[0,i], dec[0,i], save_dir)
image = '{}/{}'.format(save_dir, image_file)
iname = image.strip('.png')
