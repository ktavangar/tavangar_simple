Author__ = "Kiyan Tavnagar"
# Adapted from Sidney Mau

import glob
import numpy as np
import pandas as pd
import yaml
import fitsio as fits
import os
import query_image
import sys

mag_cut = float(sys.argv[1])
str_mag_cut = int(10*mag_cut)

save_dir = 'save_dirs/save_dir_{}'.format(str_mag_cut)

outfile = 'overdense_regions_web_{}.html'.format(str_mag_cut)
plots = os.listdir(save_dir)


INDEX = """
<html> 
  <head>
    <title>Candidates</title>
</head>
<body> 
%(table)s
</body>
</html>
"""

TABLE = """
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: center;">
      <th>Candidate</th>
      <th>DECALS-DR7 Image</th>
<!--      <th>Diagnostic Plots</th> -->
    </tr>
  </thead>
  <tbody>
%(rows)s
  </tbody>
"""

ROW = """
    <tr>
      <th>%(name)s</th>
      <!--<td> <a id="%(fname)s"></a><a href="%(fname)s.png"><img src="%(fname)s.png" width="400"></a></td>-->
      <td> <a id="%(fname)s"></a><a href="%(fname)s.png"><img src="%(fname)s.png"></a></td>
    </tr>  
"""

# different cuts need different thresholds
if mag_cut == 22.5:
        threshold = round(10**2.1, 0)
elif mag_cut == 21.5:
        threshold = round(10**2.0, 0)
elif mag_cut == 21.0:
        threshold = round(10**1.95, 0)
elif mag_cut ==20.5:
        threshold = round(10**1.9, 0)
 
print("Loading coordinates")
coord = np.load('coord_data_thrsh={}_cut={}.npy'.format(threshold, str_mag_cut))
print(coord)
ra = coord[:, :,0]
dec = coord[:,:,1]
value = coord[:,:,2]
count = len(ra)
def create_entry(ra, dec, value):
    """Create a diagnostic row for candidate"""
    ra = round(ra,2)
    dec = round(dec,2)
    image_file = 'image_{}_{}.png'.format(str(ra), str(dec))
    image = '{}/{}'.format(save_dir, image_file)
    iname = image.strip('.png')
    image_url = "http://legacysurvey.org/viewer?ra={0}&dec={1}&zoom={2}&layer=decals-dr7".format(str(ra), str(dec), 11)

    #tablerow = ROW%dict(name=plot.strip('.png'), fname=save_dir+'/'+plot.strip('.png'))
    tablerow = """
        <tr>
          <th><br>(RA, Dec) = ({}, {})<br>Value = {}<br></th>
          <td> <a id={}></a><a href={}><img src={} width="300"></a></td>
        </tr>  
    """.format(str(ra), str(dec), value, iname, image_url, image)

    return(tablerow)

def create_index_html(filename):
    """Create the index.html"""
    #tablerows = [ROW%dict(name=plot.strip('.png'), fname=save_dir+'/'+plot.strip('.png')) for plot in plots if plot.endswith('.png')]

    entries = [create_entry(ra[i,0], dec[i,0], value[i,0]) for i in range(count)]
    table = TABLE%dict(rows='\n'.join(entries))
    index = INDEX%dict(table=table)
    with open(filename,'w') as out:
        out.write(index)

################################################################################
print("Creating webpage")
create_index_html(outfile)
