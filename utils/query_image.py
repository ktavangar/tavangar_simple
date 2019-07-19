import urllib
import yaml
import os
from PIL import Image

#with open('config.yaml', 'r') as ymlfile:
#    cfg = yaml.load(ymlfile)

#save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
#if not os.path.exists(save_dir):
#    os.mkdir(save_dir)

def retrieve_image(filename, ra, dec):
    filename = '{}/image_{:0.2f}_{:0.2f}.png'.format('decals_dr7_save', ra, dec)
    url = "http://legacysurvey.org/viewer/jpeg-cutout?ra={0}&dec={1}&zoom={2}&layer=decals-dr7"
    #urllib.request.urlretrieve(url.format(ra, dec, 10), filename) #Retreaves and saves each image                                                  
    urllib.urlretrieve(url.format(ra, dec, 9), filename) #Retreaves and saves each image                                                  
    #img = Image.open(filename)
    #ax.imgshow(img)
    #ax.axis('off')

    return
