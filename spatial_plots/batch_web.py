import os
import argparse
import numpy as np

jobs = 10

def submit_webpage(mag_cut):
    str_mag_cut = str(int(10*mag_cut))

    logfile = 'log_dir/webpage_results_cut={}.log'.format(str_mag_cut)
        batch = 'csub -n {} -o {} --host all '.format(jobs, logfile) # testing condor updates
        command = 'python image_webpage.py {:0.2f}'.format(mag_cut)
        command_queue = batch + command

        print(command_queue)
        os.system(command_queue) # Submit to queue
        return

mags = [22.5, 22.0, 21.5, 21.0, 20.5]

for i in mags:
    submit_webpage(i)
