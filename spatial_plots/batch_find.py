import os
import argparse
import numpy as np

jobs = 10

# submit many different cuts at once
def submit_finding(mag_cut):
        str_mag_cut = str(int(10*mag_cut))
        logfile = 'log_dir/finding_results_cut={}.log'.format(str_mag_cut)
        batch = 'csub -n {} -o {} '.format(jobs, logfile) # testing condor updates
        command = 'python find_overdensities.py {:0.2f}'.format(mag_cut)
        command_queue = batch + command

        print(command_queue)
        os.system(command_queue) # Submit to queue
        return

mags = [22.5, 22.0, 21.5, 21.0, 20.5]

for i in mags:
        submit_finding(i)
