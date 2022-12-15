#!/usr/bin/env python
'''
MDM 1.3m Image WCS Taggging

Author: Bhagya M Subrayan

v. April 4 2022
python wcs_script.py -d <date pointing to the directory in reduced>
'''
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
import subprocess
import argparse
import sys
import shutil
from scipy.stats import mode
from astropy.stats import mad_std
import glob


def parse_args():
	'''Parse command line arguments'''
	parser= argparse.ArgumentParser()
	parser.add_argument('-r','--run',default =None,help='Name of the directory of the observring run in the format monyy (jan21)')
	parser.add_argument('-d','--date',default =None,help='Date in UTC of the observring night in the format 210209. If not specified all sub_directories processed')
	args = parser.parse_args()
	return args

args = parse_args()

reduced_path = '../reduced/' + args.run + '/' + args.date +'/'
ztf_dir = np.array(os.listdir(reduced_path))
ztf_dir = np.delete(ztf_dir,np.where(ztf_dir == 'file_log.txt'))


#print(ztf_list)

final_path = '/project/amalthea/refitt/photometry/mdm_final/'

for item in ztf_dir:
    fits_files = os.listdir(reduced_path+ item + '/')
    for ind_file in fits_files:
        cmd = ['solve-field',ind_file]
        subprocess.run(cmd, cwd = reduced_path+ item+'/')
    wcs_path = Path(reduced_path+ item+'/wcs_reduced')
    if not wcs_path.exists():
        wcs_path.mkdir()
    new_files = glob.glob(reduced_path+ item + '/*.new')
    [shutil.copy(plate_file, reduced_path+ item+'/wcs_reduced/') for plate_file in new_files]
    final_dest = Path(final_path+item)
    if not final_dest.exists():
        final_dest.mkdir()
    for filename in os.listdir(reduced_path+ item+'/wcs_reduced'):
        split = filename.split('.')
        new_name = split[0] +'.' +split[1]+'.fits'
        os.rename(reduced_path+ item+'/wcs_reduced/'+filename,reduced_path+ item+'/wcs_reduced/'+new_name)
    #final_files = glob.glob(reduced_path+item+'/wcs_reduced/'+'/*.fits')
    #print('Copying to final destination')
    #[shutil.copy(plate_file,final_path+item+'/')for plate_file in final_files]


    





    

