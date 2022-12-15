#!/usr/bin/env python
'''
MDM 1.3m Image reduction Pipeline

Author: Bhagya M Subrayan

v. April 4  2022

This script is written for reducing follow-up photometric observations for REFITT
taken from MDM 1.3m Telescopes with instruments Templeton and Echelle using latest version of astropy.

The output is: a directory with reduced images for either of the specified instruments
for the monthly observing run. Runs in astroconda environment.

The combined flats is made using all the flats of the entire observing run. If a date is specified that the script combines flats taken only on th specified date. Please make sure that there are relevant flats if a date is invoked in the command line.

Sample run if using flats from the entire run: 

python mdm_ccdproc.py -r jan21 -temp -bpc 139 140

Sample run if you want to use only the flats on a given specified date: 

python mdm_ccdproc.py -r jan21 -d 20210117 -temp -bpc 139 140

-r: Folder name correponding to a observing run
-d: Night of Observation
-temp: CCD used. -eche if echelle
-bpc: Any badpixel column that needs to be fixed (in ds9 cooredinates)
-lacosmicoff : If no lacosmic needs to be applied
'''
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
import sys
import argparse
from astropy.nddata import CCDData
from astropy import units as u
import ccdproc as ccdp
import shutil
from scipy.stats import mode
from astropy.stats import mad_std
from glob import glob


def parse_args():
    '''Parse command line arguments'''
    parser= argparse.ArgumentParser()
    parser.add_argument('-r','--run',default =None,help='Name of the directory of the observring run in the format monyy (jan21)')
    parser.add_argument('-d','--date',default =None,help='Date in UTC of the observring night in the format 210209. If not specified all sub_directories processed')
    parser.add_argument('-temp','--templeton',action = 'store_true',help='Templeton')
    parser.add_argument('-eche','--echelle',action = 'store_true',help='Echelle')
    parser.add_argument('-lacosmicoff','--no_cosmic',action='store_true', help='No la_cosmic correction will be applied if this argument is specified')
    parser.add_argument('-bpc','--badpixelcol',default =None,nargs = '+',type=str,help='Any bad pixel column in ds9 coordinates')
    args = parser.parse_args()
    return args

args = parse_args()

if args.templeton:
	inst = 'templeton'
elif args.echelle:
	inst ='echelle'
else:
	raise Exception('Templeton (-tmp) or Echelle (-eche) must be selected.')

if args.no_cosmic:
    cosmic = 0
    print('Warning: LA_Cosmic turned off')
else:
    cosmic = 1

print('\n')
print('************ CCD Reduction using Astropy ccdproc *****************')
print('\n')
print('                      Author: Bhagya. M Subrayan                             ')
print('\n')

'''Defining functions that will be helpful in reduction'''

'''Bad column'''
if args.badpixelcol is not None:
	badpixcol = [int(item) for item in args.badpixelcol]



def bad_column(image,badpix):
    for i in range(len(badpix)):
        j = badpix[i] - 1 #139 140 for Templeton
        new_col = (image.data[:,j-4]+ image.data[:,j+4])/2
        #print(new_col)
        for k in range(np.shape(image)[0]):
            image.data[k,j] = new_col[k]
    return image


'''For removing the original raw data copied from the mdm13m/raw directory.
The remaining will only contain processed .fits file'''
def remove_raw(dir):
	list = os.listdir(dir)
	for file in list:
		if file.endswith('.fit'):
			os.remove(dir + file)

'''Used for scaling when combining flats'''
def inv_median(a):
    return 1 / np.median(a)

'''Defining paths. Making a /reduced/"observing_run" directory if not already existing'''
basepath = './' + args.run +'/'
reduced_path = '../reduced/' + args.run +'/'
r = Path(reduced_path)
if not r.exists():
	r.mkdir()

'''Path to the raw directory'''
#raw_data = Path('./' + args.run)
subd_list = np.sort(os.listdir(basepath))
if 'bias' in subd_list:
    nights =  np.delete(subd_list,np.where(subd_list == 'bias'))
else:
    nights = subd_list
            
'''Making sure all the directories in the nights array does have raw files in them'''	
for item in nights:
    if (len(os.listdir(basepath+item)) == 0):
        print('Empty directory found on day:', item)
        nights = np.delete(nights,np.where(nights == item))

print('Observing Nights listed : ', nights)


if args.date is None:
    '''Filtering flats and science on all the days from the sub-directories. flat_list has useful flats, reject has rejected flats'''
    '''Useful flats available in all subdirectories of the run retained for combining flats'''

    flat_path = reduced_path + 'flats/'
    p = Path(flat_path)
    if not p.exists():
            p.mkdir()

    flat_list = []
    reject = []
    for day in nights:
            ifc = ccdp.ImageFileCollection(basepath + day)
            flat_files= ifc.files_filtered(imagetyp= 'FLAT')
            if (day == 'bias' or day =='flats'):
                    continue
            elif(len(flat_files) == 0):
                    print('No flat fields taken on :', day)
                    continue
            else:
                    for file in flat_files:
                            if (15000 < (ccdp.CCDData.read(basepath + day +'/'+file,unit='count').data.mean()) < 35000):
                                    flat_list.append(file)
                                    try:
                                            shutil.copy(basepath+day+'/'+file, flat_path)
                                    except shutil.SameFileError:
                                            pass
                            else:
                                    reject.append(file)

    '''Overscan, trim and fixpix flats'''
    ifc_flat = ccdp.ImageFileCollection(flat_path)
    for ccd, file_name in ifc_flat.ccds(imagetyp='FLAT',ccd_kwargs={'unit': 'adu'},return_fname=True):
            if(inst == 'templeton'):
                    ccd = bad_column(ccd,badpix = badpixcol)
                    ccd = ccdp.subtract_overscan(ccd,median=True, overscan=ccd[:,532:541],overscan_axis=1)
                    ccd = ccdp.trim_image(ccd[:,18:529])
            else:
                    ccd = ccdp.subtract_overscan(ccd,median=True, overscan=ccd[:,995:1055],overscan_axis=1)
                    ccd = ccdp.trim_image(ccd[:,0:993])
            new_name = file_name.split('.')[1] +'_'+ str(ccd.header['filter']) +'_flat.fits'
            ccd.write( flat_path + new_name)

    remove_raw(flat_path)
    print('\n')
    print('All calibration and science images obtained')
    print('Instrument used: ' , inst)
    if (inst == 'templeton'):
        print('Bad column in the CCD chip fixed')
    print('\n')

    ifc_flat = ccdp.ImageFileCollection(flat_path)
    flat_filters = set(h['filter'] for h in ifc_flat.headers(imagetyp='FLAT'))
    print('Flat fields of the following filters where taken during the run:', flat_filters)

    for filt in flat_filters:
            filt = str(filt)
            to_combine = ifc_flat.files_filtered(imagetyp='FLAT', filter= filt, include_path=True)
            print(to_combine)
            combined_flat = ccdp.combine(to_combine,method='median', scale=inv_median,sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,mem_limit=350e6)
            combined_flat.meta['combined'] = True
            com_name = 'combined_flat_filter_{}.fits'.format(filt.replace("''", "p"))
            combined_flat.write( flat_path + com_name)

    ifc_flat = ccdp.ImageFileCollection(flat_path)
    flats_set = {ccd.header['filter']: ccd for ccd in ifc_flat.ccds(imagetyp='FLAT', combined=True)}
    #print(flats_set)
    print('All caliberation files are overscan subtracted and trimmed')
    print('Combining flats in all the filters complete. \n')
    print('The combined flats can be checked under: mdm13m/reduced/"mmyy"/flats/')

    '''Reducing Science Images'''

    print('\n')
    print('Processing the science files....')


    for day in nights:
        if (day == 'bias' or day =='flats'):
            continue
        else:
            print('Reducing for the night : ', day)
            ifc = ccdp.ImageFileCollection(basepath + day)
            science_files = ifc.files_filtered(imagetyp= 'OBJECT')
            #print(science_files)
            obj_path = reduced_path + day+'/' #these needs to be processed just sorting now
            o_p = Path(obj_path)
            if not o_p.exists():
                o_p.mkdir()
            [shutil.copy(basepath+day+'/'+new_file,obj_path) for new_file in science_files]
            for light, file_name in ifc.ccds(imagetyp='OBJECT', return_fname=True, ccd_kwargs=dict(unit='adu')):
                if(light.header['filter'] not in flat_filters):
                    continue
                elif(inst == 'templeton'):
                    reduced = bad_column(light,badpix = badpixcol)
                    reduced = ccdp.subtract_overscan(reduced,median=True, overscan=reduced[:,532:541],overscan_axis=1)
                    reduced = ccdp.trim_image(reduced[:,18:529])
                else:
                    reduced = ccdp.subtract_overscan(light,median=True, overscan=light[:,995:1055],overscan_axis=1)
                    reduced = ccdp.trim_image(reduced[:,0:993])
                good_flat = flats_set[reduced.header['filter']]
                new_name = light.header['object']+'_mdm13m_'+ str(light.header['filter']) +'_jd_'+str(light.header['JD'])
                reduced = ccdp.flat_correct(reduced, good_flat)
                if(cosmic == 1):
                    reduced= ccdp.cosmicray_lacosmic(reduced,sigclip=8.0,sigfrac=0.75,readnoise = 5.30,objlim=10,gain = 3.47,satlevel=65535,gain_apply=True)
                    mask = CCDData(data = reduced.mask.astype('uint8'),unit = u.dimensionless_unscaled)
                    mask.header['imagetyp'] = 'la_cosmic_mask'
                    mask_name= new_name+'_mask.fits'
                    #mask.write(o_p / mask_name)
                #else:
                    #print('No la_cosmic correction made')
                reduced.header['RAW-NUM'] = file_name.split('.')[1]
                reduced.write(obj_path + new_name+ '.fits')
            remove_raw(obj_path)
else:
    date = args.date
    #print(date)
    flat_list = []
    reject = []
    flat_path = reduced_path +date+ '_flats/'
    p = Path(flat_path)
    if not p.exists():
            p.mkdir()

    ifc = ccdp.ImageFileCollection(basepath + date)
    flat_files= ifc.files_filtered(imagetyp= 'FLAT')
    for file in flat_files:
        if (15000 < (ccdp.CCDData.read(basepath + date +'/'+file,unit='count').data.mean()) < 35000):
            flat_list.append(file)
            try:
                shutil.copy(basepath+date+'/'+file, flat_path)
            except shutil.SameFileError:
                pass
        else:
            reject.append(file)


        '''Overscan, trim and fixpix flats'''
    ifc_flat = ccdp.ImageFileCollection(flat_path)
    for ccd, file_name in ifc_flat.ccds(imagetyp='FLAT',ccd_kwargs={'unit': 'adu'},return_fname=True):
            if(inst == 'templeton'):
                    ccd = bad_column(ccd,badpix = badpixcol)
                    ccd = ccdp.subtract_overscan(ccd,median=True, overscan=ccd[:,532:541],overscan_axis=1)
                    ccd = ccdp.trim_image(ccd[:,18:529])
            else:
                    ccd = ccdp.subtract_overscan(ccd,median=True, overscan=ccd[:,995:1055],overscan_axis=1)
                    ccd = ccdp.trim_image(ccd[:,0:993])
            new_name = file_name.split('.')[1] +'_'+ str(ccd.header['filter']) +'_flat.fits'
            ccd.write( flat_path + new_name)

    remove_raw(flat_path)
    print('\n')
    print('All calibration and science images obtained')
    print('Instrument used: ' , inst)
    if (inst == 'templeton'):
        print('Bad column in the CCD chip fixed')
    print('\n')

    ifc_flat = ccdp.ImageFileCollection(flat_path)
    flat_filters = set(h['filter'] for h in ifc_flat.headers(imagetyp='FLAT'))
    print('Flat fields of the following filters where taken during the run:', flat_filters)

    
    for filt in flat_filters:
            filt = str(filt)
            to_combine = ifc_flat.files_filtered(imagetyp='FLAT', filter= filt, include_path=True)
            print(to_combine)
            combined_flat = ccdp.combine(to_combine,method='median', scale=inv_median,sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,mem_limit=350e6)
            combined_flat.meta['combined'] = True
            com_name = 'combined_flat_filter_{}.fits'.format(filt.replace("''", "p"))
            combined_flat.write( flat_path + com_name)

    ifc_flat = ccdp.ImageFileCollection(flat_path)
    flats_set = {ccd.header['filter']: ccd for ccd in ifc_flat.ccds(imagetyp='FLAT', combined=True)}
    #print(flats_set)
    print('All caliberation files are overscan subtracted and trimmed')
    print('Combining flats in all the filters complete. \n')
    print('The combined flats can be checked under: mdm13m/reduced/"mmyy"/'+date+'_flats')

    '''Reducing Science Images'''

    print('\n')
    print('Processing the science files....')


    print('Reducing for the night : ', date)
    ifc = ccdp.ImageFileCollection(basepath + date)
    science_files = ifc.files_filtered(imagetyp= 'OBJECT')
    obj_path = reduced_path + date+'/' #these needs to be processed just sorting now
    o_p = Path(obj_path)
    if not o_p.exists():
        o_p.mkdir()
    [shutil.copy(basepath+date+'/'+new_file,obj_path) for new_file in science_files] 
    for light, file_name in ifc.ccds(imagetyp='OBJECT', return_fname=True, ccd_kwargs=dict(unit='adu')):
        if (light.header['filter'] not in flat_filters):
            continue
        elif (inst == 'templeton'):
            reduced = bad_column(light,badpix = badpixcol)
            reduced = ccdp.subtract_overscan(light,median=True, overscan=light[:,532:541],overscan_axis=1)
            reduced = ccdp.trim_image(reduced[:,18:529])
        else:
            reduced = ccdp.subtract_overscan(light,median=True, overscan=light[:,995:1055],overscan_axis=1)
            reduced = ccdp.trim_image(reduced[:,0:993])
        good_flat = flats_set[reduced.header['filter']]
        reduced = ccdp.flat_correct(reduced, good_flat)
        new_name = light.header['object']+'_mdm13m_'+ str(light.header['filter']) +'_jd_'+str(light.header['JD'])
        if(cosmic ==1):
            reduced = ccdp.cosmicray_lacosmic(reduced,sigclip=8.0,sigfrac=0.75,readnoise = 5.30,objlim=10,gain = 3.47,satlevel=65535,gain_apply=True)
            mask = CCDData(data = reduced.mask.astype('uint8'),unit = u.dimensionless_unscaled)
            mask.header['imagetyp'] = 'la_cosmic_mask'
            mask_name = new_name +'_mask.fits'
            #mask.write(o_p / mask_name)
        #else:
            #print('No la_cosmic correction made')
        reduced.header['RAW-NUM'] = file_name.split('.')[1]
        reduced.write(obj_path + new_name+'.fits')

    remove_raw(obj_path)


'''Sorting all the reduced files and creating an log of all the files in the respective directories. Headers with date and reducer edited'''
from datetime import date
today = date.today()

from astropy.io import ascii
from astropy.io import fits
sort_dir = os.listdir(reduced_path)
for folder in sort_dir:
    if os.path.isfile(reduced_path+folder+ '/file_log.txt'):
        print(folder)
        continue
    elif(folder == 'galex_rob'):
        continue
    else:
        #print(folder)
        for item in os.listdir(reduced_path+folder):
            #print(item)
            data, hdr = fits.getdata(reduced_path + folder+'/'+ item,header=True)
            hdr['NAME-RED'] = ('BSUBRAYA','Name of the reducer')
            hdr['DATE-RED'] = (today.strftime('%b-%d-%Y'),'Date of reduction') 
            fits.writeto(reduced_path + folder+'/'+ item,data,hdr,overwrite= True)
        ifc_fold = ccdp.ImageFileCollection( reduced_path + folder)
        table = ifc_fold.summary['file', 'imagetyp','filter','exptime','airmass']
        table.write(reduced_path + folder+'/file_log.txt',format='ascii',overwrite=True)
        if ((folder != 'flats') or (folder !='galex_rob')):
            ztf_id = np.unique(np.sort([h['OBJECT'] for h in ifc_fold.headers(imagetyp='OBJECT',ccd_kwargs = {'unit':'adu'})]))
            for iditem in ztf_id:
                if (iditem != ''):
                    ztf_path = Path(reduced_path + folder +'/' + iditem)
                    if not ztf_path.exists():
                        ztf_path.mkdir()
                    for ztf_file in glob(reduced_path + folder+'/*'+iditem+'*.fits'):
                        shutil.move(ztf_file,reduced_path + folder +'/' + iditem+'/')
                    else:
                        continue
                else:
                    print('A header is missig a ZTFID. Please check the observing log for any anomalies')


            
print('\n')
print('CCD Reduction complete.')

print('You can find the reduced files: ../reduced/"run_name"')
print('\n')
print('*******************************************************************')
print('\n')
