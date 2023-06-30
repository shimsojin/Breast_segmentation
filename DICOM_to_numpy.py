# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 18:53:13 2020

@author: Sojin Shim
"""

import os
import csv
import pydicom as dicom
import numpy as np
import itertools


#%% input
HD = 'D:'
year = '19' # 18/19/20
season = 'SS' # SS/FW
nseason = 0
num = '2' # 1/2

dirs_sel = ['Stand-L','Stand-R']
#dirs_sel = ['Stand-L','Stand-R','Stand-L-prior','Stand-R_prior']

uidcsv = HD+year+'_uid_info.csv'
fnames = ['season','dir_number','uid','dir_list','converted dir','length']

#%%
def load_scan(path):
    slices = [dicom.read_file(path + '/' + s) for s in os.listdir(path)]
    slices.sort(key = lambda x: int(x.InstanceNumber))
    return slices

def get_pixels_hu(scans):
    image = np.stack([s.pixel_array for s in scans])
    # Convert to int16 (from sometimes int16), 
    # should be possible as values should always be low enough (<32k)
    image = image.astype(np.int16)
    # Set outside-of-scan pixels to 1
    # The intercept is usually -1024, so air is approximately 0
    image[image == -2000] = 0
    # Convert to Hounsfield units (HU)
    intercept = scans[0].RescaleIntercept
    slope = scans[0].RescaleSlope
    if slope != 1:
        image = slope * image.astype(np.float64)
        image = image.astype(np.int16)   
    image += np.int16(intercept)
    return np.array(image, dtype=np.int16)

def img_process(img_org,R,C):
    default = 3 # dummy slices due to reconstruction restriction
    # remove the ring (return to air)
    temp = img_org
    r = C/2-30;
    for i, j in itertools.product(range(0,R),range(0,C)): # nested loops
    #for i, j in zip(range(x), range(y)): #simultaneous loops
        if r < np.sqrt(((i+1)-R/2)**2+((j+1)-C/2)**2):
                temp[:,i,j] = -2**10    
    # find the length(index) of the breast
    ob = temp > -500
    ob = np.mean(ob,axis=(1,2))
    l = np.sum(ob!=0)
    img=temp[default:l]    
    return img

def convert_DICOM_to_Numpy(dir_dicom): 
    dicomimg = load_scan(dir_dicom)
    imgorg = get_pixels_hu(dicomimg) 
    R = imgorg.shape[1]
    C = imgorg.shape[2]
    npimg = img_process(imgorg,R,C)
    return npimg


#%%
dir_DICOM = 'DICOM_retrieve_20'+year+'_'+season+'_anony_'+num
path_uids = os.path.join(HD,dir_DICOM)
uids = os.listdir(path_uids)
uids.sort()

output_path = path_uids

for uid in uids:
    dirlist = [] # list of dir in the uid
    blist = [] # converted
    lengths = []
    path_dirs = os.path.join(path_uids,uid)
    dirs = os.listdir(path_dirs)
    
    dlist = []
    for d in dirs:
        if os.path.isdir(os.path.join(path_dirs,d)):
            dirlist.append(d)
    
    for n,d in enumerate(dirs_sel):
        path_DICOM = os.path.join(path_dirs,d)
        if os.path.isdir(path_DICOM):
            blist.append(d)
            output_npimg = output_path+str(year)+str(nseason)+str(uid)+str(n)+'.npy' 
            
            npimg = convert_DICOM_to_Numpy(path_DICOM) # load DICOM and get numpy
            L = npimg.shape[0]
            lengths.append(L)
            
            np.save(output_npimg, npimg)
    length = np.sum(lengths)
                    
    #%% CSV file save
    
    seperator = ','
    dl = seperator.join(dirlist)
    bl = seperator.join(blist)
    
    if not os.path.isfile(uidcsv):
        f = open(uidcsv, 'w')
        with f:       
            writer = csv.DictWriter(f, fieldnames=fnames)    
            writer.writeheader()    
            
    with open(uidcsv, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fnames)
        writer.writerow({fnames[0]:season, fnames[1]:num, fnames[2]:str(uid), 
                         fnames[3]:dl, fnames[4]:bl, fnames[5]:length})

