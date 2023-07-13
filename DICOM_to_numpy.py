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


#%% Input
path_uids = 'D:/test_batch'
output_path = path_uids # default

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
    # dummy slices due to reconstruction restriction
    dummy = 3 
    # remove the cylinder: If you want to include the cylinder in the numpy images, this part has to be moved to segmentation script for cylinder segmentation.
    temp = img_org
    r = C/2-30;
    for i, j in itertools.product(range(0,R),range(0,C)): # nested loops
        if r < np.sqrt(((i+1)-R/2)**2+((j+1)-C/2)**2):
                temp[:,i,j] = -2**10    
    # save slices only upto the length(index) of the breast
    ob = temp > -500
    ob = np.mean(ob,axis=(1,2))
    l = np.sum(ob!=0)
    img=temp[dummy:l]    
    return img

def convert_DICOM_to_Numpy(dir_dicom): 
    dicomimg = load_scan(dir_dicom)
    imgorg = get_pixels_hu(dicomimg) 
    R = imgorg.shape[1]
    C = imgorg.shape[2]
    npimg = img_process(imgorg,R,C)
    return npimg

#%%

uids = os.listdir(path_uids)
uids.sort()

for uid in uids:
    path_DICOM = os.path.join(path_uids,uid)
    
    if os.path.isdir(path_DICOM):
        output_npimg = path_DICOM+'.npy'
        
        dicomimg = load_scan(path_DICOM)
        img = get_pixels_hu(dicomimg)
        R = img.shape[1]
        C = img.shape[2]
        
        npimg = img_process(img,R,C)
        np.save(output_npimg, npimg)
        
        # generate cylinder for future image conversion to DICOM
        
        if 'cylinder.npy' not in uids:
            temp = np.zeros((R,C),dtype=bool)
            
            r = C/2-30;
            for i, j in itertools.product(range(0,R),range(0,C)): # nested loops
                if r < np.sqrt(((i+1)-R/2)**2+((j+1)-C/2)**2):
                        temp[i,j] = True   
                        
            cylinder = temp&(img[3]>-500)    
            np.save(path_uids+'/cylinder.npy',cylinder)
        

