# -*- coding: utf-8 -*-
"""
Sojin SHIM, Ph.D., Instituete of Diagnostic and Interventional Radiology, University Hospital Zurich, 2019-2023
Generate 3 representative images of segmentation in 3 planes.

Input segmentation code
0: air / 1: fat / 2: glands / 3: skin / 4: muscle / 5: bone / 6: fold / 7: implant / 8: suspected region / 9: cylinder
"""

import os
import csv
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as pltimg
from scipy.ndimage import measurements
from scipy import interpolate

now = datetime.datetime.now()
print('Run at ',now.year, now.month, now.day, now.hour, now.minute, now.second)
    
input_path = "D:/test_batch"
output_path = input_path
output_glandularity = output_path + "/Glandularity.csv"

#%%
npimgs = os.listdir(input_path)
for pnpimg in npimgs:
    if 'segmentation_'+pnpimg in npimgs:    
        uid = pnpimg[:-4]
        input_img = input_path + '/'+pnpimg
        input_seg = output_path + '/segmentation_'+pnpimg
        
        img = np.load(input_img).astype(np.int16)
        mask = np.load(input_seg).astype(np.uint8)
        L,R,C = np.shape(img)
        
        thSi = 200
        offset = -2**10
        
        # Density Calculation
        HUGland = -41 # SD=10
        HUFat = -255 # SD=15
        
        HU_G_Muscle = Seg_G_Muscle = 0
        HU_G_Fold = Seg_G_Fold = 0
        
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return array[idx]
        
        def HUcalibration(HU) :
            HUmean=[-469, -255, -41, 173] 
            glandularity = [-1, 0, 1, 2]
            glandularity_near = interpolate.interp1d(HUmean,glandularity)
            
            HUint = np.linspace(-300, 0, num=3001, endpoint=True) # optimise for the glandularity in 3 digit
            HUnear = find_nearest(HUint,HU)
            
            return glandularity_near(HUnear)
        
        Fat = mask==1
        Gland = mask==2
        Skin = mask==3
        Muscle = mask==4
        Fold = mask==6
        breast_tissue = Fat|Gland
        
        C = np.sum(breast_tissue)
        
        HU_tissue = img[breast_tissue]
        meanHU = np.mean(HU_tissue)
        HU_G = HUcalibration(meanHU)
        
        Seg_G = np.sum(Gland)/C
        
        # CSV file for time record.
        fnames = ['uid', 'HU_G', 'Seg_G']
            # first time only: header
        if not os.path.isfile(output_glandularity):
            f = open(output_glandularity, 'w')
            with f:       
                writer = csv.DictWriter(f, fieldnames=fnames)    
                writer.writeheader()    
                
        with open(output_glandularity, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fnames)
            writer.writerow({fnames[0]:uid, fnames[1]:"%.4f" %HU_G, fnames[2]:"%.4f" %Seg_G})
        
        breast = mask!=0
        
        for l in range(50):
            bL = L-l
            if np.sum(breast[bL-1]):
                break
        
        tip = measurements.center_of_mass(breast[bL-1])
        tip = np.array(tip,dtype='int') #  z,x,y
        
        # save the image
        s=0
        pltimg.imsave(output_path+'/img_c_'+uid+'.png', img[s], cmap='gray', vmin=offset, vmax=thSi)
        pltimg.imsave(output_path+'/img_s_'+uid+'.png', img[:,:,tip[1]], cmap='gray',vmin=offset, vmax=thSi)
        pltimg.imsave(output_path+'/img_t_'+uid+'.png', img[:,tip[0],:], cmap='gray',vmin=offset, vmax=thSi)
        
        # save the mask
        s=0
        pltimg.imsave(output_path+'/auto_seg_c_'+uid+'.png', mask[s], cmap='gray',vmin=0, vmax=6)
        pltimg.imsave(output_path+'/auto_seg_s_'+uid+'.png', mask[:,:,tip[1]], cmap='gray',vmin=0, vmax=6)
        pltimg.imsave(output_path+'/auto_seg_t_'+uid+'.png', mask[:,tip[0],:], cmap='gray',vmin=0, vmax=6)
