# -*- coding: utf-8 -*-
"""
Sojin SHIM, M.Sc., University Hospital Zurich, 2019-2021
Segmentation of the breast CT image in HU

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
    
HD_path = "D:"

year = '18' # 18/19/20
season = 'FW' # SS/FW
#nseason = 0 # 0/1/3/4 : SS1/SS2/SS3/SS4
num = 1 # 1/2

input_path = HD_path+'/numpy_20'+year+'_'+season+'_'+str(num)
output_path = input_path
output_glandularity = output_path + "_Glandularity_comparison.csv"

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
        C_Skin = np.sum(Skin)
        C_Muscle = np.sum(Muscle)
        C_Fold = np.sum(Fold)
        
        HU_tissue = img[breast_tissue]
        meanHU = np.mean(HU_tissue)
        HU_G = HUcalibration(meanHU)
        
        HU_temp = img[breast_tissue|Skin]
        meanHU = np.mean(HU_temp)
        HU_G_Skin = HUcalibration(meanHU)
        
        Seg_G = np.sum(Gland)/C
        Seg_G_Skin = (np.sum(Gland)+C_Skin)/(C+C_Skin)
        
        HU_G_Muscle = 0
        Seg_G_Muscle = 0
        HU_G_Fold = 0
        Seg_G_Fold = 0
        
        if C_Muscle:
            HU_temp = img[(breast_tissue)|Muscle]
            meanHU = np.mean(HU_temp)
            HU_G_Muscle = HUcalibration(meanHU)
            Seg_G_Muscle = (np.sum(Gland)+C_Muscle)/(C+C_Muscle)
            
        if C_Fold:    
            HU_temp = img[(breast_tissue)|Fold]
            meanHU = np.mean(HU_temp)
            HU_G_Fold = HUcalibration(meanHU)
            Seg_G_Fold = np.sum(Gland)/(C+C_Fold)
        
        # CSV file for time record.
        fnames = ['id', 'HU_G', 'Seg_G', 'HU_G_w_Skin', 'Seg_G_w_Skin', 
                  'HU_G_w_Muscle', 'Seg_G_w_Muscle', 'HU_G_w_Fold', 'Seg_G_w_Fold',
                  'C_Breast','C_Skin','C_Muscle','C_Fold']
            # first time only: header
        if not os.path.isfile(output_glandularity):
            f = open(output_glandularity, 'w')
            with f:       
                writer = csv.DictWriter(f, fieldnames=fnames)    
                writer.writeheader()    
                
        with open(output_glandularity, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fnames)
            writer.writerow({fnames[0]:uid, fnames[1]:"%.4f" %HU_G, fnames[2]:"%.4f" %Seg_G, 
                             fnames[3]:"%.4f" %HU_G_Skin, fnames[4]:"%.4f" %Seg_G_Skin, 
                             fnames[5]:"%.4f" %HU_G_Muscle, fnames[6]:"%.4f" %Seg_G_Muscle, 
                             fnames[7]:"%.4f" %HU_G_Fold, fnames[8]:"%.4f" %Seg_G_Fold,
                             fnames[9]:"%.0f" %C,fnames[10]:"%.0f" %C_Skin,
                             fnames[11]:"%.0f" %C_Muscle,fnames[12]:"%.0f" %C_Fold,})
        
        breast = mask!=0
        
        for l in range(50):
            bL = L-l
            if np.sum(breast[bL-1]):
                break
        
        tip = measurements.center_of_mass(breast[bL-1])
        tip = np.array(tip,dtype='int') #  z,x,y
        
        # save the image
        s=0
        pltimg.imsave(output_path+'/img_%d_c.png' % int(uid),img[s], cmap='gray',vmin=offset, vmax=thSi)
#        pltimg.imsave(output_path+'/img_%d_c.tiff' % int(uid),img[s], cmap='gray',vmin=offset, vmax=thSi, dpi=600, format="tiff")
        pltimg.imsave(output_path+'/img_%d_s.png' % int(uid),img[:,:,tip[1]], cmap='gray',vmin=offset, vmax=thSi)
#        pltimg.imsave(output_path+'/img_%d_s.tiff' % int(uid),img[:,:,tip[1]], cmap='gray',vmin=offset, vmax=thSi, dpi=600, format="tiff")
        pltimg.imsave(output_path+'/img_%d_t.png' % int(uid),img[:,tip[0],:], cmap='gray',vmin=offset, vmax=thSi)
#        pltimg.imsave(output_path+'/img_%d_t.tiff' % int(uid),img[:,tip[0],:], cmap='gray',vmin=offset, vmax=thSi, dpi=600, format="tiff")
        
        # save the mask
        s=0
        pltimg.imsave(output_path+'/auto_seg_%d_c.png' % int(uid),mask[s], cmap='gray',vmin=0, vmax=6)
#        pltimg.imsave(output_path+'/auto_seg_%d_c.tiff' % int(uid),mask[s], cmap='gray',vmin=0, vmax=6, dpi=600, format="tiff")
        pltimg.imsave(output_path+'/auto_seg_%d_s.png' % int(uid),mask[:,:,tip[1]], cmap='gray',vmin=0, vmax=6)
#        pltimg.imsave(output_path+'/auto_seg_%d_s.tiff' % int(uid),mask[:,:,tip[1]], cmap='gray',vmin=0, vmax=6, dpi=600, format="tiff")
        pltimg.imsave(output_path+'/auto_seg_%d_t.png' % int(uid),mask[:,tip[0],:], cmap='gray',vmin=0, vmax=6)
#        pltimg.imsave(output_path+'/auto_seg_%d_t.tiff' % int(uid),mask[:,tip[0],:], cmap='gray',vmin=0, vmax=6, dpi=600, format="tiff")