"""
Sojin SHIM, Ph.D., University Hospital Zurich, 2019-2023
Segmentation of the breast CT image in HU

From numpy HU array to numpy segmentation array
The numpy array is already read in the preprocessing
Input numpy: HU in [-2**10, 2*10] without the cylinder in int16
Output numpy: Segmentation in uint8
No dicom images are handled in this script
Breast popluation study in another script

Output segmentation code
0: air / 1: fat / 2: glands / 3: skin / 4: muscle / 5: bone / 6: fold / 7: implant / 8: suspected region / 9: cylinder
"""
#%%
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import function_Breast_Segmentation # customized

#%% Inputs
input_path = "D:/test_batch"
output_path = input_path # default

#%% Segment and save output numpy files
npimgs = os.listdir(input_path)

for pnpimg in npimgs:
    # for numpy images that were not previously segmented excluding cylinder.npy
    if (pnpimg[-4:]=='.npy')&(pnpimg[:8]!='cylinder')&(pnpimg[:13]!='segmentation_')&('segmentation_'+pnpimg not in npimgs):
        uid = pnpimg[:-4]
        img_path = os.path.join(input_path,pnpimg)
        output_seg = output_path+'/segmentation_'+pnpimg
        output_fold = output_path+'/segmentation_fold_'+pnpimg
        
        img = np.load(img_path).astype(np.int16)
        mask, fold = function_Breast_Segmentation.Breast_Segmentation(uid, img, output_path)
        np.save(output_seg, mask)
        
        if np.sum(fold):
            temp = np.sum(fold,axis=1)
            temp = np.sum(temp,axis=1)
            ttemp = temp>0
            fl=np.sum(ttemp)
            np.save(output_fold,fold[:fl])
        
            
#%% Record segmentation time
csv_time = output_path + "/Segmentation_time.csv"

tseg = pd.read_csv(csv_time)
