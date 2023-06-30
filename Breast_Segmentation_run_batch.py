"""
Sojin SHIM, M.Sc., University Hospital Zurich, 2019-2021
Segmentation of the breast CT image in HU

From numpy HU array to numpy segmentation array
The numpy array is already read in the preprocessing
Input numpy: HU in [-2**10, 2*10] without the cylinder in int16
Output numpy: Segmentation in uint8
No dicom images are handled in this script
Breast popluation study in another script

Output segmentation code
0: air / 1: fat / 2: glands / 3: skin / 4: muscle / 5: bone / 6: fold / 7: implant / 8: suspected region / 9: cylinder
Output segmentation time
total, rib, silicone implant, muscle, skin, gland
"""
#%%
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import function_Breast_Segmentation # customized

#%% Input

HD_path = "D:"

year = '18' # 18/19/20
season = 'FW' # SS/FW
#nseason = 2 # 0/1/2/4 : SS1/SS2/SS3/SS4
num = 1 # 1/2

input_path = HD_path+'/numpy_20'+year+'_'+season+'_'+str(num)
output_path = input_path

re_fold_2018FW1 = [1810011,1810060,1810070,1810140,1810141,1810150,1810151,
                   1810201,1810221,1810251,1810280,1810331,1810350,1810351,
                   1810451,1810461,1810490,1810510,1810511,1810620,1810621,
                   1810641,1810680,1810681,1810701,1810710,1810711,1810741,
                   1810790,1810801,1810810,1810811,1810910,1810971,1811060,
                   1811171,1811211,1811220,1811221,1811250,1811251,1811330,
                   1811401,1811440,1811441,1811530,1811531,1811540,1811580,1811581]
#re_fold_2019SS1 = [1900011,1900021,1900231,1900280]
#%%
npimgs = os.listdir(input_path)

for pnpimg in npimgs:
    if ('segmentation_'+pnpimg not in npimgs)&(len(pnpimg)==11):
    #if 1: #for rerun
        uid = pnpimg[:-4]
        if 1:
        #if int(uid)-10000 in re_fold_2018FW1:# rerun for the fold segmentation modification
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
            
            
#%% read csv and rerun the segmentation for specific ones: ex) ones with tFold
# 2018FW all, 2019SS1 until uid <= 32
csv_time = output_path + "Segmentation_time.csv"

tseg = pd.read_csv(csv_time)
