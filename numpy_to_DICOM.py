# -*- coding: utf-8 -*-
"""
Sojin SHIM, M.Sc., University Hospital Zurich, 2019-2021
Dose analysis based on the MC simulation result

input: the outcome raw file of the MC simulation
output: the dose anlaysis values and 3 plane images
"""
import numpy as np
import pydicom as dicom
import os
import datetime

now = datetime.datetime.now()
print('Run at ',now.year, now.month, now.day, now.hour, now.minute, now.second)

#%%
year = '18' # 18/19/20
season = 'FW' # SS/FW
#nseason = 1 # 0,1 / 2,3
num = '1' # 1/2

dir_numpy = 'numpy_20'+year+'_'+season+'_'+num+''

path = 'D:/'
#path = '/Volumes/SSD_SS/'
#dir_DICOM = 'DICOM_retrieve_20'+year+'_'+season+'_anony'# input DICOM retrieved from PAX
path_dir_numpy = os.path.join(path,dir_numpy)
path_dir_DICOM = os.path.join(path,'InputStorageSegmentation/')

seg_numpy_flag = 'segmentation_'+year

dicom_file = os.path.join(path,'CT00000003.dcm')
cylinder = np.load(os.path.join(path,'cylinder.npy'))
offset = -2**10 # dicom_defualt.RescaleIntercept

instant_num = 3
dx = 0.3

list_npid = os.listdir(path_dir_numpy)

#%%
for npid in list_npid:
    if npid[:15]==seg_numpy_flag:
        DICOMid=npid[13:20]
        DICOM_dir = os.path.join(path_dir_DICOM,DICOMid)
        
        if not os.path.isdir(DICOM_dir):
            ds =  dicom.read_file(dicom_file)
            ds[0x20,0x32].value[1] = -3*dx
            ds[0x20,0x1041].value = 3*dx
            
            seg_img = np.load(path_dir_numpy+'/'+npid)
            L = seg_img.shape[0]

            # add the cylinder matrix on the segmented image
            temp = np.add(seg_img,cylinder*8)-offset 
            img_seg_cylinder = temp.astype(np.uint16) # np.unit16 for DICOM pixel

            os.mkdir(DICOM_dir)
            for i in range(L):
                dcm_name='CT'+'%08d'%i+'.dcm'

                ds.PixelData = img_seg_cylinder[i].tobytes()
                ds.save_as(DICOM_dir+'/'+dcm_name)
                
                ds[0x20,0x13].value += 1
                ds[0x20,0x32].value[1] -= dx
                ds[0x20,0x1041].value += dx
                
                temp = ds[0x08,0x18].value
                temp = temp[:-3]+str(int(temp[-3:])+1)
                ds[0x08,0x18].value = temp
            
            del ds
