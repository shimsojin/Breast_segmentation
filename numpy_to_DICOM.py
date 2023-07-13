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

#%% INPUTS
Input_path = 'D:/test_batch'
Output_path = Input_path # default

seg_numpy_flag = 'segmentation_'

dicom_file = os.path.join(Input_path,'CT00000003.dcm')
cylinder = np.load(os.path.join(Input_path,'cylinder.npy'))
Cylinder = 1 # 1-On /0-Off

#%%
offset = -2**10 # dicom_defualt.RescaleIntercept

instant_num = 3
dx = 0.3

list_npid = os.listdir(Input_path)

#%%
for npid in list_npid:
    if (npid[:13]==seg_numpy_flag)(npid[-4:]=='.npy'):
        DICOMid=npid[:-4]
        DICOM_dir = os.path.join(Output_path,DICOMid)
        
        if not os.path.isdir(DICOM_dir):
            ds =  dicom.read_file(dicom_file)
            ds[0x20,0x32].value[1] = -3*dx
            ds[0x20,0x1041].value = 3*dx
            
            seg_img = np.load(Input_path+'/'+npid)
            L = seg_img.shape[0]

            # add the cylinder matrix on the segmented image
            if Cylinder:
                seg_img = np.add(seg_img,cylinder*8)-offset 
                
            img_seg_cylinder = seg_img.astype(np.uint16) # np.unit16 for DICOM pixel

            os.mkdir(DICOM_dir)
            for i in range(L):
                dcm_name='CT'+'%08d'%i+'.dcm'

                ds.PixelData = img_seg_cylinder[i].tobytes()
                ds.save_as(DICOM_dir+'/'+dcm_name)    
            del ds
