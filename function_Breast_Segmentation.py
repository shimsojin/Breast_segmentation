# -*- coding: utf-8 -*-
"""
Sojin SHIM, Ph.D., Institute of Diagnostic and Interventional Radiology, University Hospital Zurich, 2019-2023
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

# %%
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
import datetime
from time import process_time
from skimage import morphology
from skimage import measure
from scipy.ndimage import gaussian_filter
from scipy.ndimage import morphology as morphologysci
from scipy.ndimage import measurements
from scipy.ndimage import distance_transform_edt
from scipy import interpolate

# %%  
def LargestConnComp(ImgBi) :
    mask = np.zeros(np.shape(ImgBi))
    ConnectedComponent = measure.label(ImgBi)
    props = measure.regionprops(ConnectedComponent)
    area = [ele.area for ele in props]
    if area:
        largest_ind = np.argmax(area)
        CompLabel = props[largest_ind].label
        PixelList = ConnectedComponent==CompLabel
        mask = np.zeros(np.shape(ImgBi))
        mask[PixelList] = 1;
    return mask

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

# %% Breast_Segmentation function: takes about 12+-10 minutes per breast
def Breast_Segmentation(id, img,output_path):
    output_time = output_path + "/Segmentation_time.csv"
    
    t_start=process_time()
    now = datetime.datetime.now()
    print(id,'Segmentation starts at ',now.year, now.month, now.day, now.hour, now.minute, now.second)
    
    # Thresholds and image processing
    thSi = 130
    thBone = 100 # 70 is challenging
    thMu_seed = 15
    thMu = -100 # -80 is too challenging
    thSkin = -180
    thBreast = -500 # set by ABCT
    thAir = offset = -2**10
    
    HUGland = -41 # SD=10
    HUFat = -255 # SD=15
    
    # when img_path is pluged: img=np.load(img_path).astype(np.int16)
    L,R,C = np.shape(img)
    mask = np.zeros(np.shape(img),dtype=np.uint8)
    
# %% Preprocessing 

    # Breast and nipple masking
    breast = img>thBreast
    breast = np.bool_(breast)
    mask[breast]=1

# %% Fold segmentation  
    tFold_start=process_time()
    tFold = 0
    
    com_chest = measurements.center_of_mass(breast[0])
    com_chest = np.array(com_chest,dtype='int') #  x,y
    
    lenfold = int(0.2*L)
    nlist = []
    fold = np.zeros(np.shape(img),dtype=np.bool_)
    
    # Investigate if there is a non-breast fold section
    for z in range(lenfold): 
        ConnectedComponent = measure.label(breast[z])
        props = measure.regionprops(ConnectedComponent)
        area = [ele.area for ele in props]
        nparea = np.array(area)
        n = len(nparea[nparea>1000])
        nlist.append(n)
    
    # Segment the fold if exists
    if any(i > 1 for i in nlist):
        breast_wo_fold = np.zeros(breast.shape)
        compmarker = np.zeros(breast[:lenfold].shape)
        
        dist = distance_transform_edt(breast[:lenfold])
        
        s = np.argmax(nlist)
        ConnectedComponent = measure.label(breast[s])
        compmarker[s] = ConnectedComponent
        
        label = morphology.watershed(-dist, compmarker, mask=breast[:lenfold])
        templabel = label[0,int(com_chest[0]),int(com_chest[1])]
        breast_wo_fold[:lenfold] = label == templabel
        fold[:lenfold] = np.logical_and(breast[:lenfold],label!=templabel) # the rest is folds
        
        #mask[fold.astype(np.bool)] = 6
        
        tFold_stop=process_time()
        tFold = tFold_stop-tFold_start
        print('Fold section segmented in %.1f s.' %tFold)
    
# %% Skin segmentation: mask = 3
    tSkin_start=process_time()  
    
    for l in range(10):
        bL = L-l
        if np.sum(breast[bL-1]):
            break
    
    tip = measurements.center_of_mass(breast[bL-1])
    tip = np.array(tip,dtype='int') #  z,x,y
    
    # breasat direction vector
    sub = tip-com_chest
    k = [bL, sub[0], sub[1]]
    k_dot = k/np.sqrt(k[0]**2+k[1]**2+k[2]**2)
    
    # Nipple segmentation: 6 mm from the nipple
    z, x, y = np.indices(breast.shape)
    r = 30 # 0.3mm*20 = 9 mm 
    
    com_sphere = [L, tip[0], tip[1]] + np.array(k_dot*r,'int')
    mask_sphere = (z-com_sphere[0])**2+(x-com_sphere[1])**2+(y-com_sphere[2])**2 <= (r*2)**2
    
    # open CV data type modification: to get the countours
    cvmask = np.array(breast * 255, dtype = np.uint8) # 3D mask
    cvimg = np.array((img-offset)/2**3) # 3D image
    cvimg = np.clip(cvimg, 0, 255).astype(np.uint8)
    
    bg = np.zeros(cvimg[0].shape, dtype=np.uint8) # 2D bg
    cvthBreast=(thBreast+2**10)/2**3
    
    n=3 # 1 for the boundary / 3 or 4 for the defualt thinnest skin thickness
    b3Skin = np.zeros(cvmask.shape)
    for s in range(L): # slice by slice
        contours, hierarchy = cv.findContours(cvmask[s], cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
        cont=cv.drawContours(bg.copy(),contours,-1,(255,255,255),2*n-1)
        if n > 1:
            cont = np.logical_and(cont, cvmask[s])
        b3Skin[s] = cont
    
    n=9 # thickest skin around nipple
    bmaxSkin = np.zeros(cvmask.shape)
    for s in range(L):
        contours, hierarchy = cv.findContours(cvmask[s], cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
        cont=cv.drawContours(bg.copy(),contours,-1,(255,255,255),2*n-1)
        if n > 1:
            cont = np.logical_and(cont, cvmask[s])
        bmaxSkin[s] = cont
        
    n=25 # thickest skin closed to the wall
    bMAXSkin = np.zeros(cvmask.shape)
    for s in range(L):
        contours, hierarchy = cv.findContours(cvmask[s], cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
        cont=cv.drawContours(bg.copy(),contours,-1,(255,255,255),2*n-1)
        if n > 1:
            cont = np.logical_and(cont, cvmask[s])
        bMAXSkin[s] = cont
    
    # Skin segmentation using watershed
    
    bSkin = b3Skin # only for patient 2 in Final selection has b1Skin. the rest b2Skin
    temp = img>thSkin
    temp = np.logical_and(temp,mask>0)
    temp = np.logical_or(temp,bSkin) # should be 2 or 3 in order not to have the gaps in the skin...
    nipple = np.logical_and(temp,mask_sphere)
    
    tempskinl = morphologysci.binary_erosion(temp, structure=np.ones((5,5,5)))
    for s in range(2):
        tempskinl[s] = morphologysci.binary_erosion(temp[s], structure=np.ones((5,5)))
    tempskin = np.logical_or(tempskinl,b3Skin) # b3 for watershed compromisation
    tempskin = np.logical_and(tempskin,np.logical_not(mask_sphere))
    
    tempskinbg = np.logical_or(tempskin,np.logical_not(breast))
    dist = distance_transform_edt(tempskinbg)
    
    
    local_maxi=dist>15 # 10 skin-gland connected. 20 gland disconnected
    markers = measure.label(local_maxi)
    
    label = morphology.watershed(-dist, markers, mask=tempskinbg)
    CompLabel =  label[0,0,0]
    skin = label==CompLabel
    skin = np.logical_and(skin,breast)
    skin = np.logical_or(skin,bSkin)
    
    # shaping after watershed
    non_skin0 = np.zeros(np.shape(skin))
    non_skin = np.zeros(np.shape(skin))
    non_non_skin = np.zeros(np.shape(skin))
    skinShape = np.zeros(np.shape(skin))
    skindiv = np.int(0.4*L)
    
    for s in range(L): 
        skinShape[s]=morphologysci.binary_dilation(skin[s], structure=np.ones((7,7)))
        skinShape[s] = np.logical_and(skinShape[s],temp[s])
        skinShape[s] = np.logical_or(skinShape[s],bSkin[s])
        
        # filter once again to remove the ditached pixels
        label = measure.label(skinShape[s])
        
        t=np.sum(skinShape[s],axis=1)
        ss=np.argmax(t)
        sss=skinShape[s,ss,:]
        
        CompLabel =  label[ss,np.argmax(sss)]
        skinShape[s] = label==CompLabel
        
        #filter once again to include the skin islands: for image 2, this should be excluded for the better matching with b1skin. otherwise, b2skin should be used once agina!!!!!!
        non_skin0[s] = np.logical_and(breast[s],np.logical_not(skinShape[s]))
        non_skin[s] = LargestConnComp(non_skin0[s])
        
        non_non_skin[s] = np.logical_and(breast[s],np.logical_not(non_skin[s]))
        
        skinShape[s] = np.logical_or(skinShape[s],non_non_skin[s])
        
        if not tFold:
            skinShape[s] = np.logical_and(skinShape[s],bMAXSkin[s])
    
    skinShape[skindiv:] = np.logical_and(skinShape[skindiv:],bmaxSkin[skindiv:])
    skinShape = np.logical_or(skinShape,bSkin)
    
    skinShape = np.logical_or(skinShape,nipple)
    skinShape = np.logical_and(skinShape,mask==1)
    mask[skinShape]=3
    
    tSkin_stop=process_time()
    tSkin = tSkin_stop-tSkin_start
    print('Skin segementation done in %.1f s.' %tSkin)

# %% Rib/Silicon implant segmentation: mask = 5
    tRS_start=process_time()
    tSi_stop = 0
    second_Rib = 0
    tRib = 0
    tSi = 0
    
    LSi = []
    Lrib = []
    
    temp = img>thBone
    temp = np.logical_and(temp,mask==1)
    temp = np.logical_and(temp,np.logical_not(bmaxSkin))
    
    if np.sum(LargestConnComp(temp)[0])>50:
        Lrib = np.argwhere(np.sum(temp,axis=(1,2))>0)
        Lrib = len(Lrib)
        rib = np.zeros(np.shape(temp))
        
        temp1 = LargestConnComp(temp)
        
        # investigate silicon implant
        if np.sum(temp1) > 150**3: # (4.5cm)^3 # Silicone exists
            tempSi = img>thSi
            setemp=morphologysci.binary_erosion(tempSi, structure=np.ones((5,5,5)))
            for s in range(2):
                setemp[s]=morphologysci.binary_erosion(tempSi[s], structure=np.ones((5,5)))
            sletemp = LargestConnComp(setemp)
            sdletemp=morphologysci.binary_dilation(sletemp, structure=np.ones((5,5,5)))
            Si=np.logical_and(sdletemp,tempSi)
            
            LSi = np.argwhere(np.sum(Si,axis=(1,2))>0)
            LSi = len(LSi)
            for s in range(LSi):
                Si[s] = morphologysci.binary_fill_holes(Si[s],np.ones((3,3)))
            
            mask[Si] = 7
            tSi_stop=process_time()
            
            # investigate additional ribs besides the implant
            etemp=morphologysci.binary_erosion(temp1, structure=np.ones((5,5,5)))
            for s in range(2):
                etemp[s]=morphologysci.binary_erosion(temp1[s], structure=np.ones((5,5)))
            letemp = LargestConnComp(etemp)
            dletemp=morphologysci.binary_dilation(letemp, structure=np.ones((5,5,5)))
            temp1=np.logical_and(dletemp,temp1)
            
            for s in range(Lrib):
                temp1[s] = morphologysci.binary_fill_holes(temp1[s],np.ones((3,3)))
        
        else: # Only Ribs
            for s in range(Lrib):
                temp1[s] = morphologysci.binary_fill_holes(temp1[s],np.ones((3,3)))
            rib = temp1
        
        # shape the additioinal ribs
        temp2 = np.logical_and(temp,np.logical_not(temp1))
        rib2 = LargestConnComp(temp2)
        if np.sum(rib2) > 15**3: # (4.5mm)^3
            second_Rib = 1
            for s in range(Lrib):
                rib2[s] = morphologysci.binary_fill_holes(rib2[s],np.ones((3,3)))
            rib = np.logical_or(rib,rib2)
        
            temp3 = np.logical_and(temp2,np.logical_not(rib))
            temp3 = LargestConnComp(temp3)
            if np.sum(temp3) > 15**3: # (4.5mm)^3
                for s in range(Lrib):
                    temp3[s] = morphologysci.binary_fill_holes(temp3[s],np.ones((3,3)))
                rib = np.logical_or(rib,temp3)
        
        rib = np.bool_(rib)
    
        mask[rib] = 5
    
        tRS_stop=process_time()
        tRS = tRS_stop-tRS_start
        
        if tSi_stop*second_Rib: # if both rib and silicone exist: Silicone is calculated seperately
            tSi = tRS_start-tSi_stop
            tRib = tRS-tSi
            print('Rib and silicon segmented in %.1f s.'%tRS)
        elif tSi_stop: # Silicon only
            tSi = tRS
            print('Silicone segmented in %.1f s.'%tSi)
        else: # Rib only
            tRib = tRS
            print('Rib segmented in %.1f s.'%tRib)
# %% Muscle segmentation: mask = 4   
    tMuscle_start=process_time()
    tMuscle = 0
    
    # Pectoralis muscle seed investigation
    seedMu = np.ones(img[0].shape)
    
    if LSi: # when silicon implant presents
        comp = img[0]>thMu
        comp = LargestConnComp(comp)
        sub = np.logical_and(comp,np.logical_not(Si[0]))
        if np.sum(sub) > 80**2: # (2.4cm)^2 Tunning!!!!
            seedMu = Si[0]
            
    elif Lrib: # when rib presents
        seedMu = temp1[0]
        
    else: # when no hard materials present
        seedMu = np.zeros(np.shape(img),dtype=np.int16)
    
        imgG = gaussian_filter(img[0], sigma=5)
        temp = np.logical_and(imgG>-5,breast[0]) # -5 for patient 7,8,10....in final list, not 15
        if np.sum(temp):
            temp = LargestConnComp(temp)
            if not np.sum(temp*mask[0]==3) :  
                seedMu = temp
                
    # Muscle segmentation
    if np.sum(seedMu):
        com_seedMu=measurements.center_of_mass(seedMu)
        com_seedMu = np.array(com_seedMu,dtype='int')
        
        temp = img>thMu
        tempf = np.zeros(np.shape(temp))
        tempef = np.zeros(np.shape(temp))
        markers = np.zeros(np.shape(temp))
    
        for s in range(L): # slice by slice due to fibroglands
            tempf[s] = morphologysci.binary_fill_holes(temp[s],np.ones((3,3))).astype(int) 
    
        tempef = morphologysci.binary_erosion(temp, structure=np.ones((5,5,5)))
        for s in range(2):
            tempef[s] = morphologysci.binary_erosion(temp[s], structure=np.ones((5,5)))
    
        for s in range(L):
            tempef[s] = morphologysci.binary_fill_holes(tempef[s],np.ones((3,3))).astype(int)     
    
        dist = distance_transform_edt(tempf)
        
        dist_mark = dist[0]
        temp_dist = dist_mark[dist_mark>0] # flattening
        temp_dist = temp_dist.tolist()
        temp_dist.sort()
        th_dist = temp_dist[-len(temp_dist)*2//5] # highest 40% of the pixels
        
        # either use the markers from slice 0 or from higher dist threshold
        local_maxi = dist_mark>th_dist # 10 did not work for forked muscle # 8 does not work for muslce having holes
        markers[0] = measure.label(local_maxi)
        
        if LSi: # In order to reduce the muscle segmentation error due to the distortion by the silicon implant
            Si_extended =  morphologysci.binary_dilation(Si, structure=np.ones((19,19,19)))
            Si_extended[:L//5] = np.ones((L//5,R,C))
            dist = dist*Si_extended
            tempef = np.logical_and(tempef,Si_extended)
            
        label = morphology.watershed(-dist, markers, mask=tempef)
    
        CompLabel =  label[0,com_seedMu[0],com_seedMu[1]]
        muscle = label==CompLabel
        
        # remove the small virture muscle forks (glands)
        for s in range(np.sum(np.sum(np.sum(muscle,1),1)>0)):
            if np.sum(muscle[s])<100:
                # from the slice after s
                muscle[s:]=np.zeros((L-s,C,C))
                break
        
        muscled = morphologysci.binary_dilation(muscle, structure=np.ones((7,7,7)))
        for s in range(L):
            muscled[s] = morphologysci.binary_fill_holes(muscled[s],np.ones((3,3))).astype(int) 
            
        muscled = np.logical_and(muscled,tempf)
    
        label = measure.label(muscled)
        CompLabel =  label[0,com_seedMu[0],com_seedMu[1]]
        muscled = label==CompLabel
    
        if Lrib:
            muscled = np.logical_and(muscled,mask==1)
        
        mask[muscled] = 4
    
        Lmu = np.argwhere(np.sum(muscled,axis=(1,2))>0)
        Lmu = len(Lmu)
    
        tMuscle_stop=process_time()
        tMuscle = tMuscle_stop-tMuscle_start
        print('Muscle segmented in %.1f s.'%tMuscle)
        
    # %% Density Calculation
    HUGland = -41 # SD=10
    HUFat = -255 # SD=15
    
    HU_tissue = img[mask==1]
    
    img_tissue = img.copy()
    img_tissue[mask!=1] = offset
    
    meanHU = np.mean(HU_tissue)
    meanGland = HUcalibration(meanHU)
    print('Glandularity: %.3f'%meanGland)
    
    # %% Glands Segmentation: Regiaon Growing: mask = 2
    tGland_start = process_time()
    
    # Finding the gland seed
    gland_list_th = []
    th_list = []
    th = -51 # 100% glandularity HU value within SD
    temp = img_tissue >= th 
    
    # if there is detected glands
    if np.sum(temp) > 0:
        lab = []
        seg_n=np.zeros(img_tissue.shape)
        seg_p=np.zeros(img_tissue.shape)
    
        label = measure.label(temp)
        props = measure.regionprops(label)
        coord_seed = [ele.coords[0] for ele in props] # 1st voxle of each element
        
        
        # Puedo region growing to boost the calculation
        gland_th = np.sum(temp)/HU_tissue.shape[0]
        gland_list_th.append(gland_th)
        th_list.append(th)
        step = 5 # 2.5% of the glandularity
        
        for i in range(20):
            th -= step
            temp = img_tissue >= th   
            gland_th = np.sum(temp)/HU_tissue.shape[0]   
            th_list.append(th)
            gland_list_th.append(gland_th)
            
            if gland_th > meanGland:
                break
        
        # Two maps to applying region growing by comparing from the seed 
            # negetive
        temp = img_tissue >= th_list[-2]
        label = measure.label(temp)
        props = measure.regionprops(label)        
        
        for e in range(len(coord_seed)): # select the label elements belonging to the seed region (a few)
            lab.append(label[coord_seed[e][0],coord_seed[e][1],coord_seed[e][2]])
        
        for m in range(np.unique(lab).shape[0]): # OR operate for all masks of the label elements
            seg_n = np.logical_or(seg_n,label==np.unique(lab)[m])
        
        gland_reg_n = sum(sum(sum(seg_n)))/HU_tissue.shape[0]
        
        coord_n = [] # 1st voxle of each element
        for m in range(np.unique(lab).shape[0]): # OR operate for all masks of the label elements
            coord_n.append(props[np.unique(lab)[m]-1].coords[0]) # list
            
            # positive
        temp = img_tissue >= th_list[-1]
        label = measure.label(temp)
        props = measure.regionprops(label)        
        
        for e in range(len(coord_n)): # select the label elements belonging to the seed region (a few)
            lab.append(label[coord_n[e][0],coord_n[e][1],coord_n[e][2]])
        
        for m in range(np.unique(lab).shape[0]): # OR operate for all masks of the label elements
            seg_p = np.logical_or(seg_p,label==np.unique(lab)[m])
        
        gland_reg_p = sum(sum(sum(seg_p)))/HU_tissue.shape[0]
        
        if abs(gland_reg_n-meanGland)>abs(gland_reg_p-meanGland):
            seg = seg_p
        else: 
            seg = seg_n
        
        mask[seg]=2
    
    tGland_stop=process_time()
    tGland = tGland_stop-tGland_start
    print('Glands segmented in %.1f s.'%tGland)
    
    tTotal = tGland_stop-t_start
    print('Segmentation completed in total %.1f s.'%tTotal)
    
    #%% CSV file for time record.
    fnames = ['id', 'Total', 'Glands', 'Skin', 'Muscle', 'Rib', 'Fold', 'Silicone']
        # first time only: header
    if not os.path.isfile(output_time):
        f = open(output_time, 'w')
        with f:       
            writer = csv.DictWriter(f, fieldnames=fnames)    
            writer.writeheader()    
            
    with open(output_time, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fnames)
        writer.writerow({fnames[0]:id, fnames[1]:"%.1f" % tTotal, fnames[2]:"%.1f" %tGland, fnames[3]:"%.1f" %tSkin, fnames[4]:"%.1f" %tMuscle, fnames[5]:"%.1f" %tRib, fnames[6]:"%.1f" %tFold, fnames[7]:"%.1f" %tSi})
        
    print('Segmentation image and csv file saved.')
    now = datetime.datetime.now()
    print(id,'Ends at ',now.year, now.month, now.day, now.hour, now.minute, now.second)
    
    return mask,fold;
