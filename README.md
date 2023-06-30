# Breast_segmentation
Fully automated breast segmentation for a novel breast CT

The 5 original scripts for the breast CT segmentation containts 5 pythong scripts.
Each script has to run in the order below as following for the described fucntion.

1. DICOM_to_numpy.py
Retrieving the DICOM images data to numpy arrays
Input files

2. Function_Breast_Segmentation.py
The function used for Breast_Segmentation_run_batch.py

3. Breast_Segmentation_run_batch.py
Segment the numpy images to numpy segmented images

4. Image_extraction.py
3 plane images extracted from the numpy images

5. Numpy_to_DICOM.py
save the numpy image to DICOM file using the default DICOM header.
