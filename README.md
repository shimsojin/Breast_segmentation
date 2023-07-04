# Breast_segmentation
Scripts written by Sojin Shim, Ph.D., Institute of Diagnostic and Interventional Radiology, University Hospital Zurich, 2018-2023

The 5 original scripts for the breast CT segmentation containts 5 python scripts.

* DICOM_to_numpy.py
Retrieving the DICOM images data to numpy arrays
Input files


* Breast_Segmentation_run_batch.py
Main script. Segment input numpy images and generate output files/
Input_path: directory including only breast CT numpy images, grey level in Hounsfield unit (HU), int16=[-2**10, 2*10] without the cylinder
	Consider pid as the unique number of the numpy file. 
Output_path: same as Input_path as default
Output files:
	segmentation_[pid]: segmented numpy images
	segmentation_fold_[pid]: additional segmented numpy images for skin fold depicted from abdominal or thoracic wall, in case exists 
	Segmentation_time.csv: csv file recording the segmentation time of each numpy file in the input directory
	Output segmentation code - 0: air / 1: fat / 2: glands / 3: skin / 4: pectoralis major muscle / 5: bone / 6: skin fold section / 7: implant / 9: cylinder

* Function_Breast_Segmentation.py
The function used for Breast_Segmentation_run_batch.py

* Image_extraction.py
Generate 3 representative images of the segmentation in coronal (c), sagittal (s),and transversal (t) planes.
Analyse the percentage breast density (glandularity, G) calculated based on HU values and segmentation ratios.
Input_path: directory including the segmented numpy images (Output_path of Breast_Segmentation_run_batch.py)
Output_path: same as Input_path as default
Output files:
	img_[pid]_(c/s/t).png: 3 plane representative breast images taken from the largest transection
	auto_seg_[pid]_(c/s/t).png: 3 plane representative segmentation images taken from the largest transection
	Glandularity_compariton.csv: csv file recording the glandularity (G) analysis for each numpy image in the input directory
		id=pid, HU_G=glandularity calculated based on HU, Seg_G=glandularity calculated based on segmentation

5. Numpy_to_DICOM.py
save the numpy image to DICOM file using the default DICOM header.

npimgs: list of nupy files' names with '.npy'
pnpimgs: numpy files' names with '.npy'
pnpimg[:-4]: numpy files' name 

