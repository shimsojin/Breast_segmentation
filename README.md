# Breast_segmentation
The breast segmentation scripts are written by Sojin Shim, Ph.D., Institute of Diagnostic and Interventional Radiology, University Hospital Zurich, 2018-2023

The scripts segment breast CT images in numpy format or convert the format to be runable by the segmentation script. The details of each script are following.

## DICOM_to_numpy.py
This script converts the DICOM to numpy with minor preprocessing - Removing the CT table cylinder and a few most posterior slices with suboptimal image quality for image analysis.
* If the entire or part of preprocessing is not desired, the code has to be modified accordingly: remove the preprocessing part and include it in the segmentation function (Function_Breast_Segmentation.py)
INPUTS::
path_uids: path to the directory(ies) including DICOM image files with unique uids. (input path)
output_path: same as the input path (puth_uids) as the default.
OUTPUTS::
Converted numpy images with the corresponding uids.

## Breast_Segmentation_run_batch.py
This is the main script. It applies the segmentation on the input numpy images and generates output files - output segmentation image and csv file indicating the time duration of the segmentation.
INPUTS::
Input_path: path to the directory including the breast CT numpy images for segmentation.
	The numpy images are in grey level in Hounsfield unit (HU), int16=[-2**10, 2*10] without the cylinder.
	Consider pid as the unique number of the numpy file. 
Output_path: same as Input_path as default.
OUTPUTS::
	segmentation_[uid]: segmented numpy images.
	segmentation_fold_[uid]: additional numpy images segmenting the skin fold depicted from abdominal or thoracic wall, in case they exist
	Segmentation_time.csv: csv file recording the segmentation time of each numpy file in second in the input directory
	Output segmentation color code - 0: air / 1: fat / 2: glands / 3: skin / 4: pectoralis major muscle / 5: bone / 6: skin fold section / 7: implant / 9: cylinder

### Function_Breast_Segmentation.py
The segmentation function applied in Breast_Segmentation_run_batch.py

## Image_extraction.py
This script generates 3 representative images of the segmentation in coronal (c), sagittal (s),and transversal (t) planes and analyze the percentage breast densities (glandularity, G) calculated based on HU values and segmentation ratios, respectively.
INPUTS::
Input_path: directory including the segmented numpy images (Output_path of Breast_Segmentation_run_batch.py).
Output_path: same as Input_path as default.
OUTPUTS::
img_[uid]_(c/s/t).png: 3 plane representative breast images taken from the largest transection.
auto_seg_[uid]_(c/s/t).png: 3 plane representative segmentation images taken from the largest transection.
Glandularity.csv: csv file recording the glandularity (G) analysis for each numpy image in the input directory.
id=uid, HU_G=glandularity [-] calculated based on HU, Seg_G=glandularity [-] calculated based on segmentation.

## Numpy_to_DICOM.py
This script save the numpy images to DICOM files using the default DICOM header, which can be replaced/customised. The converted DICOM files of segmented images can be used for Monte Carlo simulation offered by AB-CT - Advanced Breast CT GmbH. For the simulation purpose, the image of CT table cylinder is inserted back during the conversion. If the cylinder structure is not desired, it can be turned off.
INPUTS::
Input_path: directory including the segmented numpy images.
Output_path: same as Input_path as default.
dicom_file: dicom file including the meta data that is desired to be copied to the ouput dicom file
cylinder: numpy cylinder image (attached in the package)
Cylinder: switch the appearance of the cylinder structure in the final output image: 1 - On (default) / 0 - off
OUTPUTS:: 
