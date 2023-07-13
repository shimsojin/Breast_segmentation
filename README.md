# Breast_segmentation
The breast segmentation scripts are written by Sojin Shim, Ph.D., Institute of Diagnostic and Interventional Radiology, University Hospital Zurich, 2018-2023 <br />
The scripts segment breast CT images in numpy format or convert the format to be runable by the segmentation script. The details of each script are following.

## DICOM_to_numpy.py
This script converts the DICOM to numpy with minor preprocessing - Removing the CT table cylinder and a few most posterior slices with suboptimal image quality for image analysis.<br />
* If the entire or part of preprocessing is not desired, the code has to be modified accordingly: remove the preprocessing part and include it in the segmentation function (Function_Breast_Segmentation.py)<br />
INPUTS::<br />
path_uids: path to the directory(ies) including DICOM image files with unique uids. (input path)<br />
output_path: same as the input path (path_uids) as the default.<br />
OUTPUTS::<br />
Converted numpy images with the corresponding uids.<br />
cylinder.npy: saves one plane of the cylinder structure for the image conversion to DICOM by Numpy_to_DICOM.py<br />

## Breast_Segmentation_run_batch.py
This is the main script. It applies the segmentation on the input numpy images and generates output files - output segmentation image and csv file indicating the time duration of the segmentation.<br />
INPUTS::<br />
Input_path: path to the directory including the breast CT numpy images for segmentation.<br />
	The numpy images are in grey level in Hounsfield unit (HU), int16=[-2**10, 2*10] without the cylinder.<br />
	Consider pid as the unique number of the numpy file. <br />
Output_path: same as Input_path as default.<br />
OUTPUTS::<br />
	segmentation_[uid]: segmented numpy images.<br />
	segmentation_fold_[uid]: additional numpy images segmenting the skin fold depicted from abdominal or thoracic wall, in case they exist<br />
	Segmentation_time.csv: csv file recording the segmentation time of each numpy file in second in the input directory<br />
	Output segmentation color code - 0: air / 1: fat / 2: glands / 3: skin / 4: pectoralis major muscle / 5: bone / 6: skin fold section / 7: implant / 9: cylinder<br />

### Function_Breast_Segmentation.py
The segmentation function applied in Breast_Segmentation_run_batch.py

## Image_extraction.py
This script generates 3 representative images of the segmentation in coronal (c), sagittal (s),and transversal (t) planes and analyze the percentage breast densities (glandularity, G) calculated based on HU values and segmentation ratios, respectively.<br />
INPUTS::<br />
Input_path: directory including the segmented numpy images (Output_path of Breast_Segmentation_run_batch.py).<br />
Output_path: same as Input_path as default.<br />
OUTPUTS::<br />
img_[uid]_(c/s/t).png: 3 plane representative breast images taken from the largest transection.<br />
auto_seg_[uid]_(c/s/t).png: 3 plane representative segmentation images taken from the largest transection.<br />
Glandularity.csv: csv file recording the glandularity (G) analysis for each numpy image in the input directory.<br />
id=uid, HU_G=glandularity [-] calculated based on HU, Seg_G=glandularity [-] calculated based on segmentation.

## Numpy_to_DICOM.py
This script save the numpy images to DICOM files using the default DICOM header, which can be replaced/customised. The converted DICOM files of segmented images can be used for Monte Carlo simulation offered by AB-CT - Advanced Breast CT GmbH. For the simulation purpose, the image of CT table cylinder is inserted back during the conversion. If the cylinder structure is not desired, it can be turned off.<br />
INPUTS::<br />
Input_path: directory including the segmented numpy images.<br />
Output_path: same as Input_path as default.<br />
dicom_file: dicom file including the meta data that is desired to be copied to the ouput dicom file<br />
cylinder: numpy cylinder image (attached in the package)<br />
Cylinder: switch the appearance of the cylinder structure in the final output image: 1 - On (default) / 0 - off<br />
OUTPUTS::<br />
uid directories: containing the corresponding DICOM files of the numpy image segmentation_[uid].
