# CNIR-fmri_preproc_toolbox

## preprocess_single_subject.m - Main function 1. 
Makes a copy of raw data, initializes preferences for processing single subject with potentially multiple runs of fMRI data, and saves preferences as .mat

**Usage: preprocess_single_subject(ID, func_files, sess_names, procdir, tooldir, rename_task, specs_file)**, where
	
	- ID: subject identification number padded with zeros to ensure equal length in format sub-ID
	- func_files: cell array of path/filename(s) for each raw run to be processed from given subject in the order they were collected
	- sess_names: cell array of same size as func_files indicating session names for each run
	- procdir: directory where raw data will be copied and processed; function creates ID folder within procdir if it does not already exist
	- tooldir: directory where processing scripts are located; default assumes fmri_preproc_toolbox exists within present working directory
	- rename_task: e.g., task-rest_bold or task-GNG_bold; script will copy & rename func_files to procdir/ID using the following convention             sub-ID/sess_name/run-?/sub-ID_sess_name_rename_task_run-?, where run number is determined by the order of func_files
	- specs_file: path/file.m containing processing preferences; default: pwd/fmri_preprocess_par2ica_spm12.m
	
## fmri_preprocess_spm12.m - Main function 2. 
Called by preprocess_single_subject.m Performs processing steps specified in specs_file. Requires [SPM 12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) be added to the Matlab path. If starting from par/rec files, requires [dcm2niix](https://github.com/rordenlab/dcm2niix). If performing nuisance regression, requires [ALVIN lateral ventricle segmentation](https://sites.google.com/site/mrilateralventricle/). If raw data files do not contain correct TR info, must specify tr in specs_file.

**Usage: fmri_preprocess_spm12(setup_file)**, where 
	
	- setup_file: .mat created by preprocess_single_subject.m
	
Uses the following files (which you should not have to modify): 

	- template_slice_time_job_spm12.mat
	- template_realign_job_spm12.mat
	- fmri_plot_diff.m - plots 1-back and 2-back differences of 6 rigid body realignment parameters
	- QC_spm_check_reg.m - prints images of 1st & middle volume to check subject registration to MNI EPI & T1 template
	- template_smooth_job_spm12.mat
	
## fmri_preprocess_specs_par2ica_spm12.m - Example specs_file. 
Sets flags required by fmri_preprocess_spm12 to perform the following processing steps:

	- parrec2nii: Convert par/rec to nifti using dcm2niix
	- slice_time_correct: Also sets slice_order to ascending
	- motion_correct: Also sets threshold for between volume movements to 3-mm; movements above threshold will be indicated in QC.ps
	- normalize: Uses SPM's EPI template. Also sets bounding box, voxel size, and file prefix for normalized images
	- brain_mask_file: must be defined if performing temporal filtering or nuisance regression
	- hp_filter: voxelwise linear detrend data within brain_mask_file
	- ica_smooth: Use 6-mm FWHM kernel 
	
## batch_preprocess_list.m
Probably not very helpful for people outside of CNIR. Finds raw data for each subject specified in slist.txt and defines inputs required for preprocess_single_subject.m
	


