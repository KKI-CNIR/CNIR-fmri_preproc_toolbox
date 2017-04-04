
%% Steps to perform and Options for each

%Computer capability
low_ram = 0;
prefix = '';

%PAR/REC to nifti convert
parrec2nii = 1;
dcm2nii_toolbox = fullfile(cwd, 'dcm2niix', 'bin'); %only used to convert par/rec

%% Set origin - Setting this will make the program wait for user input!
%%Set origin of ACPC line to 0
set_origin = 0;

%Orient (?)
correct_orient = 0; %Only needed for ANALYZE files (not needed for NIFTI, with proper orientation information)

%Slice time correct
slice_time_correct = 1;
%slice_order = 1; %1=ascending, 2=descending; 3=interleaved; 4=custom (file?)
slice_order=1;
    %1=ascending, 
    %2=descending; 
    %3=interleaved ascending starting at 1; 
    %4=interleaved ascending starting at 2;
    %5=interleaved descending starting at end; 
    %6=interleaved descending starting at end-1
    %7=custom file - use slice_order_file variable
slice_order_file = ''; % File for slice order - ignore for 
%ref_slice=-1; % -1 will set reference slice as (spatially) middle slice

%Motion correction
motion_correct = 1;
%threshold in mm to label spikes on 1-back motion plots
motion_thresh = 3;

%Save bias-corrected version of the anatomical image
unbiased = 0;

%Coregistration
skipanat = 1; %skip all steps relating to MPRAGE
coregister = 0; %1=Anatomical to Functional, 2=Functional to Anatomical

%Segment the MPRAGE
segment_anat = 0; 

%Spatial Normalization
normalize = 2; %1=Use unified segmentation (only write), 2=Use EPI template when anatomical is not used (estimate and write)
normalize_bb = [-90 -126 -72; 90 90 108]; %this is different from default bounding box
normalize_prefix = 'wepi'; %default prefix is 'w'
normalize_vox = [2 2 2]; %default

%Create brain mask file
create_brain_mask_file = 0; %This has to be done if hp_filter, nuisance_remove or bp_filter is set
brain_mask_file=fullfile(tooldir,'brain_mask.nii');

%High pass time domain filtering
hp_filter = 0;
%hp = 0.005; %Cutoff frequency 200s 
%5 = linear detrend
band = 5;

%Smoothing for ICA - Can be skipped without affecting further steps
ica_smooth = 0;
ica_smooth_fwhm = 6; %-1 will set fwhm to twice the current voxel size, so if performed after normalization, would be 4-mm
%fwhm = 2; %twice the voxel size - check how to specify this

%% Nuisance removal - typically performed only for seed based analysis
nuisance_estimate = 1;
nuisance_regress = 1;
nuisance_file_postfix='nuisance_fdm_noGSR_pwm50_pcsf50';
nr_options.motion_params=1; %1;
nr_options.detrend_motion_params=1;
nr_options.filter_motion_params=0;
nr_options.filter_motion_params_cutoff=0.005;
nr_options.filter_motion_params_band='high';
nr_options.diff_motion_params=1; %Use motion params differential
nr_options.square_motion_params=0; %Use motion params square
nr_options.global_signal=0;
nr_options.wm_signal=2; %1=mean, 2=pca, 0=not included
nr_options.wm_pca_percent=.50;
nr_options.csf_signal=2; %1=mean, 2=pca, 0=not included
nr_options.csf_pca_percent=.50;
nr_options.ALVIN_HOME = fullfile(cwd, 'ALVIN_v1');

%Smoothing
rs_smooth = 1; %For FC analysis
smooth_fwhm = 6; % -1 will set fwhm to twice the voxel size
%fwhm = 2; %twice the voxel size

%Band-pass time domain filtering
bp_filter = 1; %For FC analysis
bp = [0.01 0.1]; %Filter cut-offs in Hz
%filter_order=2; %Filter order


