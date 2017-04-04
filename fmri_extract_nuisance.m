function fmri_extract_nuisance(func_files, tr, rp_file, wm_file, csf_file, brain_mask_file, options, nuisance_file)
%Function to extract nuisance time courses from the data and use PCA to
%reduce dimensions. These components are combined with the motion
%parameters and, optionally, global signal in to a combined nuisance matrix
%
%Usage
%   fmri_extract_nuisance(func_files,tr,rp_file,wm_file,csf_file,brain_mask_file, ...
%                         options,nuisance_file)
%Input variables
%   func_files - cell array of functional files (from one run)
%   tr         - Repetition time (time taken for one dynamic/volume acquisition
%   rp_file    - realignment parameters, i.e., file with 6 columns of realignment
%                parameters (such as what is obtained from spm_realign or mcflirt)
%   wm_file    - White matter segmentation result in same space as functional
%                (usually got from spm_segment or FSL's fast)
%   csf_file   - CSF segmentation result in same space as functional
%   brain_mask_file - filename of the brain_mask image (usually got from
%                fmri_create_brain_mask)
%   options    - a structure with fields listed below which specify what
%                nuisances to use (skip those that are not needed - atleast
%                one option is necessary!)
%                options.motion_params=1; 0 or 1 - to specify if motion parameters should be considered nuisance
%                options.filter_motion_params=0; 0 or 1 - to specify if
%                   motion_paraters need to be filtered - this is needed
%                   only when the data have been filtered along time in between
%                   realignment and nuisance extraction/removal steps
%                options.filter_motion_params_cutoff=0.005; Cut-off
%                   frequency in Hertz for filtering (use the same cut-off
%                   used for the data)
%                options.filter_motion_params_band='high'; 'high' or
%                   'low' or 'band' or 'stop' - use the same as was used for the data
%                options.detrend_motion_params=0; 0 or 1 - to speficy if
%                   motion parameters need to be detrended - this is needed
%                   only when the data have been detrended along time between
%                   realignment and nuisance extraction/removal
%                options.diff_motion_params=1; %Use one-back difference motion params?
%                options.square_motion_params=0; %Use quadratic motion params?
%                options.global_signal=0; %Use global signal?
%                options.wm_signal=2; %1=mean of WM regions, 2=PCA of WM voxel timecourses
%                options.wm_pca_percent=0.50; % percent signal to use when
%                   using pca on WM; if not specified, will default to top 5
%                   principal components
%                options.csf_signal=2; %1=mean of CSF region, 2=PCA of CSF voxel timecourses
%                options.csf_mask=1; %1=bwperim erode, 2=imerode, 3=just threshold
%                options.csf_pca_percent=0.50; % percent signal to use when
%                   using pca; if not specified, will default to top 5
%                   principal components
%                options.ALVIN_HOME - path to the ALVIN toolbox
%   nuisance_file - the output filename where the nuisance mat file will be written
%MB Nebel Sep 19, 2014 - No major changes; removed obsolete information
%from help section
%MB Nebel Nov 20, 2013 - Saves percent variance explained by white matter
%and csf regressors to nuisance .mat
%MB Nebel July 16, 2013 - Added option to include first difference of mean
%signals from white matter and csf
%MB Nebel March 5, 2013 - if wm_mask or csf_mask contain fewer than 5
%voxels, force to ignore during nuisance regression
%MB Nebel Jan 8, 2013 - if wm_pca_percent and csf_pca_percent are not
%defined in options, select first 5 principle components of each as
%recommended by Behzadi et al 2007 & Whitfield-Gabrieli et al 2012.
%MB Nebel Dec 6, 2012 - fixed error in the way covariance matrices and principal component variances are calculated; 
%   limited csf_mask to voxels within ALVIN ventricle mask regardless of mask type chosen
%MB Nebel Sept 20, 2012 - Added a new option - csf_mask to indicate type of erosion to use and changed output type of wm and csf masks to .nii
%Suresh E Joel Aug 22, 2011 - Added a few more options - detrend mp,
%   ability to not specify fields in option, if it is not needed
%Suresh E Joel Jul 28, 2011 - Cleaned code and commented v
%Suresh E Joel May 20, 2011 - Toolbox (fmri_preproc_toolbox) compatibility
%Suresh E Joel, Mar 7 2011 - added options (select which nuisances to extract)
%Suresh E Joel, Mar 2, 2011 - cleaned up code, no major changes. Removed creation of brain mask (reads already created file).
%Suresh E Joel Dec 21, 2010 - check if CSF (or WM) mask is empty - ignore
%   empty masks and write nuisance names
%Suresh E Joel Dec 21, 2010 - Changed to 4D and separate CSF and WM PCA;
%   changed CSF threshold to 99% instead of 99.5%
%MB Nebel, Aug 23, 2010 - Added line to create a brain mask; changed the way the script looks for tissue segmentation files
%Suresh E Joel, Aug 19, 2010 - Changed for using in RS_LDDMM study
%Suresh E Joel, July 10, 2009, Back up for Change orient, COMPCOR MASKS (Erosion and Thresholding change.)
%Priti, Last Modified 05/26/09, Principal Components explaining 95% variance
%Priti, Last edited on 03/05/09, Added PCA and Detrending and Combined all the Nuisance Parameters
%Suresh E Joel, Feb 11, 2008; Priti Srinivasan Jun, 2009

disp(['Extracting Nuisances: ',fileparts(func_files{1})]);

P=strvcat(func_files);%#ok
so = get_orient(P(1,:));
nui.names=[];
nui.tc=[];

if(~isfield(options, 'ALVIN_HOME'))
    options.ALVIN_HOME = fullfile(fileparts(pwd), 'ALVIN_v1');
end

alvin_mask = spm_read_vols(spm_vol(fullfile(options.ALVIN_HOME, 'ALVIN_mask_v1.img')));

%% Motion parameters
if(isfield(options,'motion_params') && options.motion_params==1)
    nm=load(rp_file);
    if(isfield(options,'detrend_motion_params') && options.detrend_motion_params==1)
        nm=detrend(nm);
    end;
    if(isfield(options,'filter_motion_params') && options.filter_motion_params==1)
        [b,a]=butter(2,2*tr*options.filter_motion_params_cutoff, ...
            options.filter_motion_params_band); %High pass filter
        nm=filtfilt(b,a,nm);
    end
    nm=detrend(nm,'constant'); %demeaning - done in all cases
    nui.names={'mp_x','mp_y','mp_z','mp_yaw','mp_pitch','mp_roll'};
    nui.tc=nm;
    if(isfield(options,'diff_motion_params') && options.diff_motion_params==1)
        nui.tc=[nui.tc, [0,0,0,0,0,0; diff(nui.tc)]];
        n=1;for i=(length(nui.names)+1):(length(nui.names)+6)
            nui.names{i}=['d',nui.names{n}];n=n+1; end;
    end;
    if(isfield(options,'square_motion_params') && options.square_motion_params==1)
        nui.tc=[nui.tc,nm.^2];
        n=1;for i=(length(nui.names)+1):(length(nui.names)+6)
            nui.names{i}=['s',nui.names{n}];n=n+1; end;
    end;
end; %End motion params


%% Read White matter segmentation result and make eroded mask
if(isfield(options,'wm_signal') && options.wm_signal~=0)%Make mask for either mean or PCA
    if(~strcmp(get_orient(wm_file), so))% If reorienting is necessary
        warning('Mywarn:CheckOrientation','Orientations of the Functional and WM mask are not the same');
        change_orient(wm_file,so);
    end;
    Vm = spm_vol(wm_file);
    wm_mask = spm_read_vols(Vm);
    wm_mask = (wm_mask/max(wm_mask(:)))>0.99;
    wm_mask = double(imerode(logical(wm_mask), ones(3,3,3)));
    if(nnz(wm_mask)<5), %IF WM mask has fewer than 5 voxels, force ignore options.wm_signal
        options.wm_signal=0;
        warning('MyWarn:MaskEmpty','WM mask is empty');
    end;
    % Write out the eroded WM mask
    Vm.fname = fullfile(fileparts(Vm.fname),'wm_erode_mask.nii');
    Vm.private.dat.fname = Vm.fname;
    spm_write_vol(Vm,wm_mask);
    clear Pm Vm
end;

%% Read CSF segmentation result and make eroded mask
if(isfield(options,'csf_signal') && options.csf_signal~=0) %Make mask for either mean or PCA
    if(~strcmp(get_orient(csf_file),so)) % If reorienting is necessary,
        warning('Mywarn:CheckOrientation','Orientations of the Functional and CSF mask are not the same');
        change_orient(csf_file,so);
    end;
    Vm = spm_vol(csf_file);
    csf_mask = spm_read_vols(Vm);
    csf_mask = (csf_mask/max(csf_mask(:)))>0.99;% Threshold the CSF mask at 99%
    csf_mask = (csf_mask .* alvin_mask) > 0; %multiply by alvin_mask to exclude everything other than ventricles
    if(isfield(options, 'csf_mask') && options.csf_mask==1)
        csf_mask = csf_mask-bwperim(csf_mask,6); %Erode
    elseif(isfield(options, 'csf_mask') && options.csf_mask==2)
        csf_mask = double(imerode(logical(csf_mask), ones(3,3,3)));
    end
    if(nnz(csf_mask) < 5) %If CSF mask is empty force ignore options.csf_signal
        options.csf_signal=0;
        warning('MyWarn:MaskEmpty','CSF mask contains < 5 voxels');
    end;
    % Write out eroded CSF mask
    Vm.fname = fullfile(fileparts(Vm.fname),'csf_erode_mask.nii');
    Vm.private.dat.fname = Vm.fname;
    spm_write_vol(Vm,csf_mask);
    clear Vm;
end;


%% If using wm, csf or global signal
if((isfield(options,'wm_signal') && options.wm_signal ~=0) || ...
        (isfield(options','csf_signal') && options.csf_signal ~= 0) || ...
        (isfield(options,'global_sigal') && options.global_signal ~= 0))
    
    % Read the functional image headers
    V = spm_vol(P); 
    
    % Read brain_mask image
    Vm=spm_vol(brain_mask_file);
    brain_mask=logical(spm_read_vols(Vm));
    clear Vm;
    
    %% Extract Nuisance Time courses from Data
    Y=spm_read_vols(V);
    if(options.wm_signal),
        tc_wm=zeros(size(Y,4),sum(wm_mask(:)));
        wm_mask=logical(wm_mask);
    end;
    if(options.csf_signal),
        tc_csf=zeros(size(Y,4),sum(csf_mask(:)));
        csf_mask=logical(csf_mask);
    end;
    if(options.global_signal),
        mY=zeros(size(Y,4),1);
    end;
    for i=1:size(Y,4),
        tY=Y(:,:,:,i);
        tY(isnan(tY))=0;
        if(options.global_signal~=0), mY(i) = mean(tY(brain_mask))'; end;
        if(options.wm_signal), tc_wm(i,:) = tY(wm_mask); end;
        if(options.csf_signal), tc_csf(i,:) = tY(csf_mask); end;
    end;
    if(options.global_signal), mY=detrend(mY,'constant'); end; %demean
    clear tY;
    
   
    %% Compute Mean/PCA of the WM and CSF separately
    % WM signals
    if(options.wm_signal > 0)
        %center tc_wm by subtracting off column means
        tc_wm = bsxfun(@minus, tc_wm, mean(tc_wm, 1));
        if(options.wm_signal==1) %Mean
            npca_wm=mean(tc_wm,2);
            if (isfield(options,'wm_diff') && options.wm_diff ~=0)
                wm_diff = [0; diff(npca_wm, 1, 1)];
                npca_wm = [npca_wm, wm_diff];
            end
            
        else %PCA
            [u,s] = svd(cov(tc_wm'));
            eigen_values = diag(s);
            if(isfield(options,'wm_pca_percent'))
                ncomp=find(cumsum(eigen_values)/sum(eigen_values) > options.wm_pca_percent, 1, 'first' );
            else %take 1st five as recommended by Behzadi & Whitfield-Gabrielli
                ncomp = 5;  
            end
            npca_wm = u(:,1:ncomp);
            nui.wm_percent = sum(eigen_values(1:ncomp))/sum(eigen_values);
                
            clear u s ncomp;
        end;
        nui.tc=[nui.tc,npca_wm];
        for i=1:size(npca_wm,2)
            nui.names{1+length(nui.names)}=['wm',num2str(i)];
        end;
    end;
    
    
    % CSF signals
    if(options.csf_signal > 0)
        %center tc_csf by subtracting off column means
        tc_csf = bsxfun(@minus, tc_csf, mean(tc_csf, 1));
        if(options.csf_signal==1) %Mean
            npca_csf=mean(tc_csf,2);
            if (isfield(options,'csf_diff') && options.csf_diff ~=0)
                csf_diff = [0; diff(npca_csf, 1, 1)];
                npca_csf = [npca_csf, csf_diff];
            end
            
        else % PCA
            [u,s] = svd(cov(tc_csf'));
            eigen_values = diag(s);
            if(isfield(options,'csf_pca_percent'))
                ncomp=find(cumsum(eigen_values)/sum(eigen_values) > options.csf_pca_percent, 1, 'first' );
            else
                ncomp = 5;
                
            end
            npca_csf = u(:,1:ncomp);
            nui.csf_percent = sum(eigen_values(1:ncomp))/sum(eigen_values);
            clear u s ncomp
        end;
        nui.tc=[nui.tc,npca_csf];
        for i=1:size(npca_csf,2)
            nui.names{1+length(nui.names)}=['csf',num2str(i)];
        end;
    end;
    
    % Global signal
    if(isfield(options,'global_signal') && options.global_signal==1)
        nui.tc=[nui.tc,mY];
        nui.names{1+length(nui.names)}='global';
        if (isfield(options,'global_diff') && options.global_diff ~=0)
            mY_diff = [0; diff(mY, 1, 1)];
            nui.tc = [nui.tc, mY_diff];
            nui.names{1+length(nui.names)}='global_diff';
        end
    end;
    
end; %End wm || csf || gs ~= 0

%% Intensity Normalize the nuisance regressors (otherwise matrix inversion fails in Matlab
for i=1:size(nui.tc,2),
    nui.tc(:,i)=nui.tc(:,i)./max(squeeze(nui.tc(:,i)));
end;

%% Save the file
save(nuisance_file,'nui');
