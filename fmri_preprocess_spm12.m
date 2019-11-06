function fmri_preprocess_spm12(setup_file)

%Function to preprocess a single subject with multiple fmri runs.
%The steps and specifications for what is to be done is from the set-up
%file (created from preprocess_single_subject.m)

%Suresh E Joel Dec 22, 2010, Modified Mar 2, 2011
%MB Nebel modified Jan 17, 2012: added flag to compile motion parameters
%MB Nebel modified Feb 19, 2013: fmri_regress_nuisance returns a variable
%that indicates the file prefix for nuisance-regressed data
%MB Nebel modified March 20, 2017: use dcm2niix to convert par/rec to nifti
%to avoid NIFTI_PAIR format; no longer cd to Template directory

%% Initializing
if(isunix)
    load(setup_file);
else
    load(char(setup_file));
end

%disp(['Resting State Preprocessing on: ',fileparts(func_files{1})]);

if ~exist('QC_File', 'var')
    old_wd = pwd;
    cd(fileparts(fileparts(func_files{1}{1})));
    QC_File = fullfile(pwd, 'QC.ps');
    cd(old_wd);
end

%% PAR/REC to NIFTI using r2agui
if(parrec2nii)
    
    %%have not updated anatomical conversion to use dcm2nii because we
    %%stopped using MPRAGE for spatial normalization
    if(~skipanat)
        % supply Raw folder
        % rec_convert needs an output folder for Anatomical
        [apath,fname,e]=fileparts(anat_file); %#ok<NODEF>
        if(~strcmp(e,'.rec'))
            error('Please input rec file to convert PAR/REC');
        end;
        anat_file=rec_convert(anat_file,3);
        anat_file=anat_file{1};
        [~,fname2,ext]=fileparts(anat_file);
        anat_file = fullfile(apath, [strrep(fname2, '_001', '') ext]);
        clear fname2;
        
        if(isunix)
            system(sprintf('mv "%s" "%s"', ...
                fullfile(apath, [fname, '_001.nii']), ...
                fullfile(apath, [fname, '.nii'])));
            system(sprintf('mv "%s" "%s"', ...
                fullfile(apath, [fname, '_001.img']), ...
                fullfile(apath, [fname, '.img'])));
            system(sprintf('mv "%s" "%s"', ...
                fullfile(apath, [fname, '_001.hdr']), ...
                fullfile(apath, [fname, '.hdr'])));
            % 	system('rm _001.*');
        else
            if(~strcmp(ext, '.nii'))
                [success, message, ~] = movefile(fullfile(apath, [fname, '_001.hdr']), fullfile(apath, [fname '.hdr']));
            end
            [success, message, ~] = movefile(fullfile(apath, [fname, '_001' ext]), fullfile(apath, [fname ext]));
        end %isunix
    end %~skipanat
    
    %use dcm2niix to convert functional par/rec to .nii
    for ir=1:length(func_files)
        [~,~,e]=fileparts(func_files{ir}{1});
        if(strcmp(e,'.rec') || strcmp(e, '.dcm'))
            %use dcm2nii to convert par/rec to .nii
            convert_str = sprintf('%s %s', fullfile(dcm2nii_toolbox, 'dcm2niix'), func_files{ir}{1});
            [status, result] = system(convert_str);
            
            %haven't been able to figure out how to specify the
            %output.nii that I want so move it after creation
            niilist = dir(fullfile(func_dir{ir},'*.nii'));
            %should only be one .nii if starting from raw data
            [success, message, ~] = movefile(fullfile(func_dir{ir}, niilist.name), fullfile(func_dir{ir}, nii_names{ir}));
            func_files{ir}{1} = fullfile(func_dir{ir}, nii_names{ir});
            
            %if don't have dcm2nii, can try rec_convert but then must save
            %data in NIFTI_PAIR format
            %             func_files{ir}=rec_convert(func_files{ir}{1},3);
        end;
    end;
end %if parrec2nii


%% Set origin
if(set_origin)
    %Anat file
    spm_image('init',anat_file);
    uiwait;
    %Func file
    spm_image('init',func_files{1}{1});
    uiwait;
    
    %     %Reorient the rest of the func files (only from run1)
    %     reorient_files=dir(fullfile(fileparts(func_files{1}{1}),'*_reorient.mat'));
    %     for i=1:length(reorient_files),
    %         filedate(i,:)=datestr(reorient_files(i).date,'yyyymmddHHMMSS');
    %     end;
    %     [jnk,ii]=sortrows(filedate);
    %     reorient_file=fullfile(fileparts(func_files{1}{1}),reorient_files(ii(end)).name);
    %     p=load(reorient_file,'-ascii');
    %     P=strvcat(func_files{1});
    %     P(1,:)=[];
    %     reorient_nifti(P,p);
end;

%% Open a sample image and try to get tr (assume same for all runs)
vs=spm_vol(func_files{1}{1});
%if you did not define tr in specs file, try to figure it out
if(~exist('tr', 'var'))
    tr=vs(1).private.timing.tspace;
    if(tr>10) %units are probably in msec
        tr=tr/1000;
        warning(['Verify TR is ',num2str(tr),' seconds']);
    elseif(tr<0.5)
        tr=tr.*1000;
        warning(['Verify TR is ',num2str(tr),' seconds']); %#ok<*WNTAG>
    end;
end;

% double check TR
if(tr<0.5)
    error(['WRONG TR: TR is ',num2str(tr),' seconds'])
end

%% Slice time correction
if(slice_time_correct)
    load(fullfile(tooldir, 'template_slice_time_job_spm12.mat')); % Load template
    for irun=1:length(func_files)
        vs=spm_vol(func_files{irun}{1});
        for iframe=1:length(vs)
            frames{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe)); %#ok<*AGROW>
        end;
        matlabbatch{1}.spm.temporal.st.scans{irun} = frames';
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{1}.spm.temporal.st.nslices=vs(1).dim(3); %Assumes that the 3rd dimension is slices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{1}.spm.temporal.st.tr=tr;
    matlabbatch{1}.spm.temporal.st.ta=matlabbatch{1}.spm.temporal.st.tr - ...
        (matlabbatch{1}.spm.temporal.st.tr/matlabbatch{1}.spm.temporal.st.nslices);
    switch slice_order
        case 1
            matlabbatch{1}.spm.temporal.st.so = 1:matlabbatch{1}.spm.temporal.st.nslices;
        case 2
            matlabbatch{1}.spm.temporal.st.so = matlabbatch{1}.spm.temporal.st.nslices:-1:1;
        case 3
            matlabbatch{1}.spm.temporal.st.so = [1:2:matlabbatch{1}.spm.temporal.st.nslices, ...
                2:2:matlabbatch{1}.spm.temporal.st.nslices];
        case 4
            matlabbatch{1}.spm.temporal.st.so = [2:2:matlabbatch{1}.spm.temporal.st.nslices, ...
                1:2:matlabbatch{1}.spm.temporal.st.nslices];
        case 5
            matlabbatch{1}.spm.temporal.st.so = [matlabbatch{1}.spm.temporal.st.nslices:-2:1, ...
                matlabbatch{1}.spm.temporal.st.nslices-1:-2:1];
        case 6
            matlabbatch{1}.spm.temporal.st.so = [matlabbatch{1}.spm.temporal.st.nslices-1:-2:1, ...
                matlabbatch{1}.spm.temporal.st.nslices:-2:1];
        case 7
            matlabbatch{1}.spm.temporal.st.so = load(slice_order_file,'-ascii');
        otherwise
            error('Slice order has to be 1 - 7')
    end;
    if (~exist(ref_slice, 'var'))
    matlabbatch{1}.spm.temporal.st.refslice = round(matlabbatch{1}.spm.temporal.st.nslices/2);
    else
        matlabbatch{1}.spm.temporal.st.refslice = ref_slice;
    end
    save(fullfile(fileparts(func_dir{1}),'01_slice_time_correct_job.mat'),'matlabbatch');
    spm('defaults', 'FMRI');
    
    disp(func_files{1}{1});
    spm_jobman('run', matlabbatch);
    prefix = strcat(prefix, matlabbatch{1}.spm.temporal.st.prefix);
    clear matlabbatch vs frames;
end

%% Realignment
if(motion_correct)
    load(fullfile(tooldir,'template_realign_job_spm12.mat')); % Load template
    
    for irun=1:length(func_files) % Each run
        matlabbatch{1}.spm.spatial.realign.estimate.data{irun} = [];
        vs=spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix)));
        [~,fname,~]=fileparts(vs(1).fname);
        for iframe=1:length(vs) %Each dynamic
            frames{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end
        if(size(frames, 2) > size(frames, 1))
            matlabbatch{1}.spm.spatial.realign.estimate.data{irun} = frames';
        else
            matlabbatch{1}.spm.spatial.realign.estimate.data{irun} = frames;
        end
        %         matlabbatch{1}.spm.spatial.realign.estimate.data{ir} = matlabbatch{1}.spm.spatial.realign.estimate.data{ir}';
        rp_file{irun}=fullfile(func_dir{irun}, strcat('rp_', fname,'.txt'));%If nuisance remove is selected
    end
    save(fullfile(fileparts(func_dir{irun}),'02_realign_job.mat'),'matlabbatch');
    
    disp(matlabbatch{1}.spm.spatial.realign.estimate.data{1}{1});
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
    for irun=1:length(func_files) % Each run
        fmri_plot_diff(QC_File, rp_file{irun}, motion_thresh, irun);
    end
end

%%Unbiased Anatomical
if(unbiased)
    load(fullfile(tooldir,'template_segment_job_unbiasedonly_spm12.mat'))
    matlabbatch{1}.spm.spatial.preproc.channel.vols=[];
    matlabbatch{1}.spm.spatial.preproc.channel.vols={anat_file};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,1')};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,2')};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,3')};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,4')};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,5')};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,6')};
    save(fullfile(fileparts(func_files{1}{1}),'03_unbiased_anat_job.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
    anat_file = fullfile(anat_dir, strcat('m', anat_file_nodir, '.nii'));
end

%% Coregistration
if(coregister)
    load(fullfile(tooldir,'template_coregister_job_spm8.mat')); % Load template
    [func_dir,func_file,ext]=fileparts(func_files{1}{1});
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1}=fullfile(func_dir,[prefix,func_file,ext]);
    %     matlabbatch{1}.spm.spatial.coreg.estimate.source{1}=anat_file;
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1}=fullfile(anat_dir, strcat('c1', anat_file_nodir, '.nii'));
    matlabbatch{1}.spm.spatial.coreg.estimate.other{1}=anat_file;
    save(fullfile(fileparts(func_files{1}{1}),'04_coregister_job.mat'),'matlabbatch');
    
    disp(func_files{1}{1});
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
    imgs = [];
    for irun=1:length(func_files)
        imgs = [{func_files{irun}{1}}]; %#ok<*CCAT1,AGROW>
    end;
    imgs = [imgs; {anat_file}];
    fmri_spm_check_reg(imgs, '-dpsc2', QC_File, '-append');
    imgs = [];
end;

%% Segmentation
if(segment_anat)
    load(fullfile(tooldir,'template_segment_job_spm12.mat')); % Load template
    matlabbatch{1}.spm.spatial.preproc.channel.vols=[];
    matlabbatch{1}.spm.spatial.preproc.channel.vols={anat_file};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,1')};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,2')};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,3')};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,4')};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,5')};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm={fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,6')};
    
    save(fullfile(fileparts(func_files{1}{1}),'05_segment_job.mat'),'matlabbatch');
    
    disp(func_files{1}{1});
    spm_jobman('run',matlabbatch);
    
    [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
    wm_file=fullfile(anat_dir,['wc2',anat_file_nodir,ext]);
    csf_file=fullfile(anat_dir,['wc3',anat_file_nodir,ext]);
    clear matlabbatch;
    imgs = [{fullfile(anat_dir, ['wc1',anat_file_nodir,ext])}, {wm_file}, {csf_file}];
    fmri_spm_check_reg(imgs, '-dpsc2', QC_File, '-append');
end


%% Spatial normalization
if(normalize==1) %just write - BROKEN (need to update the way that file list is generated)
    load(fullfile(tooldir,'template_normalize_job_spm12.mat')); % Load template
    n=1;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample=[];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def=[];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def={fullfile(anat_path, strcat('y_', anat_file_nodir, ext))};
    for irun=1:length(func_files) %Each run
        for i=1:length(func_files{irun}) %Each dynamic
            [func_dir,func_file,ext]=fileparts(func_files{irun}{i});
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n} = ...
                fullfile(func_dir,[prefix,func_file,ext]);
            n=n+1;
        end;
    end;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n} = ...
        anat_file;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n+1} = ...
        fullfile(anat_dir, ['c1',anat_file_nodir,'.nii']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n+2} = ...
        fullfile(anat_dir, ['c2',anat_file_nodir,'.nii']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n+3} = ...
        fullfile(anat_dir, ['c3',anat_file_nodir,'.nii']);
    [anat_dir,anat_file_nodir]=fileparts(anat_file);
    
    save(fullfile(fileparts(func_files{1}{1}),'06_normalize_job.mat'),'matlabbatch');
    
    disp(func_files{1}{1});
    spm_jobman('run',matlabbatch);
    prefix=['w',prefix];
    
    [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
    wimgs = [];
    %wimgs is for QC plotting
    for irun=1:length(func_files) %Each session
        [f_dir,f_file_nodir,f_ext]=fileparts(func_files{irun}{1});
        wimgs = [wimgs; ...
            {fullfile(f_dir, [prefix, f_file_nodir, f_ext])}]; %#ok<AGROW>
    end;
    % add anatomical
    wimgs = [wimgs; ...
        {fullfile(anat_dir, ['wc1' anat_file_nodir ext])}; ...
        {fullfile(fileparts(which('spm')),'toolbox', 'OldSeg', 'grey.nii')}];
    QC_spm_check_reg(wimgs, QC_File);
    clear matlabbatch;
    
elseif(normalize==2) %estimate and write using SPM's EPI template
    normPath = fullfile(fileparts(which('spm')), 'toolbox', 'OldNorm');
    clear vs frames
    for irun = 1:nruns
        vs=spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix)));
        for iframe=1:length(vs)
            frames{iframe}=strcat(vs(iframe).fname, ',', num2str(iframe));
        end;
        
        %use middle volume to estimate normalization parameters;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {strcat(vs(ceil(length(vs)/2)).fname, ',', num2str(ceil(length(vs)/2)))};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
        
        if(size(frames, 2) > size(frames, 1))
            matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = frames';
        else
            matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = frames;
        end
        
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {fullfile(normPath, 'EPI.nii,1')};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 45;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = normalize_bb;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = normalize_vox;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = normalize_prefix;
        
        save(fullfile(fileparts(func_dir{irun}), strcat('03_epi_normalize_job_run', num2str(irun), '.mat')),'matlabbatch');
        
        disp(vs(1).fname);
        spm_jobman('run',matlabbatch);
        
        %print out images of 1st & middle volume to check subject registration to MNI EPI & T1 template
        wimgs = [];
        template_EPI_file = fullfile(normPath, 'EPI.nii, 1');
        template_T1_file=fullfile(normPath, 'T1.nii,1');
        full_func_file1=strcat(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, [normalize_prefix prefix])), ',1');
        full_func_file2=strcat(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, [normalize_prefix prefix])), ',', num2str(ceil(length(vs)/2)));
        wimgs=[{full_func_file1};{full_func_file2};{template_EPI_file};{template_T1_file}];
        QC_spm_check_reg(wimgs, QC_File);
        clear job_file jobs matlabbatch vs frames;
        
    end %irun
    
    prefix=strcat(normalize_prefix, prefix);
    
end; %if normalize


%% Create brain_mask - This is possible only if we had segmented and is needed only if we run one of the following!! - FIX THIS
if(create_brain_mask_file==1)
    anat_file
    fmri_create_brain_mask(anat_file, brain_mask_file);
end;

disp(['Prefix ', prefix]);
for irun=1:length(func_files) %The following needs to be done session by session
    prefix_loop=prefix;
    disp(['Prefix_loop ', prefix_loop]);
    clear curr_func_files;
    %% HP filter or detrend
    if(hp_filter)
        vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix)));
        for iframe=1:length(vs)
            curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end;
        disp(curr_func_files{1});
        cutoffs = hp_filter;
        brain_mask_filename = brain_mask_file;
        drop_band = 0;
        if ~exist('band','var')
            % detrend is default option
            band = 5;
            drop_band = 1;
        end
        disp(['HP/Linear Detrend; band = ', num2str(band)])
        fmri_time_filt(curr_func_files, tr, hp_filter, band, brain_mask_file, low_ram);
        if (drop_band), clear band, end
        prefix_loop = strcat('f', prefix_loop);
        clear curr_func_files vs frames;
    end;
    
    %% Smooth the data - This is in case we want to run ICA
    if(ica_smooth)
        load(fullfile(tooldir,'template_smooth_job_spm12.mat')); % Load template
        
        vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix_loop)));
        
        for iframe=1:length(vs)
            curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end
        
        if(ica_smooth_fwhm==-1) %use 2*current voxel size, which =4 mm after normalization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ica_smooth_fwhm=round(abs(vs(1).mat(1,1).*2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if(length(ica_smooth_fwhm)==1)
            ica_smooth_fwhm=repmat(ica_smooth_fwhm,1,3);
        end
        matlabbatch{1}.spm.spatial.smooth.fwhm = ica_smooth_fwhm;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.data=curr_func_files';
        save(fullfile(fileparts(fileparts(func_files{1}{1})), strcat('07_ica_smooth_job_run', num2str(irun), '.mat')), 'matlabbatch');
        spm_jobman('run', matlabbatch);
        clear curr_func_files vs frames;
        %don't change prefix here - no further processing necessary on these files
    end;
    
    %% Estimate nuisances
    if(nuisance_estimate)
        [~, func_file, fext] = fileparts(func_files{irun}{1});
        
        %get list of volumes from which nuisances will be estimated and
        %removed
        vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix_loop)));
        for iframe=1:length(vs)
            curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end;
        
        % define .mat where estimated nuisance time series will be saved
        nuisance_file = fullfile(func_dir{irun}, strcat(prefix_loop, func_file, '_', nuisance_file_postfix, '.mat'));
        
        % Find realignment parameter file
        if(~motion_correct) %rp_file variable not defined if started script after motion correction
            if(~exist('rp_file', 'var')) %might have defined rp file in batch
                rplist = dir(fullfile(fileparts(curr_func_files{1}),'rp*.txt'));
                rpf = {rplist.name};
                if(length(rpf)>1)
                    warning('Multiple rp files') %sometimes rerun rigid body realignment after removing high motion volumes
                    disp(rpf)
                    which_file = input('Which rp file?');
                    rp_file_nm = rpf{which_file};
                else
                    rp_file_nm = rpf{1};
                end %length>1
                rp_file{irun} = fullfile(fileparts(curr_func_files{1}), rp_file_nm);
            end %~exist rp_file
        end %~motion_correct
        
        % if anat_file is not defined, use SPM tissue priors to create
        % white matter and CSF masks
        if ~exist('wm_file','var') && ~exist('anat_file', 'var')
            which_spm = fileparts(which('spm'));
            wm_file = fullfile(which_spm, 'toolbox', 'OldSeg', 'white.nii');
        elseif ~exist('wm_file','var')
            [anat_dir,anat_file_nodir,ext] = fileparts(anat_file);
            wm_file = fullfile(anat_dir, strcat('wc2', anat_file_nodir, ext));
        end
        
        if nr_options.wm_signal == 0
            wm_file='';
        end
        
        if ~exist('csf_file','var') && ~exist('anat_file', 'var')
            csf_file = fullfile(which_spm, 'toolbox', 'OldSeg', 'csf.nii');
        elseif ~exist('csf_file','var')
            csf_file = fullfile(anat_dir,['wc3',anat_file_nodir,ext]);
        end
        
        if nr_options.csf_signal== 0
            csf_file='';
        end
        if ~exist('npref','var')
            npref = 'n';
        end
        fmri_extract_nuisance(curr_func_files, tr, rp_file{irun}, wm_file, csf_file, brain_mask_file, nr_options, nuisance_file);
        
    end %if nuisance_estimate
    
    if(nuisance_regress)
        if(~nuisance_estimate) %if nuisance_estimate, already have a list of curr_func_files & nuisance_file defined
            vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix_loop)));
            for iframe=1:length(vs)
                curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
            end
            nuisance_file = fullfile(func_dir{irun}, strcat(prefix_loop, func_file, '_', nuisance_file_postfix, '.mat'));
        end
        
        if(low_ram~=1)
            [npref] = fmri_regress_nuisance(curr_func_files, brain_mask_file, nuisance_file, nr_options);
            
        else
            fmri_regress_nuisance_1D(curr_func_files,brain_mask_file,nuisance_file, nr_options);
            
        end
        prefix_loop=[npref,prefix_loop];
        clear curr_func_files vs nuisance_file;
    end; %if nuisance_regress
    
    %% Smooth the residuals
    if(rs_smooth)
        load(fullfile(tooldir,'template_smooth_job_spm12.mat')); % Load template
        vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix_loop)));
        for iframe=1:length(vs)
            curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end
        if(smooth_fwhm==-1), smooth_fwhm=round(abs(vs.mat(1,1).*2)); end;
        if(length(smooth_fwhm)==1)
            smooth_fwhm=repmat(smooth_fwhm,1,3);
        end
        matlabbatch{1}.spm.spatial.smooth.fwhm = smooth_fwhm;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.data=curr_func_files';
        save(fullfile(fileparts(fileparts(func_files{1}{1})),'08_smooth_residuals_job.mat'), 'matlabbatch');
        spm_jobman('run', matlabbatch);
        prefix_loop=['s',prefix_loop];
        clear curr_func_files vs;
    end
    
    %% Time domain filtering
    if(bp_filter)
        vs = spm_vol(strrep(func_files{irun}{1}, strcat('run-', num2str(irun, '%02d'), filesep), strcat('run-', num2str(irun, '%02d'), filesep, prefix_loop)));
        for iframe=1:length(vs)
            curr_func_files{iframe} = strcat(vs(iframe).fname, ',', num2str(iframe));
        end
        fmri_time_filt(curr_func_files,tr,bp,3,brain_mask_file,low_ram);
        disp('BP')
        prefix_loop=['f',prefix_loop];
        clear curr_func_files vs
    end
    
end

if(slice_time_correct || motion_correct || normalize)
    % Make QC File a pdf
    if(~ispc && ~ismac)
        system(sprintf('ps2pdfwr "%s" "%s"', QC_File, ...
            strrep(QC_File, '.ps', '.pdf')));
    end
end

end


