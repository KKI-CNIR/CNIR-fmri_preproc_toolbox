%procdir: directory where raw data will be copied and processed
%tooldir: directory where processing scripts are located

if(ismac)
    frawdir = '/Volumes/fmri/DATA/Children/3T';
    procdir = '/Users/nebel/Projects/preprocessing/processed_data';
    tooldir = '/Users/nebel/Projects/preprocessing/fmri_preproc_toolbox';
elseif(isunix)
    addpath('/matlab/spm12/');
    frawdir = '/mnt/data/fmri/DATA/Children/3T';
    procdir = '/mnt/bisivity/fmri/RestingState/Kids/rsData';
    tooldir = '/mnt/bisivity/fmri/RestingState/Jobs/sj_scripts/fmri_preproc_toolbox';
else
    frawdir = '\\kki-gspnas1\DCN_data$\fmri\DATA\Children\3T';
    procdir = '\\kki-gspnas1\LNIR_DATA$\fmri\RestingState\Kids\rsData';
    tooldir = '\\kki-gspnas1\LNIR_DATA$\fmri\RestingState\Jobs\sj_scripts\fmri_preproc_toolbox\';
end

%slist: path/name.txt containing list of IDs to be processed
%istart: line of slist to start on
%task: what to look for in each subject's raw data folder
%rename_task: e.g., task-rest_bold or task-GNG_bold; processing script will copy
%       and rename func_files to procdir/ID using the following convention -
%       sub-ID/sess-DateOfSession/func/run-??/sub-ID_sess-DateOfSession_rename_task_run-??.nii,
%       where run number is determined by the order of func_files
%specs_file: path/file.m containing processing preferences

slist = fullfile(tooldir, 'subjects2process.txt');
istart = 3;
iend = 9;
task = '*RestState*';
rename_task = 'task-rest_bold';
specs_file = 'fmri_preprocess_specs_par2cluster_spm12.m';

addpath(tooldir)

%start with raw data
prefixes = {''};

%open slist and get functional IDs (SM numbers)
fid = fopen(slist);
subjects = textscan(fid, '%s');
fclose(fid);
SMs = strtrim(subjects{1});

if(~exist('iend', 'var'))
iend = length(SMs);
end

cwd = pwd;

%keep track of IDs in list missing data
fid = fopen(fullfile(cwd, 'ID_DOS_missing_parrec.txt'), 'w');

%loop over subjects in slist; figure out what raw files to use for each ID
for isub = istart:iend
    rawSM = SMs{isub};
    
    %add leading zeros to make all IDs the same length
    ID = strcat('sub-', num2str(str2num(rawSM), '%04d'));
    
    disp(strcat('Working on ID: ', ID))
    
    %define folders where subject's raw data are stored
    fSource = fullfile(frawdir, rawSM);
    
    %get list of directories in fSource
    fdates = dir(fSource);
    dirs = [fdates.isdir]';
    dates = {fdates(dirs).name}';
    
    %find which subdirectories have dates as names
    dates = regexp(dates, '^[0-9]*_[0-9]*_[0-9]*', 'match');
    dates(cellfun('isempty', dates)) = [];
    
    %check # of folders with dates as names
    if(length(dates) > 1)
        warning('Multiple scan dates')
        disp(dates)
        disp_str = sprintf('Which DOS? Enter 1 for %s, 2 for %s, etc:', char(dates{1}), char(dates{2}));
        which_file = input(disp_str);
        fdate = char(dates{which_file});
    else
        fdate = char(dates{1}');
    end %more than one scan date
    
    fprintf('Using data from %s\n', fdate)
    
    %get list of folders from DOS with task in their names
    tasklist = dir(fullfile(fSource, fdate, task));
    n_fold = length(tasklist);
    
    if(isempty(tasklist))
        warning(strcat('No ', task,' data on ', char(match)))
        disp(dates)
        disp_str = sprintf('Try a different date? Enter 0 to quit');
        which_fdate = input(disp_str);
        if(which_fdate)
            fdate = dates{which_fdate};
            tasklist = dir(fullfile(fSource, match, task));
            n_fold = length(tasklist);
        else
            error('No data was collected on that date')
        end
    end %no task folders on DOS
    
    if(n_fold > 1)
        warning('Multiple task folders')
        disp(tasklist(:).name)
        which_file = input('Which task?');
        taskfold = tasklist(which_file).name;
    else
        taskfold = tasklist.name;
    end
    
    %convert functional date from MM_DD_YYYY to YYYYMMDD format
    matchStr = regexp(fdate, '_', 'split');
    DOS = strcat('ses-', matchStr{3}, num2str(str2num(matchStr{1}), '%02d'), num2str(str2num(matchStr{2}), '%02d'));
    
    %get list of .rec files in task folder that match the task name
    freclist = dir(fullfile(fSource, fdate, taskfold, strcat(rawSM, '*.rec')));
    %# of functional runs
    nruns = length(freclist);
    fprintf('Found %s run(s) matching %s\n', num2str(nruns), task)
    
    if(nruns)
        runs = {freclist(:).name}';
        
        prefixes = repmat(prefixes, nruns, 1);
        
        %initialize empty cell variables; format of cells is holdover from
        %when data saved as collection of 3D files
        func_files = cell(nruns, 1);
        sess_names = cell(nruns, 1);
        
        for irun = 1:nruns
            
            func_files{irun}{1} = fullfile(fSource, fdate, taskfold, runs{irun});
            sess_names{irun}{1} = DOS;
            
        end; %irun
        
        preprocess_single_subject(ID, func_files, sess_names, procdir, tooldir, rename_task, specs_file)
        
        cd(cwd)
        
        clear QC_File func_files sess_names full_func_file1 full_func_file2;
        
    else %no task.rec from DOS; needs to be copied from godzilla
        warning(strcat('Missing ', task,' .rec from ', DOS))
        fprintf(fid, '%s %s\n', rawSM, DOS);
        continue
    end
    
    fclose(fid);
    
end

