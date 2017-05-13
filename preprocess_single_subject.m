function [ ] = preprocess_single_subject(ID, func_files, sess_names, procdir, tooldir, rename_task, specs_file, prefixes)
%Inputs
%   ID: subject identification number padded with zeros to ensure equal length in format sub-ID
%   func_files: cell array of path/filename(s) for each raw run to be processed
%       from given subject in the order they were collected
%   sess_names: cell array of same size as func_files indicating session
%       names for each run
%   procdir: directory where raw data will be copied and processed; script
%       creates ID folder within procdir if it does not already exist
%   tooldir: directory where processing scripts are located; default
%       assumes fmri_preproc_toolbox exists within present working directory
%   rename_task: e.g., task-rest_bold or task-GNG_bold; script will copy
%       and rename func_files to procdir/ID using the following convention -
%       sub-ID/sess_name/run-?/sub-ID_sess_name_rename_task_run-?,
%       where run number is determined by the order of func_files
%   specs_file: path/file.m containing processing preferences; default:
%       fmri_preprocess_par2ica_spm12.m
%NOTE: assumes 4D files

%check for missing input arguments
if(~exist('ID', 'var'))
    error('Need to specify ID number')
end

nruns = length(func_files);
if ~isequal(nruns, length(sess_names))
    error('lengths of func_files and sess_names do not match')
end

if(~exist('procdir', 'var'))
    error('Must specify procdir: directory where raw data will be copied and processed')
end

if(~exist('tooldir', 'var'))
    tooldir = fullfile(pwd, 'fmri_preproc_toolbox');
end

if(~exist('rename_task', 'var'))
    rename_task = 'task-rest_bold';
end

if(~exist('specs_file', 'var'))
    specs_file = 'fmri_preprocess_specs_par2ica_spm12.m';
end

if(~exist('prefixes', 'var'))
    prefixes = cell(nruns, 1);
end

cwd = pwd;

addpath(tooldir)

disp(strcat('Working on ID: ', ID))

%folder where subject's data will be processed
sub_path = fullfile(procdir, ID);

%check if this folder exists; if not, create it
if(~exist(sub_path, 'dir'))
    [SUCCESS, MESSAGE, MESSAGEID] = mkdir(sub_path);
end

func_dir = cell(nruns, 1);
sess_dir = fullfile(procdir, ID, sess_names{1}{1});

%% Make a copy of raw data first
for irun = 1:nruns
    
    %where to copy raw data before processing it
    func_dir{irun} = fullfile(procdir, ID, sess_names{irun}{1}, 'func', strcat('run-', num2str(irun, '%02d')));
    
    %if ftarget doesn't exist, create it
    if(~exist(func_dir{irun}, 'dir'))
        [SUCCESS, MESSAGE, MESSAGEID] = mkdir(func_dir{irun});
    end
    
    %figure out raw data file type (assumes 4D)
    [fSource, sourceName, sourceExt] = fileparts(func_files{irun}{1});
    
    if(sourceExt=='.rec')
        %copy functional par & rec to func_dir
        func_files{irun}{1}=fullfile(func_dir{irun}, strcat(sourceName, '.rec'));
        [success, message, ~] = copyfile(fullfile(fSource, strcat(sourceName, '.rec')), func_files{irun}{1});
        [success, message, ~] = copyfile(fullfile(fSource, strcat(sourceName, '.par')), fullfile(func_dir{irun}, strcat(sourceName, '.par')));
        nii_names{irun} = strcat(ID, '_', sess_names{irun}{1}, '_', rename_task, '_run-', num2str(irun, '%02d'), '.nii');
        
    else
        func_files{irun}{1} = fullfile(func_dir{irun}, strcat(prefixes{irun}, ID, '_', sess_names{irun}{1}, '_', rename_task, '_run-', num2str(irun, '%02d'), '.nii'));
        [success, message, ~] = copyfile(fullfile(fSource, strcat(sourceName, '.nii')), func_files{irun}{1});
    end %if .rec
    
end; %irun

%specify path/name.ps of file to print quality info to
QC_File = fullfile(sess_dir, strcat(ID, '_', sess_names{irun}{1}, '_', rename_task, '_QC.ps'));

%load preferences into matlab memory
run(specs_file)

%save current matlab workspace
[~, sname, ~] = fileparts(specs_file);
setup_file=fullfile(sess_dir, ...
    strcat(ID, '_', sess_names{irun}{1}, '_', rename_task, '_', sname, '_', datestr(now,'yyyymmdd'), '.mat'));
save(setup_file);

cd(sess_dir)

%process data
fmri_preprocess_spm12(setup_file);
