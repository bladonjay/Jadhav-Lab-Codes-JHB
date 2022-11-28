function [filtered]= jhb_LFPdayprocess(rawDir, dataDir, animID, sessionNum, varargin)
% This function creates theta filtered EEGs from non- referenced EEG files.

%rawDir -- the directory where the raw dat folders are located
%dataDir -- the directory where the processed files should be saved
%animID -- a string identifying the animal's id (appended to the
%beginning of the files)
%sessionNum -- the session number (in chronological order for the animal)
%
% Varargin Options
%-----------------
%		'f', matfilename
%			specifies the name of the mat file containing the
%			theta filter to use
%			(default filters/thetafilter.mat).
%			Note that the filter must be called 'thetafilter'.
%		'ref', 0 or 1
%			specifies whether to use eegref (1) or eeg (0)	
% TO DO: Add ability to only run thetadayprocess on specific tetrodes,
% currently runs through all tetrodes and epochs on a given day
% RN EDIT 8/8/17: Skips filtering if file already exists
% RN EDIT 11/7/17: Accepts varargin 'ref' can be set to 1 or 0. If 1 uses eegref and savefile is gammaref. If 0 uses eeg and save is gamma.
 
% Default for theta filtering is to apply filter to eeg (referenced to ground)
eegStr = 'eeg';
 
% Default theta filter
f = [fileparts(mfilename('fullpath')) filesep 'filters/thetafilter.mat'];
band='theta';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'f'
            f = varargin{option+1};
        case 'ref'
             if varargin{option+1}==1
                 eegStr = 'eegref';
             elseif varargin{option+1}==0
                 eegStr = 'eeg';
             else
                 error('Invalid ref option. Must be 0 or 1')
             end
        case 'band'
            band=varargin{option+1};
    end
end

if isempty(f)
    error('Filter not found!');
else
    myFilt=load(f);
end
filtField=fieldnames(myFilt);
filterSpec=myFilt.(filtField{1});

currDir = pwd;
cd(dataDir);

if isempty(dir(['*EEG']))
    error('No EEG folder found!');
end

cd('EEG');
sString= sprintf('%02i',sessionNum);

relFiles= dir([animID eegStr sString '*.mat']);
for i=1:length(relFiles)

    x = sscanf(relFiles(i).name,[animID eegStr '%d-%d-%d.mat']);
    e = x(2);
    t = x(3);
    saveFile = strrep(relFiles(i).name,'eeg',band);
    if ~isempty(dir(saveFile))
        fprintf('LFP filtered file already exists for day %02d. Skipping...\n',sessionNum);
        continue;
    else
        fprintf('%s Filtering LFP for %02i-%02i-%02i\n',band,sessionNum,e,t);
    end

    eeg=load(relFiles(i).name);
    eeg=eeg.(eegStr);
    filtered{sessionNum}{e}{t} = filtereeg2(eeg{sessionNum}{e}{t}, filterSpec, 'int16', 1);
    filtered{sessionNum}{e}{t}.voltage_scaling    =eeg{sessionNum}{e}{t}.voltage_scaling;
    filtered{sessionNum}{e}{t}.low_pass_filter    =eeg{sessionNum}{e}{t}.low_pass_filter;
    filtered{sessionNum}{e}{t}.referenced         =eeg{sessionNum}{e}{t}.referenced;
    filtered{sessionNum}{e}{t}.data_voltage_scaled=eeg{sessionNum}{e}{t}.data_voltage_scaled;
    filtered{sessionNum}{e}{t}.nTrodeChannel      =eeg{sessionNum}{e}{t}.nTrodeChannel;
    filtered{sessionNum}{e}{t}.nTrode             =eeg{sessionNum}{e}{t}.nTrode;
    filtered{sessionNum}{e}{t}.endtime            =eeg{sessionNum}{e}{t}.endtime;
    filtered{sessionNum}{e}{t}.clockrate          =eeg{sessionNum}{e}{t}.clockrate;
    filtered{sessionNum}{e}{t}.timerange          =eeg{sessionNum}{e}{t}.timerange;

    % save the resulting file
    if strcmp(eegStr,'eegref')
        filteredref = filtered;
        save(saveFile,'filteredref');
    else
        save(saveFile,'filtered');
    end
    clear thetaref theta eeg
end
cd(currDir)