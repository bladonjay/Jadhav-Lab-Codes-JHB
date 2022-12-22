%% dataPreprocessing SOP
%{

Welcome to preprocessing

This will outline the steps required to preprocess data in broad sweeps,
there are nested scripts and functions in this SOP that will guide you
towards details


The general outline of data preprocessing:
A. Behavioral data processing
    1. you can pull data from any folder, the filenames must follow the
        convention: logMM-DD-YYYY(1-Rat1-Rat2.stateScriptLog, an example is...
        log10-12-2020(1-XFB1-XFB2).stateScriptLog (from the date 10/12/2022,
        run 1, rats XFB1 (on 123 side) and XFB2(on abc side).
    2. run SocialTaskBehaviorAnalysis script to generate dataset of n session
        struct that will have fields like file, folder, rat data, date, run
        number, and rat names


B. Neural Data processing
    1. organize your files and folders in a very specific manner
        some of this is in run_ML_Jay and in JHB_neural_data_preprocessing

    2. extract spike data from one DAY of data across recording epochs using
        MountainSort
        This is in run_ML_Jay
        This is validated with JHB_neural_data_preprocessing

C. Camera Data
    1. convert .h264 files into mp4 for fast reading using ffmpeg
        use fileRenamingScript for now, but this will change, this includes
        ffmpegtranscode, a function that requires an ffmpeg download adn
        the sdk download in matlab

    2. extract position data from videos via DeepLabCut
       DeepLabCut can be run from matlab scripting but i need to build it...
        this converts fileformats without compression:

    3. make sure all of those datasets have the same clock (havent done
    yet)

    4. you have a dataset now!
%}
%% Pre-workup. organizing file names and generating file conventions


% the folder convention for neural data example is here:
% for neural data there are two conventions, and i had to switch so both
% are still around...
% my convention:
% socialTask (folder)
%   XFB3(animal raw data)
%       8-27-2021 (1st recording day)
%           JHB_XFB3_8-27-2021 (run folder)
%               JHB_XFB3_8-27-2021_run2.rec (rec file for run 2)
%
% Jadhav convention:
%           
%   XFB3_direct (Processed data)
%       XFB3_04_20210830 (day 4)
%           XFB3_04_20210830_run2 (run 2 on that day)
%           XFB3_04_20210830_Session.mat % full data
%           XFB3_04_20210830_spikes.mat
%           Mountainsort (MS data, from which offlinesorter data was
%           generated)
%           OfflineSorter (ofs data, from which spikes.mat was generated)


% the folder convention for behavioral data only is here:
% BehaviorOnly/completeDatasets
%   CohortXFB
%       6-4-2021
%           log06-04-2021(1-XFB1-XFB2).stateScriptLog
%           log06-04-2021(1-XFB1-XFB2).mp4
%% A. Behavioral data Processing

edit SocialTaskBehaviorAnalysis


%% B. Neural Data processing



%
%
% for spike data use gatherMS_to_FFunits in parent folder
%
%

edit gatherMS_to_FFunits.m

edit JHB_neural_data_preprocessing

cd Preprocessing
edit run_ML_Jay.m


%% C. 

%% 2. extract spike data, LFP from one day of data
cd 'E:\GithubCodeRepositories\Jadhav-Lab-Codes\FMR1\Preprocessing'
edit run_ML_Jay.m

cd 'E:\GithubCodeRepositories\Jadhav-Lab-Codes\FMR1'
edit JHB_neural_data_preprocessing





%% parsing behavioral data

% lets start with CohortB because we have neural data here


cohortDir=uigetdir();
dayDirs=dir(cohortDir);
dayDirs=dayDirs(3:end);
% i'm going to bracket it to between 6-2 and 6-30 2021

datePat='(?<month>[0-9]+)-(?<day>[0-9]+)-(?<year>[0-9]+)';
for i=1:length(dayDirs)
    newdate=dayDirs(i).name;
    newdate(newdate=='-')='/';
    realDate=regexp(dayDirs(i).name,datePat,'names');
    dayDirs(i).datenum=sprintf('%02d%02d%02d',str2double(realDate.year(end-2:end)),str2double(realDate.month),str2double(realDate.day));
end

dayDirs([dayDirs.datenum]<datenum('6/1/2021') | [dayDirs.datenum]>datenum('6/30/2021'))=[];