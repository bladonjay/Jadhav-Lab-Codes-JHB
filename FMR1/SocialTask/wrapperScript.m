%% dataPreprocessing SOP
%{

Welcome to preprocessing

This will outline the steps required to preprocess data in broad sweeps,
there are nested scripts and functions in this SOP that will guide you
towards details


The general outline of data preprocessing:

1. organize your files and folders in a very specific manner

2. extract spike data from one DAY of data across recording epochs using
MountainSort

3. extract position data from videos via DeepLabCut

4. extract event data from the DIO's or the statescript ledger using
probably custom code

5. make sure all of those datasets have the same clock

6. you have a dataset now!
%}



%% 1. organizing file names and generating file conventions


% the folder convention example is here:
% socialTask (folder)
%   XFB3(animal raw data)
%       01_210105 (1st recording on 01/05/2021)
%           
%   XFB3_direct (Processed data)
%       


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