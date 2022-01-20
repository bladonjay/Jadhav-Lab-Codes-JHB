%JHB_data_preprocessing
%% this is a general script to explore data preprocessing and then double check the nueral data


% this will get you spike and LFP data
% run_ML_jay


%% load spike data

load('J:\My Drive\Work\XFB3_04_210827_spikes.mat')
load('D:\My Drive\Work\XFB3_04_210827_spikes.mat')

units=units([units.noise_overlap]<.015 & [units.isolation]>.95); % only take well isolated units

%% load LFP data

addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\TrodesFullPackage\Trodes_2-2-3_Windows64'));
[lfpRaw]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\XFB3_03_20210827.LFP\XFB3_03_20210827.LFP_nt9ch1.dat');
rmpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\TrodesFullPackage\Trodes_2-2-3_Windows64'));


addpath('E:\GithubCodeRepositories\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab')
[lfpRaw]=readTrodesExtractedDataFile('J:\My Drive\Work\XFB3_03_20210827.LFP_nt9ch1.dat');
rmpath('E:\GithubCodeRepositories\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab');

lfpData=struct;
lfpData.descript=lfpRaw.description;
lfpData.voltageScaling=lfpRaw.voltage_scaling;
lfpData.original_file=lfpRaw.original_file;
lfpData.data=double(lfpRaw.fields.data);
lfpData.rawTime=[1:length(lfpData.data)]+lfpRaw.first_timestamp;
lfpData.startBin=lfpRaw.first_timestamp;
lfpData.sampRate=1500; % default
lfpData.starttime=lfpRaw.first_timestamp/1500;
lfpData.samprate=lfpData.sampRate;
lfpData.realTime=lfpData.rawTime/lfpData.sampRate; % i think its 1.5khz...

tet9=lfpData;
%tet39=lfpData;


%% make sure we are seeing spikes, ripples, theta

% this will generate a spectrogram, a ripple band filtered signal above a
% theta filtered signal, and then spikes below.  This will be for the whole
% session but you can zoom in in the x axis

% so 128/1500=85 msec window with a roughly 40 msec crawl
[s,w,t]=spectrogram(tet9.data,100,25,1:2:200,1500,'yaxis','power');

figure; 
sp=subplot(3,1,1);
imagesc(t+lfpData.starttime,w,SmoothMat2(zscore(log(abs(s)),1,2),[10 10])); set(gca,'ydir','normal');
box off; ylabel('Frequency, Hz'); set(gca,'XTick',[])

sp(2)=subplot(6,1,3);
% first theta power
[~,thetaAmp,thetaFilt] = GetLFPBand(tet9.data,tet9.rawTime); % defaults to theta
load('ripplefilter.mat');
ripEEG=filtereeg2(tet9,ripplefilter,'int16',0);
%{
plot(tet9.realTime,thetaFilt); ylabel('Theta filtered signal');
yyaxis right;
plot(tet9.realTime,ripEEG.data(:,1)); ylabel('ripple filtered signal')
legend('Theta (5-12 Hz)','Ripple (150-200 hz)')
%}

plot(tet9.realTime,tet9.data);
box off; axis off;
sp(3)=subplot(6,1,4);
plot(tet9.realTime,zscore(thetaAmp)); hold on;
plot(tet9.realTime,zscore(smoothdata(ripEEG.data(:,3),20)),'k')
box off; ylim([-2 10])

sp(4)=subplot(3,1,3);
% load spike data

% now go for all of them;
mycolors=parula(length(units));

% now go through and create units;


%% need to parse pfc vs hpc units (probably two color scales, red and blue

it=0.5;
for i=1:length(units)
    okspikes=(units(i).spikes(:,2)>lfpData.realTime(1) & units(i).spikes(:,2)<lfpData.realTime(end));
    if sum(okspikes)>50 && units(i).meanrate<20
        newx=repmat(units(i).spikes(okspikes,2)',3,1);
        newy=repmat([it it+1 nan]',sum(okspikes),1);
        plot(newx(:),newy(:),'Color',mycolors(i,:),'LineWidth',2);
        hold on;
        it=it+1;
    end
end
linkaxes(sp,'x'); box off

