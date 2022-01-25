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
% with each lfp file you have to get the timestamps too
[lfpRaw]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\03_20210827.LFP\03_20210827.LFP_nt9ch1.dat');

[lfptimestamps]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\03_20210827.LFP\03_20210827.timestamps.dat');
[spikebandtimestamps]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\XFB3_03_210827.time.dat');



[lfpRaw]=readTrodesExtractedDataFile('D:\My Drive\Work\XFB3_03_20210827.LFP_nt9ch1.dat');
rmpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\TrodesFullPackage\Trodes_2-2-3_Windows64'));


addpath('E:\GithubCodeRepositories\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab')
[lfpRaw]=readTrodesExtractedDataFile('J:\My Drive\Work\XFB3_03_20210827.LFP_nt9ch1.dat');
rmpath('E:\GithubCodeRepositories\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab');

lfpData=struct;
lfpData.descript=lfpRaw.description;
lfpData.voltageScaling=lfpRaw.voltage_scaling;
lfpData.original_file=lfpRaw.original_file;
lfpData.data=double(lfpRaw.fields.data);
lfpData.rawTime=lfptimestamps.fields.data;
lfpData.startBin=lfpRaw.first_timestamp;
lfpData.sampRate=1500; % default
lfpData.starttime=lfpRaw.first_timestamp/1500;
lfpData.samprate=lfpData.sampRate;
lfpData.realTime=lfpData.rawTime/lfptimestamps.clockrate; % i think its 1.5khz...

tet9=lfpData;
%tet39=lfpData;


%% make sure we are seeing spikes, ripples, theta

% this will generate a spectrogram, a ripple band filtered signal above a
% theta filtered signal, and then spikes below.  This will be for the whole
% session but you can zoom in in the x axis

% so 128/1500=85 msec window with a roughly 40 msec crawl
%[s,w,t]=spectrogram(tet9a.data,100,25,1:2:200,1500,'yaxis','power');

figure; 
%sp=subplot(3,1,1);
%imagesc(t+lfpData.starttime,w,SmoothMat2(zscore(log(abs(s)),1,2),[10 10])); set(gca,'ydir','normal');
%box off; ylabel('Frequency, Hz'); set(gca,'XTick',[])

sp(2)=subplot(6,1,3);
% first theta power
[~,thetaAmp,thetaFilt] = GetLFPBand(tet9.data,tet9.rawTime,[5 12],0); % defaults to theta
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
ylabel('Raw Ca1 LFP');
sp(3)=subplot(6,1,4);
plot(tet9.realTime,zscore(thetaAmp)); hold on;
plot(tet9.realTime,zscore(smoothdata(ripEEG.data(:,3),20)),'k')
box off; ylim([-2 10]); axis off;
legend('Theta power','Ripple Power');
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
linkaxes(sp,'x'); box off;

%% testing my timing issue

[lfptimeRaw]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\03_20210827.LFP\03_20210827.timestamps.dat');


[spikeRaw]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\03_20210827.time\03_20210827.continuoustime.dat');
[epochData]=readTrodesExtractedDataFile('F:\SocialData\Neural\XFB3\03_20210827\03_20210827.time\03_20210827.time.dat');



epochTS=epochData.fields.data;
lfpTimestamps=lfpRaw.fields.data;
spikeTimestamps=spikeRaw.fields.data;

% my guess is that it is handling the recording lags inproperly.

spikeEpoch=spikeTimestamps(spikeTimestamps>epochTS(1) & spikeTimestamps<epochTS(1,2));
lfpEpoch=lfpT

figure; plot(lfpTimestamps);
hold on;
plot(spikeTimestamps);

%%

figure; subplot(2,1,1);
plot(tet9.rawTime,tet9.data);
subplot(2,1,2);


it=0.5;
for i=1:length(units)
    okspikes=(units(i).spikes(:,1)>lfpData.realTime(1) & units(i).spikes(:,1)<lfpData.rawTime(end));
    if sum(okspikes)>50 && units(i).meanrate<20
        newx=repmat(units(i).spikes(okspikes,1)',3,1);
        newy=repmat([it it+1 nan]',sum(okspikes),1);
        plot(newx(:),newy(:),'Color',mycolors(i,:),'LineWidth',2);
        hold on;
        it=it+1;
    end
end





%% plotting units only
figure;
it=0.5;

mycolors=parula(length(units));

for i=1:length(units)
   okspikes=(units(i).spikes(:,2)>5400 & units(i).spikes(:,2)<6800);
    if sum(okspikes)>50 && units(i).meanrate<20
        newx=repmat(units(i).spikes(okspikes,2)',3,1);
        newy=repmat([it it+1 nan]',sum(okspikes),1);
        plot(newx(:),newy(:),'Color',mycolors(i,:),'LineWidth',2);
        hold on;
        it=it+1;
    end
end
xlabel('Time, seconds'); ylabel('Unit number');
box off;








