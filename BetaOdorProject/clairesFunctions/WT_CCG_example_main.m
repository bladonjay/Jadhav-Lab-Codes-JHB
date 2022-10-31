%clc;
%clear all;
%close all;
%%
%-----------load the data----------%
%load('Spiketrains_running.mat')
%load('Spiketrains_SWRs.mat')
%%
%-----------CCG params----------%
% for theta CCG
bin = 0.01; % 10 ms
tmax = 0.5; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln. 
% for SWR CCG
bin_swr = 0.002; % 2 ms
sw1_swr = 0.01; % 10 ms smoothing
%%
%-----------Theta CCG----------%
[timebase_theta, rawcorr_theta, corr_sm_theta] = spiketrainxcorr(spikes1,spikes2,bin,tmax,sw1);
set(0,'defaultaxesfontsize',16);
figure('color','w')
plot(timebase_theta.*1000,corr_sm_theta,'linewidth',2)
hold on
plot([0,0],[0,0.35],'k--')
xlabel('Time lag (ms)')
ylabel('Correlation')
title('Crosscorrelogram during active running')
box off
%%
%------convert raster to histogram-----%
Cell1_awake_hist=rast2mat(Cell1_SWRawake_raster);
Cell2_awake_hist=rast2mat(Cell2_SWRawake_raster);
Cell1_sleep_hist=rast2mat(Cell1_SWRsleep_raster);
Cell2_sleep_hist=rast2mat(Cell2_SWRsleep_raster);
awakeSWRnum = length(Cell1_awake_hist(:,1));
sleepSWRnum = length(Cell1_sleep_hist(:,1));
%%
%-----------Raster Plots----------%
figure('color','w','position',[100 100 500 500])
for i = 1:awakeSWRnum
    subplot(221)
    plot(Cell1_SWRawake_raster{i},i.*ones(1,length(Cell1_SWRawake_raster{i})),...
        's','markersize',2,'MarkerEdgeColor','none','MarkerFaceColor','r')
    hold on
    plot([0,0],[0,awakeSWRnum],'--','color','cyan','linewidth',2)
    axis([-500 500 0 awakeSWRnum])
    xlabel('Time from SWR onset (ms)')
    ylabel('SWR number')
    title('Cell 1 (awake)')
    box off
    subplot(223)
    plot(Cell2_SWRawake_raster{i},i.*ones(1,length(Cell2_SWRawake_raster{i})),...
        's','markersize',2,'MarkerEdgeColor','none','MarkerFaceColor','r')
    hold on
    plot([0,0],[0,awakeSWRnum],'--','color','cyan','linewidth',2)
    axis([-500 500 0 awakeSWRnum])
    xlabel('Time from SWR onset (ms)')
    ylabel('SWR number')
    title('Cell 2 (awake)')
    box off
end
for i = 1:sleepSWRnum
    subplot(222)
    plot(Cell1_SWRsleep_raster{i},i.*ones(1,length(Cell1_SWRsleep_raster{i})),...
        's','markersize',2,'MarkerEdgeColor','none','MarkerFaceColor','k')
    hold on
    plot([0,0],[0,sleepSWRnum],'--','color','cyan','linewidth',2)
    axis([-500 500 0 sleepSWRnum])
    xlabel('Time from SWR onset (ms)')
    ylabel('SWR number')
    title('Cell 1 (sleep)')
    box off
    subplot(224)
    plot(Cell2_SWRsleep_raster{i},i.*ones(1,length(Cell2_SWRsleep_raster{i})),...
        's','markersize',2,'MarkerEdgeColor','none','MarkerFaceColor','k')
    hold on
    plot([0,0],[0,sleepSWRnum],'--','color','cyan','linewidth',2)
    axis([-500 500 0 sleepSWRnum])
    xlabel('Time from SWR onset (ms)')
    ylabel('SWR number')
    title('Cell 2 (sleep)')
    box off
end
%%
%-----------PSTH Plots----------%
% bin into 10-ms 
Cell1_awake_hist_bin = squeeze(sum(reshape(Cell1_awake_hist,awakeSWRnum,10,100),2)).*10;
Cell2_awake_hist_bin = squeeze(sum(reshape(Cell2_awake_hist,awakeSWRnum,10,100),2)).*10;
Cell1_sleep_hist_bin = squeeze(sum(reshape(Cell1_sleep_hist,sleepSWRnum,10,100),2)).*10;
Cell2_sleep_hist_bin = squeeze(sum(reshape(Cell2_sleep_hist,sleepSWRnum,10,100),2)).*10;

%smooth window
nstd=2;
x = 1:nstd+1;
mu = mean(x);
g0 = exp( (-(x-mu).^2)./(2*(nstd^2)));
s = sum(g0);
g1 = g0/s;
filterlen = length(g1);
startp = round(filterlen / 2);
timevec = -495:10:495;

% psth plots
figure('color','w','position',[200 200 800 400])
subplot(121)
psth1 = mean(Cell1_awake_hist_bin);
sm = conv(psth1, g1);
endp = startp + length(psth1) - 1;
psth1_sm = sm(startp:endp);
plot(timevec,psth1_sm,'linewidth',2)
hold on
psth2 = mean(Cell2_awake_hist_bin);
sm = conv(psth2, g1);
endp = startp + length(psth2) - 1;
psth2_sm = sm(startp:endp);
plot(timevec,psth2_sm,'linewidth',2)
plot([0,0],[0,1],'--','color','cyan','linewidth',2)
xlabel('Time from SWR onset (ms)')
ylabel('Firing rate (Hz)')
title('Awake SWRs')
legend('Cell 1','Cell 2')
box off
subplot(122)
psth1 = mean(Cell1_sleep_hist_bin);
sm = conv(psth1, g1);
endp = startp + length(psth1) - 1;
psth1_sm = sm(startp:endp);
plot(timevec,psth1_sm,'linewidth',2)
hold on
psth2 = mean(Cell2_sleep_hist_bin);
sm = conv(psth2, g1);
endp = startp + length(psth2) - 1;
psth2_sm = sm(startp:endp);
plot(timevec,psth2_sm,'linewidth',2)
plot([0,0],[0,1],'--','color','cyan','linewidth',2)
xlabel('Time from SWR onset (ms)')
ylabel('Firing rate (Hz)')
title('Sleep SWRs')
legend('Cell 1','Cell 2')
box off
%%
%-----------SWR CCG----------%
Cell1hist2 = Cell1_awake_hist';
Cell2hist2 = Cell2_awake_hist';
Cell1h = Cell1hist2(:);
Cell2h = Cell2hist2(:);
[~, rawcorr_swr_awake, corr_sm_swr_awake] = spiketrainxcorr(find(Cell1h)./1000,find(Cell2h)./1000,bin_swr,tmax,sw1_swr);
Cell1hist2=Cell1_sleep_hist';
Cell2hist2=Cell2_sleep_hist';
Cell1h=Cell1hist2(:);
Cell2h=Cell2hist2(:);
[timebase_swr, rawcorr_swr_sleep, corr_sm_swr_sleep] = spiketrainxcorr(find(Cell1h)./1000,find(Cell2h)./1000,bin_swr,tmax,sw1_swr);

figure('color','w')
plot(timebase_swr.*1000,corr_sm_swr_sleep,'k','linewidth',2)
hold on
plot(timebase_swr.*1000,corr_sm_swr_awake,'r','linewidth',2)
plot([0,0],[0,0.05],'k--')
xlabel('Time lag (ms)')
ylabel('Correlation')
title('Crosscorrelogram during SWRs')
legend('Sleep','Awake')
box off
