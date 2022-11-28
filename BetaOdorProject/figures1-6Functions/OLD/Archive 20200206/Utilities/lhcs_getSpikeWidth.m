function spikewidth = lhcs_getSpikeWidth(topDir, animal, day, tet, cellnum)
animdir = [topDir, animal,'Expt\',animal,'\'];
daystr = getTwoDigitNumber(day);

dayfolder = dir([animdir,daystr,'*']);
dayfolder = dayfolder.name;
datestr = erase(dayfolder,[daystr,'_']);

load([animdir,dayfolder,'\',animal,'_',datestr,'.matclust\','matclust_param_nt',num2str(tet),'.mat'])% load matclustparams and waves for the tetrode of interest
load([animdir,dayfolder,'\',animal,'_',datestr,'.matclust\','waves_nt',num2str(tet),'.mat'])% load matclustparams and waves for the tetrode of interest

wavesIdx= clustattrib.clusters{cellnum}.index; % get idx of spikes in clust
curWaves= waves(:,:,wavesIdx);  % get waveforms of spike
curWavesCh1=curWaves;
curWavesCh1(:,2)=[]; %delete middle field to get amplitudes for channel 1
[peak,peakIdx]=(max(curWavesCh1(:,:))); 
    for w=1:length(peakIdx)
        [trough(w), trouIdx(w)]=min(curWavesCh1(peakIdx(1,w):end,w));
        tempid = find(curWavesCh1(:,w)==trough(w));
        trouIdx(w)= tempid(end);
    end
spikewidth= mean(double(trouIdx-peakIdx)*(1/30000)*1000); % in ms
