clear
close all
region = 'CA1';

[topDir,figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];

load([dataDir, 'PSDdata_run',region]);
run = PSDdata.sessions_post_norm;
odor = PSDdata.sessions_pre_norm;
fpass = PSDdata.fpass;
stepsize = (fpass(2)-fpass(1))/size(run,2);
freq = fpass(1):stepsize:fpass(2)-stepsize;

new_odor = [];
new_run = [];
newfreq = freq(1):stepsize/5:freq(end);
for t = 1:size(odor,1) 
    new_odor = [new_odor;smoothdata(interp1(freq,odor(t,:),newfreq),'gaussian',15)];
    new_run = [new_run;smoothdata(interp1(freq,run(t,:),newfreq),'gaussian',15)];
end

%remove CS41/CS42, and combine epochs for each day
load([topDir,'AnalysesAcrossAnimals\allEps_odorplace']);
uniquedays = unique(allEps(:,[1,2]),'rows');
rm = (uniquedays(:,1)==5 | uniquedays(:,1)==6 | uniquedays(:,1) == 7);
uniquedays = uniquedays(~rm,:);

tmp_odor = [];
tmp_run = [];
for d = 1:length(uniquedays)
    %find anim day index that match from all epochs
    ind= find(ismember(allEps(:,[1,2]),uniquedays(d,:),'rows'));
    tmp_odor = [tmp_odor;mean(new_odor(ind,:),1)];
    tmp_run = [tmp_run;mean(new_run(ind,:),1)];
end

new_odor = tmp_odor;
new_run = tmp_run;

[~,peakodorind] = max(new_odor(:,43:71),[],2);
[~,peakrunind] = max(new_run(:,43:71),[],2);

fr = newfreq(43:71);
peakodor = fr(peakodorind);
peakrun = fr(peakrunind);

p_CA1odorRun = signrank(peakodor,peakrun);

load([dataDir, 'PSDdata_runTC']);
odor = PSDdata.sessions_pre_norm;
new_odor = [];
for t = 1:size(odor,1)
        new_odor = [new_odor;smoothdata(interp1(freq,odor(t,:),newfreq),'gaussian',15)];
end

%combine epochs for each day
load([topDir,'AnalysesAcrossAnimals\allEps_odorplace']);
uniquedays = unique(allEps(:,[1,2]),'rows');
good = (uniquedays(:,1)==6 | uniquedays(:,1) == 7);
uniquedays = uniquedays(good,:);
gooddays = allEps(allEps(:,1) == 6 | allEps(:,1) == 7,:);

tmp_odor = [];
for d = 1:length(uniquedays)
    %find anim day index that match from all epochs
    ind= find(ismember(gooddays(:,[1,2]),uniquedays(d,:),'rows'));
    tmp_odor = [tmp_odor;mean(new_odor(ind,:),1)];
end
new_odor = tmp_odor;
[~,peakodorind] = max(new_odor(:,36:71),[],2);
peakodor_TC = fr(peakodorind);

cs_boxplot([peakodor_TC';peakodor';peakrun'],[zeros(length(peakodor_TC),1);ones(length(peakodor),1);ones(length(peakrun),1)+1])

ylim([5 10])
ylabel('Peak Frequency (Hz)')
xticklabels({'TCodor','CA1odor','CA1run'})

p_CA1odorTCodor = ranksum(peakodor_TC,peakodor);
 

nl = newline;
text(1,5,['p CA1odorRun= ',num2str(p_CA1odorRun),nl,'p CA1odorTCodor= ',num2str(p_CA1odorTCodor)]);


figfile = [figDir,'EEG\powerOdorVrun_boxplot'];

print('-dpdf', figfile);
print('-djpeg',figfile);

