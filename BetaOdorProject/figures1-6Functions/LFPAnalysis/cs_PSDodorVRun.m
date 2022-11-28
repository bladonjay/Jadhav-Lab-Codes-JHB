clear
close all
region = 'CA1';

[topDir,figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];

load([dataDir, 'PSDdata_run',region]);
run = PSDdata.sessions_post;
odor = PSDdata.sessions_pre;
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


psd_odor = mean(new_odor,1);
psd_run = mean(new_run,1);

%only use freq band around RR (5-10 Hz)
newfreq = newfreq(36:71);
psd_odor = psd_odor(36:71);
psd_run = psd_run(36:71);

norm_odor = rescale(psd_odor);
norm_run = rescale(psd_run);

figure,
hold on
%plot(freq,psd_odor)
plot(newfreq,psd_odor);
%plot(freq,psd_run)
plot(newfreq,psd_run);

[maxodor,ind] = max(psd_odor);
plot([newfreq(ind) newfreq(ind)],[0 maxodor],'k--');

[maxrun,ind] = max(psd_run);
plot([newfreq(ind) newfreq(ind)],[0 maxrun],'k--');


[~,peakodorind] = max(new_odor(:,36:71),[],2);
[~,peakrunind] = max(new_run(:,36:71),[],2);

peakodor = newfreq(peakodorind);
peakrun = newfreq(peakrunind);

p = signrank(peakodor,peakrun);
text(5.5,200,['p = ',num2str(p)]);
ylabel('Power')
xlabel('Frequency (Hz)')


% regions = {'OB'};
% for r = 1:length(regions)
%     region = regions{r};
%     load([dataDir, 'PSDdata_run',region]);
%     odor = PSDdata.sessions_pre;
%     new_odor = [];
%     
%     newfreq = freq(1):stepsize/5:freq(end);
%     for t = 1:size(odor,1)
%         new_odor = [new_odor;smoothdata(interp1(freq,odor(t,:),newfreq),'gaussian',15)];
%     end
%     
%     psd_odor = mean(new_odor,1);
%     
%     newfreq = newfreq(36:71);
%     psd_odor = psd_odor(36:71);
%     
%     norm_odor = rescale(psd_odor);
%     
%     plot(newfreq,norm_odor);
%     
% end
figfile = [figDir,'IgarashiReplication\powerOdorVrun_',region];

print('-dpdf', figfile);
print('-djpeg',figfile);

