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


psd_odor = mean(new_odor,1);
psd_run = mean(new_run,1);
sem_odor = stderr(new_odor);
sem_run = stderr(new_run);



%only use freq band around RR (5-10 Hz)
newfreq = newfreq(36:71);
psd_odor = psd_odor(36:71);
psd_run = psd_run(36:71);
sem_odor = sem_odor(36:71);
sem_run = sem_run(36:71);

norm_odor = rescale(psd_odor);
norm_run = rescale(psd_run);

figure,
hold on
%plot(freq,psd_odor)
plot(newfreq,norm_odor,'k-');
%plot(freq,psd_run)
plot(newfreq,norm_run,'k-');

[maxodor,ind] = max(norm_odor);
plot([newfreq(ind) newfreq(ind)],[0 maxodor],'k--');
patch([newfreq fliplr(newfreq)], [norm_odor+sem_odor, fliplr(norm_odor-sem_odor)],'k','FaceAlpha',0.25)


[maxrun,ind] = max(norm_run);
plot([newfreq(ind) newfreq(ind)],[0 maxrun],'k--');
patch([newfreq fliplr(newfreq)], [norm_run+sem_run, fliplr(norm_run-sem_run)],'k','FaceAlpha',0.25)


[~,peakodorind] = max(new_odor(:,36:71),[],2);
[~,peakrunind] = max(new_run(:,36:71),[],2);

peakodor = newfreq(peakodorind);
peakrun = newfreq(peakrunind);

peakCA1 = peakodor;

p_CA1odorRun = signrank(peakodor,peakrun);
%text(5.5,0.1,['p = ',num2str(p_CA1odorRun)]);
ylabel('Normalized Power')
xlabel('Frequency (Hz)')

colors = {'r','g','m'};
regions = {'TC'};
for r = 1:length(regions)
    region = regions{r};
    load([dataDir, 'PSDdata_run',region]);
    odor = PSDdata.sessions_pre_norm;
    run = PSDdata.sessions_post_norm;
    new_odor = [];
    new_run = [];
    
    newfreq = freq(1):stepsize/5:freq(end);
    for t = 1:size(odor,1)
        new_odor = [new_odor;smoothdata(interp1(freq,odor(t,:),newfreq),'gaussian',15)];
        new_run = [new_run;smoothdata(interp1(freq,run(t,:),newfreq),'gaussian',15)];
    end
    
    psd_odor = mean(new_odor,1);
    sem_odor = stderr(new_odor);
    psd_run = mean(new_run,1);
    sem_run = stderr(new_run);
    
    newfreq = newfreq(36:71);
    psd_odor = psd_odor(36:71);
    psd_run = psd_run(36:71);
    sem_odor = sem_odor(36:71);
    sem_run = sem_run(36:71);

    
    norm_odor = rescale(psd_odor);
    norm_run = rescale(psd_run);
%     sem_odor = rescale(sem_odor);
%     sem_run = rescale(sem_run);
    
   
    plot(newfreq,norm_odor,[colors{r},'-']);
    %plot(newfreq,norm_run,[colors{r},'-']);
    [maxodor,ind] = max(norm_odor);
    plot([newfreq(ind) newfreq(ind)],[0 maxodor],[colors{r},'--']);
    hold on
    patch([newfreq fliplr(newfreq)], [norm_odor+sem_odor, fliplr(norm_odor-sem_odor)],'r','FaceAlpha',0.25)

    [~,peakodorind] = max(new_odor(:,36:71),[],2);
    %[~,peakrunind] = max(new_run(:,36:71),[],2);

    peakodor = newfreq(peakodorind);
    %peakrun = newfreq(peakrunind);
    
    p_TCodorCA1un = ranksum(peakodor,peakrun);
    
    p_CA1odorTCodor = ranksum(peakodor,peakCA1);

    
end
nl = newline;
text(5.5,0.1,['p CA1odorRun= ',num2str(p_CA1odorRun),nl,'p TCodorCA1un= ',num2str(p_TCodorCA1un),nl,'p CA1odorTCodor= ',num2str(p_CA1odorTCodor)]);


figfile = [figDir,'EEG\powerOdorVrun'];

print('-dpdf', figfile);
print('-djpeg',figfile);

