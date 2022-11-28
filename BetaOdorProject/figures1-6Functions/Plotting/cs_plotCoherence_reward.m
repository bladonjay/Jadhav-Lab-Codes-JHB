%get coherogram using pre-calculated coherence files. simply load coherence
%files, take times within trial windows, and average. Should already be a
%zscore.

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

regions = {'CA1-PFC','CA1-OB','PFC-OB'};

trialtypes = {'odorplace'};

%freqband = 'floor'; maxfreq = 15;
freqband = 'low'; maxfreq = 40;

timewin = [1.5 0.5];

for r = 1:length(regions)
    region = regions{r};
Coh_allanimals = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    runeps = [];
    for t = 1:length(trialtypes)
        eps = cs_getRunEpochs(animDir, animal, trialtypes{t});
        runeps = [runeps;eps];
    end
    
    runeps = sortrows(runeps);
    
    days = unique(runeps(:,1));
    
    Coh = [];
    for d = days'
        daystr = getTwoDigitNumber(d);
        eps = runeps(runeps(:,1) == d,2);
        load([animDir, animal, 'coherence', region, daystr]);
        triggers = loaddatastruct(animDir, animal, 'rewards',d);
        for ep = eps'
            cohfull = coherence{d}{ep}.Coh;
            eptime = coherence{d}{ep}.time;
            trigs = sort([triggers{d}{ep}.leftWindows(:,2); triggers{d}{ep}.rightWindows(:,2)]);
            wins = [trigs - timewin(1), trigs + timewin(2)];
            
            for w = 1:size(wins,1)
                win = wins(w,:);
                inds = isExcluded(eptime, win);
                if win(2) <= eptime(end)
                    if sum(inds) >0
                        coh = cohfull(:,inds);
                        
                        try
                            Coh = cat(3,Coh,coh);
                        catch
                            coh = coh(:,1:size(Coh,2));
                            Coh = cat(3,Coh,coh);
                        end
                    end
                end
            end
        end
    end
    Coh_allanimals = cat(3,Coh_allanimals,mean(Coh,3));
end
freq = coherence{d}{ep}.freq;
ind = freq <= maxfreq;
freq = freq(ind);
Coh_allanimals = Coh_allanimals(ind,:,:);
Coh = mean(Coh_allanimals,3);

Coh = interp2(Coh,3);
freq = 1:(maxfreq-1)/size(Coh,1):maxfreq-((maxfreq-1)/size(Coh,1));
times = -timewin(1):(timewin(2)+timewin(1))/size(Coh,2):timewin(2)-(timewin(2)+timewin(1))/size(Coh,2);



figure,
imagesc(times,freq,Coh)
colormap(jet);
set(gca,'YDir','normal')

hold on
plot([0 0],[1 max(freq)], 'k--')
xlabel('Time from odor onset (seconds)');
ylabel('Frequency (Hz)');
colorbar

figfile = [figDir,'Cohgrams\',region,' coherence_',freqband,'_reward'];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end