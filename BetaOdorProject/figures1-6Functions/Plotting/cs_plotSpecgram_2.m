%get coherogram using pre-calculated coherence files. simply load coherence
%files, take times within trial windows, and average. Should already be a
%zscore.
clear
[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS41','CS42'};
%regions  = {'CA1','PFC','OB'};
regions = {'TC'};
trialtypes = {'odorplace'};
freqband = 'floor'; maxfreq = 15;
%freqband = 'low'; maxfreq = 40;

timewin = [0.5 1.5];

for r = 1:length(regions)
    region = regions{r};
Spec_allanimals = [];
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
    
    Spec = [];
    for d = days'
        daystr = getTwoDigitNumber(d);
        eps = runeps(runeps(:,1) == d,2);
        load([animDir, animal, 'spec', freqband, region, daystr]);
        odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',d);
        for ep = eps'
            specfull = spec{d}{ep}.Spec';
            eptime = spec{d}{ep}.time;
            trigs = odorTriggers{d}{ep}.allTriggers;
            wins = [trigs - timewin(1), trigs + timewin(2)];
            
            for w = 1:size(wins,1)
                win = wins(w,:);
                inds = isExcluded(eptime, win);
                if win(2) <= eptime(end)
                    if sum(inds) >0
                        s = specfull(:,inds);
                        
                        try
                            Spec = cat(3,Spec,s);
                        catch
                            s = s(:,1:size(Spec,2));
                            Spec = cat(3,Spec,s);
                        end
                    end
                end
            end
        end
    end
    Spec_allanimals = cat(3,Spec_allanimals,mean(Spec,3));
end
freq = spec{d}{ep}.freq;
Spec = mean(Spec_allanimals,3);

Spec = interp2(Spec,3);
freq = 1:(maxfreq-1)/size(Spec,1):maxfreq-((maxfreq-1)/size(Spec,1));
times = -timewin(1):(timewin(2)+timewin(1))/size(Spec,2):timewin(2)-(timewin(2)+timewin(1))/size(Spec,2);



figure,
imagesc(times,freq,Spec)
colormap(hot);
set(gca,'YDir','normal')

hold on
plot([0 0],[1 max(freq)], 'k--')
xlabel('Time from odor onset (seconds)');
ylabel('Frequency (Hz)');
colorbar

figfile = [figDir,'Specgrams\',region,'_specgram_',freqband];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end