%cs_npWinCoherence

%load coherence files and np windows
%separate epoch-long coherence matirices into average coh on each trial at
%each freq:
%row = freq band, column = trial

%append into existing coherence files

clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31'};
%animals = {'CS33'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

for a = 1:length(animals)
    animal = animals{a};
     animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
     npWindow= loaddatastruct(animDir,animal,'nosepokeWindow');
     for r = 1:length(regions)
         region = regions{r};
     cohFiles = dir([animDir, animal,'coherence',region,'0*']);
     for f = 1:length(cohFiles)
         load(cohFiles(f).name);
         day = find(~cellfun(@isempty, coherence));
         
         eps = find(~cellfun(@isempty, coherence{day}));
         
         for ep = eps
             if isempty(coherence{day}{ep}.Coh)
                 continue
             end
             coh = coherence{day}{ep}.Coh;
             times = coherence{day}{ep}.time;
             
             windows = npWindow{day}{ep};
             
             npWinCoh = [];
             for w = 1:size(windows,1)
                 win = windows(w,:);
                 trialind = isExcluded(times,windows);
                 npWinCoh = [npWinCoh, mean(coh(:,trialind),2)];
             end
             
             coherence{day}{ep}.npWinCoh = npWinCoh;
         end
         save([animDir,cohFiles(f).name],'coherence');
         clear coherence
     end
     
     end
end