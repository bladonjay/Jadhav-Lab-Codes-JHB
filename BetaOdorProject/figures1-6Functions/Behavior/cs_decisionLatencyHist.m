%cs_decisionLatencyHist

%finds nosepoke duration across all animals/trials and compares between
%correct and incorrect trials. 

[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

duration = []; 
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    
    for day = days'
        epochs = runeps(runeps(:,1)==day,2);
        
        for ep = epochs'
            
            %correct trials
            
            nps = nosepokeWindow{day}{ep};
            dur = nps(:,2)-nps(:,1);
            
            duration = [duration; dur];
           
        end
    end
end
figure
[x, b] = histcounts(duration,[0.5:0.05:1.5,inf]);
histogram(duration,b)
box off
xlim([0.5 1.55])
xlabel('Nosepoke Duration (s)');
ylabel('Number of Trials');
mn = mean(duration);
sd = std(duration);

figfile = [figDir,'Behavior\NosepokeDurationHist'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    
