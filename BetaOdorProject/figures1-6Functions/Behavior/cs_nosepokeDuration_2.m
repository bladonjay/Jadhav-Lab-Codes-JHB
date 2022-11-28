%cs_nosepokeDuration_2

%finds nosepoke duration across all animals/trials and compares between
%correct and incorrect trials. 

[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

duration_c = []; duration_i = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    
    
    for day = days'
        epochs = runeps(runeps(:,1)==day,2);
        
        dur_c = [];
        dur_i = [];
        for ep = epochs'
            
            %correct trials
            [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
            nps = nosepokeWindow{day}{ep}([cl;cr],:);
            duration = nps(:,2)-nps(:,1);
            
            dur_c = [dur_c; duration];
            
            %incorrect trials
            [~,~,il,ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
            nps = nosepokeWindow{day}{ep}([il;ir],:);
            duration = nps(:,2)-nps(:,1);
            
            dur_i = [dur_i; duration];
            
        end
        duration_c = [duration_c; mean(dur_c)];
        duration_i = [duration_i; mean(dur_i)];
    end
    
end

save([topDir,'AnalysesAcrossAnimals\npDurations'],'duration_c','duration_i');


groups = [ones(length(duration_c),1); zeros(length(duration_i),1)];

cs_boxplot([duration_c;duration_i],groups);
% figure,
% boxplot([duration_c;duration_i],groups);
%axis([0 2.75 0 1.4])
%yticks([0:0.2:1.4])
ylabel('Nosepoke Duration (s)');
xticklabels({'Correct','Incorrect'});



[p] = ranksum(duration_c,duration_i);
text(1.5,0.2,['p = ',num2str(round(p,2))]);

figfile = [figDir,'Behavior\NosepokeDuration'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
