function [meanNPDuration]=cs_meanNosepokeDuration(animals)

%This function finds the average length of nosepokes for each animal on each day. Uses
%odorTriggers file for NP start, then finds the following NP downstate for
%NP off. 

%topDir should be TOP directory containing all animal expt folders
[topDir, figDir] = cs_setPaths();
meanNPDuration = zeros(length(animals),1);
allnps = [];
for a = 1:length(animals)
    animal = animals{a} ;
    
    dataDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    cd(dataDir)
    
    npWindowfiles = dir([animal,'nosepokeWindow*']);
    
    gather = [];
    
    for f = 1:length(npWindowfiles)
        
        file = npWindowfiles(f).name;
        load(file)
        day = find(~cellfun(@isempty,nosepokeWindow));
        epochs = find(~cellfun(@isempty,nosepokeWindow{day})); 
        
        for e = 1:length(epochs)
            epoch = epochs(e);
            durations = nosepokeWindow{1,day}{1,epoch}(:,2) - nosepokeWindow{1,day}{1,epoch}(:,1);
            
            gather = [gather; durations];
            
        end
    end
    allnps = [allnps; gather];
    mn = mean(gather);
    
    meanNPDuration(a,1) = mn;
end

figure, 
histogram(allnps);
meanall = mean(meanNPDuration);
hold on
plot([meanall meanall],[0 200], 'r--');
xlabel('Nosepoke Duration (seconds)')
ylabel('Count')
figtitle = 'Nosepoke Duration Distribution';
figfile = [figDir,'Behavior\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
%disp(meanNPDuration)
