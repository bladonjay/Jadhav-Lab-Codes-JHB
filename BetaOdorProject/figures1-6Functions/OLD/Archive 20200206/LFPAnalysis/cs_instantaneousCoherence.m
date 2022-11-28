function cs_instantaneousCoherence(animals, topDir, regions, figDir)

%for each pair of regions, calculate instantaneous coherence 200ms after
%trigger- across all epochs and animals
figure, hold on
for r = 1:length(regions)
    region = regions{r};
    
    allData = [];
    
    for a = 1:length(animals)
        animal = animals{a};
        
        coherencefiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'coherence',region,'*']);
        odorfiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'odorTriggers*']);
        
        for d = 1:length(coherencefiles)
            day = d; 
            load([topDir,animal,'Expt\',animal,'_direct\',coherencefiles(d).name])
            load([topDir,animal,'Expt\',animal,'_direct\',odorfiles(d).name])
            
            coherence = coherence{1,d};
%             pos = pos{1,d};
            odorTriggers = odorTriggers{1,d};
            epochs = find(~cellfun(@isempty, coherence));
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                cohtimes = coherence{1,epoch}.time;
                freqs = coherence{1,epoch}.freqs;
                coherence_tmp = coherence{1,epoch}.coherence;
                trigs = odorTriggers{1,epoch}.allTriggers;
                trigs = trigs+0.2; %200ms after trigger
                trigsinds = find(cohtimes == trigs);
               
                %find index of nearest coherence times to trigger times
                TMP = bsxfun(@(x,y) abs(x-y), trigs(:), reshape(cohtimes,1,[]));
                [~, indx] = min(TMP,[],2) ;
                
                newCoh = coherence_tmp(:,indx);
                
                allData = [allData, newCoh];
            end
        end
                
    
    end
    meanAllAnimals(r,:) = mean(allData,2)';
    semAllAnimals(r,:) = std(allData,0,2)/sqrt(size(allData,2));
end

plot(freqs, meanAllAnimals(1,:),'-b');
plot(freqs, meanAllAnimals(2,:),'-r');
plot(freqs, meanAllAnimals(3,:),'-g');

uE = meanAllAnimals + semAllAnimals;
lE = meanAllAnimals - semAllAnimals;
yP = [lE,fliplr(uE)];
xP=[freqs,fliplr(freqs)];

patch(xP,yP(1,:),1, 'facecolor', 'b', 'edgecolor','none', 'facealpha',0.2)
patch(xP,yP(2,:),1, 'facecolor', 'r', 'edgecolor','none', 'facealpha',0.2)
patch(xP,yP(3,:),1, 'facecolor', 'g', 'edgecolor','none', 'facealpha',0.2)


% shadedErrorBar(freqs,meanAllAnimals(1,:),semAllAnimals(1,:), 'lineprops','-b');
% shadedErrorBar(freqs,meanAllAnimals(2,:),semAllAnimals(2,:), 'lineprops','-r');
% shadedErrorBar(freqs,meanAllAnimals(3,:),semAllAnimals(3,:), 'lineprops','-g');

xlabel('Frequency (Hz)')
ylabel('Coherence')
legend(regions)


figfile = [figDir,'1i_MeanCohAsFuncOfFreq'];

print('-djpeg', figfile);
print('-dpng', figfile);
saveas(gcf,figfile,'fig');

