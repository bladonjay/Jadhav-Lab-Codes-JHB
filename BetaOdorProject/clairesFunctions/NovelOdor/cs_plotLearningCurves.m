%cs_plotLearningCurves

[topDir,figDir] = cs_setPaths();
animals = {'CS39','CS41','CS42','CS44'};
trialtypes = {'novelodor','novelodor2'};
threshold = 0.60;

for a = 1:length(animals)
    animal = animals{a};
    
    for t = 1:length(trialtypes)
        trialtype = trialtypes{t};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        dayepochs = cs_getRunEpochs(animDir, animal, trialtype);
        %add buffer to fill in with number of trials per epoch
        dayepochs(:,3) = NaN;
        
        %concatenate binary perf
        perf = [];
        for ep = 1:size(dayepochs,1)
            daystr = getTwoDigitNumber(dayepochs(ep,1));
            
            load([animDir,'BinaryPerf\',animal,'BinaryPerf',daystr])
            perf = [perf;BinaryPerf{dayepochs(ep,1)}{dayepochs(ep,2)}];
            dayepochs(ep,3) = length(BinaryPerf{dayepochs(ep,1)}{dayepochs(ep,2)});
        end
        
        if ~isempty(perf)
        [pc] = getestprobcorrect_niceplot(perf, 0.5, 2, 1);
        hold on
        %first time when learning probability crosses threshold
        crossing = find(pc(:,2) >= threshold,1,'first');
        
        plot([1 size(pc,1)], [threshold threshold],'r.-.')
        
        figfile = [figDir,'Behavior\',animal,'_',trialtype];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        saveas(gcf,figfile,'fig');
        close all
        end
        
    end
end

