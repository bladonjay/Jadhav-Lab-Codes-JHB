clear
animals = {'CS39','CS41','CS42','CS44'};
trialtypes = {'novelodor','novelodor2','noodor'};

%Perform state space analysis on novel odor sessions to determine when the
%animal learned the new odor pairs.

%Then add fields to odorTriggers structure to separate trials into
%pre-learning and post-learnring.

[topDir,figDir] = cs_setPaths();
threshold = 0.60;

for a = 1%1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    for tr = 1:length(trialtypes)
        trialtype = trialtypes{tr};
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
        
        if isempty(perf)
            continue
        end
        [pc] = getestprobcorrect_niceplot(perf, 0.5, 0, 1);
        hold on
        plot([1 length(pc)], [threshold threshold], 'r:');
        figfile = [figDir, 'Behavior\',animal,'_',trialtype];
        print('-dpdf',figfile);
        print('-djpeg',figfile);
        close all
        
        %first time when learning probability crosses threshold
        crossing = find(pc(:,2) >= threshold,1,'first');
        
        %add fields to odorTriggers structures
        tottrials = 0;
        for ep = 1:size(dayepochs,1)
            daystr = getTwoDigitNumber(dayepochs(ep,1));
            
            load([animDir,animal,'odorTriggers',daystr])
            tottrials = tottrials + dayepochs(ep,3);
            if  tottrials < crossing %rat does not learn yet
                prelearn = odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.allTriggers;
                odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.prelearn = prelearn;
                odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.postlearn = [];
            else
                %find ind corresponding to crossing
                ind = dayepochs(ep,3) - (tottrials- crossing);
                if ind > 0 % learning happens in this epoch
                    prelearn = odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.allTriggers(1:ind);
                    postlearn = odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.allTriggers(ind+1:end);
                    odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.prelearn = prelearn;
                    odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.postlearn = postlearn;
                else % learning already occured
                    postlearn = odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.allTriggers;
                    odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.prelearn = [];
                    odorTriggers{dayepochs(ep,1)}{dayepochs(ep,2)}.postlearn = postlearn;
                end
            end
            save([animDir,animal,'odorTriggers',daystr],'odorTriggers');
        end
    end
end