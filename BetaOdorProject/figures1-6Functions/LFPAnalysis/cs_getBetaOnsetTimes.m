%beta onset times?
%separate beta into trials

topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC','OB'};

bufferwin = [0.2 0.2];
%extends window around NP, since beta can start early. Make sure this is
%the same as in cs_getHighBetaTimes

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
    load([animDir,animal,'highBeta.mat']);
    
    days = 1:length(highBeta);
    
    for d = 1:length(days)
        day = days(d);
        daystr = getTwoDigitNumber(day);
        
        load([animDir, animal, 'nosepokeWindow',daystr,'.mat']);
        
        epochs = find(~cellfun(@isempty, nosepokeWindow{day}));
        
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            
            windows = nosepokeWindow{day}{epoch};
            beta = highBeta{day}{epoch};
            epochwinbeta = zeros(size(windows,1),2);
            for w = 1:size(windows,1)
                window = [windows(w,1)- bufferwin(1)  windows(w,2)+bufferwin(2)];
                winbeta = [nan, nan];
                if exist('trigbeta','var')
                    clear trigbeta
                end
                for r = 1:length(regions)
                    region = regions{r};
                    regionbeta = beta.(region);
                    if ~isempty(regionbeta)
                        %                          regionbeta = regionbeta(:,1);
                        trigbetainds = find(regionbeta(:,1) >= window(1) & regionbeta(:,1) < window(2));
                        if ~isempty(trigbetainds)
                            trigbeta.(region) = [regionbeta(trigbetainds(1),1), regionbeta(trigbetainds(end),2)];
                        else
                            clear trigbeta
                            break
                        end
                    else
                        clear trigbeta
                        break
                    end
                    
                end
                if exist('trigbeta','var')
                    alltrigbeta = cell2mat(struct2cell(trigbeta));
                    winbeta = [max(alltrigbeta(:,1)), min(alltrigbeta(:,2))];
                    %check to make sure this works
                    
                end
                epochwinbeta(w,(1:2)) = winbeta;
            end
            betaWindows{day}{epoch} = epochwinbeta; %#ok<SAGROW>
        end
    end
    save([animDir, animal,'betaWindows.mat'],'betaWindows');
    clear betaWindows
end

