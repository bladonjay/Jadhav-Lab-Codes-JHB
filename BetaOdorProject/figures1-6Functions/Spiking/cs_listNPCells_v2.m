% find all cells that are odor-responsive
% uses 2 critera:
% - must spike at least once per trial on average (i.e. total # spikes >= #
%   trials
% - must have a significant change in FR from pre-odor period (using sign rank
%   test) for at least one odor


clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
    cellct=0; activect=0;
    region = regions{r};
    
    npCells =[]; activeCells=[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        % get all cells for animal
        load([animDir,animal,'cellinfo.mat'])
        
        % this is for pyramidal cells in the area with more than 100 spikes
        cellfilter = ['isequal($area,''',region,''') && isequal($type, ''pyr'') && $numspikes > 100'];
        cells = evaluatefilter(cellinfo,cellfilter);
        
        % remove epochs
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        runEps = cs_getRunEpochs(animDir,animal,'odorplace');
        days = unique(runEps(:,1));
        
        spikes = loaddatastruct(animDir, animal, 'spikes',days);
        nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
        odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
       
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            
            eps = runEps((runEps(:,1) == day),2);
            
            %loop through cells on each day
            for c = 1:size(daycells,1)
                
                npspikes = [];
                prespikes = [];
                totaltrigs = 0;
                cell = daycells(c,:);
                triallabels = [];

                 
                for ep = 1:length(eps)
                    epoch = eps(ep);
                    
                    % find spikes during odor and time-matched pre-odor
                    % period
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                            
                            trigs = nosepokeWindow{day}{epoch};
                            totaltrigs = totaltrigs + size(trigs,1);
                            
                            winlength = nosepokeWindow{day}{epoch}(:,2)- nosepokeWindow{day}{epoch}(:,1);
                            pretrigs = [nosepokeWindow{day}{epoch}(:,1)-winlength, nosepokeWindow{day}{epoch}(:,1)];
                            
                            
                            for t = 1:size(trigs,1)
                                
                                %gather spikes
                                winspikes = epspikes(logical(isExcluded(epspikes, trigs(t,:))));
                                npspikes = [npspikes; length(winspikes)];
                                
                                pre = epspikes(logical(isExcluded(epspikes,pretrigs(t,:))));
                                prespikes = [prespikes;length(pre)];
                            end
                            
                            %gather trial labels
                            labels = zeros(size(trigs,1),2);
                            [cr,cl] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            labels([cl;cr],1) = 1;
                            labels(cl,2) = 1;
                            triallabels = [triallabels; labels];
                        end
                        
                    end
                    
                end
                cellct=cellct+1;
                
                if ~isempty(npspikes) || ~isempty(prespikes)
                    correctleft = find(triallabels(:,1) == 1 & triallabels(:,2) == 1);
                    correctright = find(triallabels(:,1) == 1 & triallabels(:,2) == 0);
                    
                    lspikes = [npspikes(correctleft), prespikes(correctleft)];
                    rspikes = [npspikes(correctright), prespikes(correctright)];
                    allspikes = [lspikes;rspikes];
                    %do signrank tests
                    p1 = signrank(lspikes(:,1),lspikes(:,2));
                    p2 = signrank(rspikes(:,1),rspikes(:,2));
                    p3 = signrank(allspikes(:,1),allspikes(:,2)); % or both combined...
                    
                    if sum(npspikes)/totaltrigs >=1 %at least 1 spike per trial during NP
                        activeCells=[activeCells; a, cell];
                        if p1 < 0.05 || p2 <0.05 || p3 < 0.05%at least one odor response must be significant
                            npCells = [npCells; a, cell]; %if it meets criteria, add the cell to the list
                        end
                    end
                end
                
               
            end
        end
    end
    fprintf('%d total pyrams in %s\n', cellct,regions{r})
    fprintf('%d active pyrams in %s\n',length(activeCells),regions{r});
    fprintf('%d np pyrams in %s\n',length(npCells),regions{r});
    save([dataDir,'npCells_',region,'.mat'], 'npCells','activeCells');
end
%clear