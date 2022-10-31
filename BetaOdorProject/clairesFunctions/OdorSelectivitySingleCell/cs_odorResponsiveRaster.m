clear
[topDir, figDir] = cs_setPaths();

win = [0.5 1];

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};

load([figDir, 'cmap_selectivity.mat']);

%%

for r = 1:length(regions)
    region = regions{r};
    %find cells that meet selection criteria (pyramidal cells in area of
    %interest)
    load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    data = {};
    inds = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        cells = npCells((npCells(:,1) == a),2:end);
         load([animDir,animal,'cellinfo.mat'])
%         
%         %Take previously defined selective cells
%         cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region, ...
%             '''))'];
%         
%         cells = evaluatefilter(cellinfo,cellfilter);
%         
%         noeps = cells(:,[1 3 4]);
%         cells = unique(noeps,'rows');
%         
        days = unique(cells(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            %runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            
            for c = 1:size(daycells,1)
                
                correctleftspikes = [];
                correctrightspikes = [];
                incorrectleftspikes = [];
                incorrectrightspikes = [];
                
                cell = daycells(c,:);
                runeps = cs_getRunEpochs(animDir,animal,'odorplace',day);
                runeps = runeps(:,2);
                goodeps = cs_findGoodEpochs(cellinfo{day}, {'meanrate'},cell(2:3));
                epochs = intersect(runeps,goodeps);
                if all(cell == [3 1 7])
                    keyboard
                end
                for epoch = epochs'
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        
                        [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                        %triginds = sort([correct_left; correct_right]);
                        
                        trigs = odorTriggers{day}{epoch}.allTriggers;
                        %% Find Spikes on each trial
                        for t = 1:length(trigs)
                            trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                            winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                            
                            winspikes = winspikes -trigs(t);
                            
                            %classify into trial type
                            if ismember(t, correct_left)
                                correctleftspikes{end+1,1} = winspikes;
                            elseif ismember(t, correct_right)
                                correctrightspikes{end+1,1} = winspikes;
                            elseif ismember(t, incorrect_left)
                                incorrectleftspikes{end+1,1} = winspikes;
                            elseif ismember(t, incorrect_right)
                                incorrectrightspikes{end+1,1} = winspikes;
                            end
                            
                        end
                        
                    end
                    
                end
                
                data{end+1,1}.cleft = correctleftspikes;
                data{end,1}.cright = correctrightspikes;
                data{end,1}.ileft = incorrectleftspikes;
                data{end,1}.iright = incorrectrightspikes;
                
                inds = [inds; a, cell];
                
                
            end
            
        end
        
        
    end
    
    rasterData.data = data;
    rasterData.inds = inds;
    save([topDir,'AnalysesAcrossAnimals\rasterData_np_',region,],'rasterData');
    
end
