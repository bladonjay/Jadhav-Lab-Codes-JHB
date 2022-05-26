%cs_preferredOdorPhaseLocking

%compare phase locking (kappa?) of selective cells on preferred vs
%non-preferred odor trials
clear
close all
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35'};

cellregions = {'CA1', 'PFC'};
eegregions = {'CA1','PFC','OB'};

for er = 1:length(eegregions)
    eegregion = eegregions{er};
    
    for r = 1:length(cellregions)
        cellregion = cellregions{r};
        
        load([topDir, 'AnalysesAcrossAnimals\selectiveCells_',cellregion,'.mat']);
        aninds = unique(selectivecells(:,1));
        allcells = [];
        for a = 1:length(aninds)
            animal = animals{aninds(a)};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            animcells = selectivecells(selectivecells(:,1) == aninds(a),2:4);
            days = unique(animcells(:,1));
            
            load([animDir, animal,'cellinfo.mat']);
            load([animDir, animal,'betaWindows.mat']);
            load([animDir, animal,'tetinfo.mat']);
            tetfilter = ['((strcmp($area, ''',eegregion,''')) && (strcmp($descrip2, ''betatet'')))'];
            
            
            epochmat = cs_getRunEpochs(animDir, animal, 'odorplace');
            
            for d = 1:length(days)
                day = days(d);
                daystr = getTwoDigitNumber(day);
                epochs = epochmat(epochmat(:,1) == day, 2);
                daycells = animcells(animcells(:,1) == day,2:3);
                
                load([animDir, animal,'odorTriggers',daystr,'.mat']);
                load([animDir, animal,'spikes',daystr,'.mat']);
                
                for c = 1:size(daycells,1)
                    cell = daycells(c,:);
                   
                    %Determine odor preference of cell
                    odorPref = cellinfo{day}{epochs(1)}{cell(1)}{cell(2)}.selectivity;
                    day_kappa = [];
                    
                    celleps = cs_getCellEpochs(cellinfo,day, cell(1),cell(2),epochs);
                    
                    
                    for ep = 1:length(celleps)
                        epoch = celleps(ep);
                        epstr = getTwoDigitNumber(epoch);
                        %Get spiking on preferred trials, and on opposite trials
                        betawins = betaWindows{day}{epoch};
                        testbetawins = find(~isnan(betawins));

                        if ~isempty(testbetawins)
                        
                             %Get high beta periods for each of those trial groups
                            [cl, cr, ~,~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            leftwins = betawins(cl,:);
                            rightwins = betawins(cr,:);

                            stime = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                            
                            betatet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);
                            tetstr = getTwoDigitNumber(betatet);
      
                            load([animDir, 'EEG/', animal, 'beta',daystr,'-',epstr,'-',tetstr,'.mat'],'beta');
                            t = geteegtimes(beta{day}{epoch}{betatet});
                            betaphase = beta{day}{epoch}{betatet}.data(:,2);
                            sph = betaphase(lookup(stime, t));
                            
                             %calc phase locking - kappa value
                            leftwins = leftwins(~any(isnan(leftwins),2),:);
                            if ~isempty(leftwins)
                                    goodspikes = isExcluded(stime, leftwins);
                                    newsph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                                    
                                    if length(newsph)> 100
                                        %[m, ph] = modulation(sph);
                                         [prayl, zrayl] = circ_rtest(newsph);
                                        kappa_L = prayl;
%                                          [~, kappa_L] = circ_vmpar(newsph);
                                    end
                            end
                            
                             %calc phase locking - kappa value
                             rightwins = rightwins(~any(isnan(rightwins),2),:);
                            if ~isempty(rightwins)
                                    goodspikes = isExcluded(stime, rightwins);
                                    newsph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                                    
                                    if length(newsph)> 100
                                        %[m, ph] = modulation(sph);
                                         [prayl, zrayl] = circ_rtest(newsph);
                                         kappa_R = prayl;
                                         %[~, kappa_R] = circ_vmpar(newsph);
                                         
                                    end
                            end
                            
                            switch odorPref
                                case 'leftSelective'
                                    epoch_kappa = [kappa_L, kappa_R];
                                case 'rightSelective'
                                    epoch_kappa = [kappa_R, kappa_L];
                            end
                            
                        end
                           day_kappa = [day_kappa;epoch_kappa]; 
                    end
                    mn = mean(day_kappa,1);
                    allcells = [allcells;mn];
                   
                    %save matrix of values for correct/incorrect
                    
                end
            end

        end
        %plot change in kappa - 
        toPlot{r} = allcells;
    end
    figure, hold on
    for r = 1:length(cellregions)
        plot([r, r+0.5], toPlot{r})
        [~,p] = ttest(toPlot{r}(:,1),toPlot{r}(:,2))
    end
    
end