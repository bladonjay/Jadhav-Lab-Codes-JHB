%cs_dynamicPhase
%for each selective neuron, plot its preferred phase of beta in 100ms bins
%during NP (0-1.5 s)
clear
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42'};
regions = {'CA1','PFC'};
betaregions = {'CA1','PFC','OB'};

binsize = 0.1; %100ms bins
win = 1; %1 second after np
%get cells
%load spikes

for br = 1:length(betaregions)
    betaregion = betaregions{br};
    for r = 1:length(regions)
        region = regions{r};
        
        %allphases = [];
        
        %figure, hold on
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
            tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
            
            %find cells
            filt = ['strcmp($area,''',region,''') && ~isempty($SI)'];
            cellsfull = evaluatefilter(cellinfo,filt);
            cells = unique(cellsfull(:,[1 3 4]),'rows');
            
            %days = unique(cells(:,1));
            
            for c = 1:size(cells,1)
                allphases = [];
                allspikes = [];
                day = cells(c,1);
                tet = cells(c,2);
                cell = cells(c,3);
                
                odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
                spikes = loaddatastruct(animDir, animal, 'spikes',day);
                
                epochs = cellsfull(ismember(cellsfull(:,[1 3 4]), [day tet cell],'rows'),2);
                ph = [];
                for ep = epochs'
                    
                    %find beta tet
                    filt = ['strcmp($area,''',region,''') && strcmp($descrip2,''betatet'')'];
                    betatet = evaluatefilter(tetinfo{day}{ep},filt);
                    if isempty(betatet)
                        switch region
                            case 'CA1'
                                ref = 'hpcRef';
                            case 'PFC'
                                ref = 'pfcRef';
                        end
                        filt = ['strcmp($descrip,''',ref,''') && strcmp($descrip2,''betatet'')'];
                        betatet = evaluatefilter(tetinfo{day}{ep},filt);
                    end
                    
                    %get beta
                    beta = loadeegstruct(animDir, animal,'beta',day,ep,betatet);
                    betaphase = beta{day}{ep}{betatet}.data(:,2);
                    betatime  = geteegtimes(beta{day}{ep}{betatet});
                    
                    
                    %create matrix of time bins for each trial
                    npWin = odorTriggers{day}{ep}.allTriggers;
                    npWin = [npWin, npWin+win];
                    %                     timebins = repmat(npWin,1,win/binsize+1);
                    %                     for i = 1:(win/binsize)
                    %                         timebins(:,i+1) = timebins(:,i)+binsize;
                    %                     end
                    
                    %get spike times
                    
                    if ~isempty(spikes{day}{ep}{tet}{cell}.data)
                        spiketimes = spikes{day}{ep}{tet}{cell}.data(:,1);
                        
                        %find spikes that fall during time bins, and the beta phase of those spikes
                        for t = 1:size(npWin)
                            trialspikes = spiketimes(isExcluded(spiketimes,npWin(t,:)));
                            phases = double(betaphase(lookup(trialspikes, betatime)))/10000;
                            times = trialspikes - npWin(t,1);
                            allspikes = [allspikes;times];
                            allphases = [allphases;phases];
                        end
                        
                        %                         binphase = zeros(1,size(timebins,2))-1;
                        %                         for b = 1:size(timebins,2)-1
                        %                             spikesinbin = spiketimes(isExcluded(spiketimes,[timebins(:,b),timebins(:,b+1)]));
                        %                             phases = double(betaphase(lookup(spikesinbin, betatime)))/10000;
                        %
                        %                             allspikes = [allspikes;spikesinbin-timebins(:,1)];
                        %                             allphases2 = [allphases2;phases];
                        %
                        %                             binmean = mean(phases);
                        %                             binphase(b) = binmean;
                        %                         end
                        %                         plot(binphase);
                        %                         pause
                        %                         close all
                        %                         ph = [ph;binphase];
                    else
                        continue
                    end
                    
                    
                end
                %                 mn = nanmean(ph,1);
                %                 allphases = [allphases;mn];
                figure, plot(allspikes,allphases,'k.');
                pause
                close all
            end
            
        end
        meanphases = nanmean(allphases,1);
        ranovatbl = ranova(allphases);
        pause
    end
end