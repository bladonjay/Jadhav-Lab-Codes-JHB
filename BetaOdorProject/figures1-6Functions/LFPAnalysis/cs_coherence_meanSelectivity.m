%cs_coherence_meanSelectivity

%get correlation between coherence and the selectivity in the population
%use all odor responsive cells

clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
freqs = {'beta','resp'};
cellregions = {'CA1','PFC'};

for d = 1:length(freqs)
    
    freq = freqs{d};
    switch freq
        case 'beta'
            bandpass = [15 30];
        case 'resp'
            bandpass = [7 8];
    end
    
    for c = 1:length(cellregions)
        region = cellregions{c};
        
        meanSelectivity = [];
        meanCoh = [];
        for a = 1:length(animals)
            animal = animals{a};
            disp(['Doing ',animal,' ',region])
            animDir = [topDir, animal,'Expt\',animal,'_direct\'];
            
            cellinfo = loaddatastruct(animDir, animal,'cellinfo');
            odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
            dayEps = cs_getRunEpochs(animDir, animal,'odorplace');
            days = unique(dayEps(:,1));
            
            for day = days'
                daystr = getTwoDigitNumber(day);
                eps = cs_getRunEpochs(animDir, animal,'odorplace',day);
                eps = eps(:,2);
                
                %get cell selectivity
                cellfilt = ['(strcmp($area, ''',region,''')) && strcmp($type,''pyr'') && (~isempty($SI))'];
                cells = evaluatefilter(cellinfo{day},cellfilt);
                cells = unique(cells(:,[2,3]),'rows');
                
                SI = [];
                for c = 1:size(cells,1)
                    cell = cells(c,:);
                    eps = cs_findGoodEpochs(cellinfo{day},{'SI'},cell);
                    si = abs(cellinfo{day}{eps(1)}{cell(1)}{cell(2)}.SI);
                    if isnan(si) || isempty(si)
                        keyboard
                    end
                    SI = [SI;si];
                    
                end
                if isnan(nanmean(SI))
                    keyboard
                end
                meanSelectivity = [meanSelectivity;nanmean(SI)];
                
                %get coherence
                load([animDir, animal,'coherenceCA1-PFC',daystr,'.mat']);
                
                dayCoh = [];
                for ep = eps'
                    [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    coh = coherence{day}{ep}.npWinCoh(:,[cl;cr]);
                    
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    coh = coh(goodrows,:);
                    coh = mean(coh,1);
                    dayCoh = [dayCoh,coh];
                end
                
                meanCoh = [meanCoh;mean(dayCoh)];
            end
            
        end
        figure,
        plot(meanCoh, meanSelectivity,'k.');
        [CC,p] = corrcoef(meanCoh, meanSelectivity);

        
    end
end

