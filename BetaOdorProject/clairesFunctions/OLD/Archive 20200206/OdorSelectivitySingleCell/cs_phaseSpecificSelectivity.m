%% Params
%topDir = 'F:\Data\OdorPlaceAssociation\';
topDir = 'D:\OdorPlaceAssociation\';

animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC'};

% animals = {'CS31','CS33'};
% regions = {'CA1'};
win = [0 1];
iterations = 1000;
winsize = win(2) + win(1);

for r = 1:length(regions)
    region = regions{r};
    
    %% Get preferred phase
    load([topDir,'AnalysesAcrossAnimals\populationPhaseLocking_',region,'.mat'])
    
    
         phaserange = [prefPhase- pi/4, prefPhase+ pi/4];
         nprefphaserange = phaserange + pi;
         phaserange = wrapTo2Pi(phaserange);  
         nprefphaserange = wrapTo2Pi(nprefphaserange);
         
    %must be greater than first value and less than second value. Can wrap
    %around now. 
    %%
    selmatrix = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        load([animDir, animal,'cellinfo.mat']);
        load([animDir, animal,'tetinfo.mat']);
        
        cellfilter = ['((isequal($area,''',region, ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
        
        cells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            
            trigs_all= []; correct_left_all = []; correct_right_all = [];
            for ep = 1:length(runeps)
                epoch = runeps(ep);
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                trigs = odorTriggers{day}{epoch}.allTriggers;
                
                correct_left_all = [correct_left_all; trigs(correct_left)];
                correct_right_all = [correct_right_all; trigs(correct_right)];
                trigs_all = [trigs_all ; trigs];
            
            end
            
            
            for c = 1:size(daycells,1)
                
                disp(['Doing ', animal,' day ',daystr,' cellnum ',num2str(c)]);
                leftspikespref = []; 
                rightspikespref = [];
                leftspikesNpref = []; 
                rightspikesNpref = [];
%                 incorrectleftspikes = [];
%                 incorrectrightspikes = [];
                
                cell = daycells(c,:);
                
                runspikes = []; phases = [];
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    epstr = getTwoDigitNumber(epoch);
                    
                    
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        
                        %lookup spike phases
                        tetfilter = ['isequal($descrip,''hpcRef'')'];
                        reftet = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                        tetstr = getTwoDigitNumber(reftet);
                        load([animDir,'EEG\',animal,'beta',daystr,'-',epstr,'-',tetstr,'.mat'])

                        btime = geteegtimes(beta{day}{epoch}{reftet});
                        bph = beta{day}{epoch}{reftet}.data(:,2);% beta phase
                        %bph = bph/10000;
                    
                        s = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        sph = bph(lookup(s, btime)); 
                        sph = double(sph)/10000;
                        sph = wrapTo2Pi(sph);
                        
                        
                        runspikes = [runspikes; s];
                        phases = [phases; sph];
                                                    
                    end
                end
 %% Only take spikes in win that fall within preferred phase range
 
                for t = 1:length(correct_left_all)
                    trig = correct_left_all(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikephase = phases(runspikes > trigwin(1) & runspikes <= trigwin(2));
                        if phaserange(2) > phaserange(1) 
                            prefspikes = (winspikephase > phaserange(1) & winspikephase < phaserange(2));
                        else
                            prefspikes = (winspikephase < phaserange(1) | winspikephase > phaserange(2));
                        end
                        
                        if nprefphaserange(2) > nprefphaserange(1)
                            nprefspikes = (winspikephase > nprefphaserange(1) & winspikephase < nprefphaserange(2));
                        else
                            nprefspikes = (winspikephase < nprefphaserange(1) | winspikephase > nprefphaserange(2));
                        end
                        
                        
                        prefwinspikes = sum(prefspikes); 
                        nprefwinspikes = sum(nprefspikes);
                        leftspikespref = [leftspikespref; prefwinspikes]; 
                        leftspikesNpref = [leftspikesNpref; nprefwinspikes]; 
                end
                
                for t = 1:length(correct_right_all)
                    trig = correct_right_all(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikephase = phases(runspikes > trigwin(1) & runspikes <= trigwin(2));
                        if phaserange(2) > phaserange(1) 
                            prefspikes = (winspikephase > phaserange(1) & winspikephase < phaserange(2));
                            nprefspikes = (winspikephase < phaserange(1) | winspikephase > phaserange(2));
                        else
                            nprefspikes = (winspikephase > phaserange(1) & winspikephase < phaserange(2));
                            prefspikes = (winspikephase < phaserange(1) | winspikephase > phaserange(2));
                        end
                        
                        prefwinspikes = sum(prefspikes); 
                        nprefwinspikes = sum(nprefspikes);
                        rightspikespref = [rightspikespref; prefwinspikes]; 
                        rightspikesNpref = [rightspikesNpref; nprefwinspikes]; 
                end

 %% Calculate Selectivity 

                %Preferred Phase
                    prefleftfr = mean(leftspikespref)/winsize;
                    prefrightfr = mean(rightspikespref)/winsize;
                    
                    prefphaseSelectivity = (prefleftfr - prefrightfr)/(prefleftfr + prefrightfr);
                    %selectivity = cellinfo{cell(1)}{2}{cell(2)}{cell(3)}.SI;
                    
                 %Non-Preferred Phase
                    nprefleftfr = mean(leftspikesNpref)/winsize;
                    nprefrightfr = mean(rightspikesNpref)/winsize;
                    
                    nprefphaseSelectivity = (nprefleftfr - nprefrightfr)/(nprefleftfr + nprefrightfr);
                    
                    selmatrix(end+1,:) = [prefphaseSelectivity, nprefphaseSelectivity];
                    
%           
                    end
                end
        
        
    end
    selmatrix = abs(selmatrix);
    figure, hold on
    for p = 1:size(selmatrix,1)
        plot([1,2], selmatrix(p,:), 'k-')
    end
    axis([0.5 2.5 0 1])

    p = signrank(selmatrix(:,1),selmatrix(:,2))
end
    