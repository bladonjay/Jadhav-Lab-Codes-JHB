function columnVectors = cs_columnVectors_randsamp(animals, region, win, nsamp, selectiveonly, savedata) 

%This function creates column vectors for PCA, but selects a random sample
%of trials (with replacement) for each epoch. This results in the same
%number of trials across all cells, so that vectors can be concatenated
%easily 

[topDir,~] = cs_setPaths();
dataDir = [topDir, 'AnalysesAcrossAnimals\'];
winstr = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
winsize = win(2)+win(1);        
columnVectors = [];
for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        %cellfilter = ['isequal($area,''',region,'''))'];
        
        if selectiveonly == 1
            cellfilter = ['((isequal($area,''',region, ... 
            ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
            
            selstr = 'selective';
        else
            cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region,''')))'];
            
            selstr = '';
        end
        
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
            
            
            %Gather trig times 
            %trigs_all= []; 
            correct_left_all = []; correct_right_all = [];
            for ep = 1:length(runeps)
                epoch = runeps(ep);
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                trigs = odorTriggers{day}{epoch}.allTriggers;
                
            correct_left_all = [correct_left_all; trigs(correct_left)];
            correct_right_all = [correct_right_all; trigs(correct_right)];
            %trigs_all = [trigs_all ; trigs];
            
            end
            
            %select random sample of trials
            samp = datasample(1:size(correct_left_all,1),nsamp);
            trials_L = correct_left_all(samp);
            
            samp = datasample(1:size(correct_right_all,1),nsamp);
            trials_R = correct_right_all(samp);
            
            for c = 1:size(daycells,1)

                leftspikes = []; 
                rightspikes = [];              
                
                cell = daycells(c,:);

                runspikes = []; 
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];                             
                    end
                    
                    
                end
                
                %Find spikes in trigger windows
                for t = 1:length(trials_L)
                        trig = trials_L(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikes =runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
                        leftspikes = [leftspikes;length(winspikes)];
                end
                
                for t = 1:length(trials_R)
                        trig = trials_R(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikes =runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
                        rightspikes = [rightspikes;length(winspikes)];
                end
                
                
                sumspikes = sum(sum([leftspikes; rightspikes]));
                
                %if sumspikes >= 100 %only take the cell if its total spikes during all trig windows is greater than total number of trials. (avg one spike per trial)
                    leftfr = leftspikes/winsize;
                    rightfr = rightspikes/winsize;
                    
                    fr = [leftfr;rightfr];
                    
                    columnVectors = [columnVectors,fr];
%                     totalcells = totalcells+1;
%                     cellinds{totalcells,1} = [a, cell];
                    
%                     lefttrials{totalcells,1}= leftspikes;
%                     righttrials{totalcells,1} = rightspikes;
                %end
                
            end
        end
end
 
if savedata == 1
     filename = ['columnVectors_randsamp_',region, '_', winstr,'_',selstr];
    
     save([dataDir,filename],'columnVectors')
end