%cs_FullEpochCoherence
clear
topDir = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
animals = {'CS41','CS42'};

%regions = {'CA1-PFC','CA1-OB','PFC-OB'};
%regions = {'V1-OB'};
regions = {'OB-TC','CA1-TC','PFC-TC'};
%---- Filters -----%
%runepochfilter = 'isequal($environment, ''odorplace'') || isequal($environment, ''novelodor'') || isequal($environment, ''novelodor2'')';

runepochfilter = 'isequal($environment, ''odorplace'')';
freqband = 'low';
%----- Params -----%
movingwin = [1000 20]/1000;
params.Fs = 1500;
params.err = [2 0.05];
switch freqband
    case 'low'
        params.fpass = [0 40];
    case 'floor'
        params.fpass = [0 15];
end

params.tapers = [2 3]; %DO NOT GO LOWER THAN [2 3]- MESSES UP

for r = 1:length(regions)
    region = regions{r};
    
    %----- Tetrode Filter -----%
    
    switch region
        case 'CA1-PFC'
            %tetfilter = {'(strcmp($area, ''CA1'' ))','(strcmp($area, ''PFC'' ))'};
            tetfilter = {'strcmp($area, ''CA1'' ) && $numcells>0' ,'strcmp($area, ''PFC'' ) && $numcells>0'};
        case 'PFC-OB'
            %tetfilter = {'(strcmp($area, ''PFC''))', '(strcmp($area, ''OB'' ))'};
            tetfilter = {'strcmp($area, ''PFC'' ) && $numcells>0',...
                'strcmp($area, ''OB'')'};
        case 'CA1-OB'
            %tetfilter = {'(strcmp($area, ''CA1'' ))','(strcmp($area, ''OB''))'};
            tetfilter = {'strcmp($area, ''CA1'' ) && $numcells>0 ',...
                'strcmp($area, ''OB'')'};
            
        case 'CA1-TC'
            tetfilter = {'(strcmp($area, ''CA1'' )) && $numcells >0','(strcmp($area, ''TC''))'};
        case 'PFC-TC'
            tetfilter = {'(strcmp($area, ''PFC'' )) && $numcells >0','(strcmp($area, ''TC''))'};
        case 'OB-TC'
            tetfilter = {'(strcmp($area, ''OB'' ))','(strcmp($area, ''TC''))'};
    end
    
    for a = 1:length(animals)
        animal = animals{a};
        %disp(['Doing ',animal]);
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        tetinfo = loaddatastruct(animDir,animal,'tetinfo');
        rewards = loaddatastruct(animDir, animal,'rewards');
        
        taskfiles = dir([animDir,animal,'task*']);
        epMatrix = [];
        for f = 1:length(taskfiles)
            load([animDir,taskfiles(f).name]);
            
            ep = evaluatefilter(task, runepochfilter);
            epMatrix = [epMatrix;ep];
        end
        
        days = unique(epMatrix(:,1));
        
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = epMatrix(epMatrix(:,1) == day,2);
            
            try
                load([animDir,animal,'coherence',region,daystr])
            catch
                coherence = [];
            end
            
            tets1 = evaluatefilter(tetinfo{day},tetfilter{1});
            tets1 = unique(tets1(:,2));
            
            %if no tets with cells, just use all tets in region?
            if isempty(tets1)
                newstr = [erase(tetfilter{1}, '$numcells>0'), 'strcmp($descrip2,''betatet'')'];
                tets1 = evaluatefilter(tetinfo{day},newstr);
                tets1 = unique(tets1(:,2));
            end
            
            tets2 = evaluatefilter(tetinfo{day},tetfilter{2});
            tets2 = unique(tets2(:,2));
            if isempty(tets2)
                newstr = [erase(tetfilter{2}, '$numcells>0'), 'strcmp($descrip2,''betatet'')'];
                tets2 = evaluatefilter(tetinfo{day},newstr);
                tets2 = unique(tets2(:,2));
            end
            
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                %                 if isempty(tets1)
                %                     tetfilter{1} = '(strcmp($area, ''CA1'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''hpcRef'' ) && strcmp($descrip2, ''betatet'' ))';
                %                     tets1 = evaluatefilter(tetinfo{day}{ep},tetfilter{1});
                %                 end
                %
                %                 if isempty(tets2)
                %                     tetfilter{2} = '(strcmp($area, ''PFC'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''pfcRef'' ) && strcmp($descrip2, ''betatet'' ))';
                %                     tets2 = evaluatefilter(tetinfo{day}{ep},tetfilter{2});
                %                 end
                
                %                 if isempty(tets1) || isempty(tets2)
                %                     continue
                %                 end
                
                [A,B] = meshgrid(tets1,tets2);
                pairs=reshape(cat(2,A',B'),[],2);
                
                %if isempty(coherence{day}{ep}.Coh)
                allPairsCoh = [];
                for p = 1:size(pairs,1)
                    tet1 = pairs(p,1);
                    tet2 = pairs(p,2);
                    
                    %disp(['Doing ', num2str(tet1), '-',num2str(tet2)]);
                    
                    eeg1 = loadeegstruct(animDir, animal, 'eeg', day, ep, tet1);
                    time1 = [eeg1{day}{ep}{tet1}.starttime:(1/eeg1{day}{ep}{tet1}.samprate):eeg1{day}{ep}{tet1}.endtime]';
                    eeg1 = eeg1{day}{ep}{tet1}.data;
                    
                    eeg2 = loadeegstruct(animDir, animal, 'eeg', day, ep, tet2);
                    time2 = [eeg2{day}{ep}{tet2}.starttime:(1/eeg2{day}{ep}{tet2}.samprate):eeg2{day}{ep}{tet2}.endtime]';
                    eeg2 = eeg2{day}{ep}{tet2}.data;
                    
                    
                    if ~all(time2 == time1)% check if time vectors are matched
                        eeg2 = interp1(time2,eeg2,time1,'nearest');
                    end
                    
                    %----- Do the cohereogram calc -----%
                    disp(['Doing ',animal,' day ', num2str(day), ' epoch ', num2str(ep), ' tets ',num2str(tet1), '/', num2str(tet2)]);
                    [Coh,Phi,~,~,~,t,freq] = cohgramc(eeg1,eeg2,movingwin,params);
                    Coh = Coh';
                    t = t + time1(1);
                    fs_c = round(1/(t(2)-t(1))); %length of time bins
                    
                    %----- Find the mean and std for whole epoch for zscore later -----%
                    meanCohEpoch = mean(Coh,2);
                    sdCohEpoch = std(Coh,0,2);
                    
                    %----- Normalize with Zscore -----%
                    epochCoh = (Coh - meanCohEpoch)./sdCohEpoch;
                    allPairsCoh = cat(3,allPairsCoh,epochCoh);
                end
                
                Coh = mean(allPairsCoh,3);
                %out.index = [tet1,tet2];
                
                out.Coh = Coh;
                out.time = t;
                out.freq = freq;
                out.mean = meanCohEpoch;
                out.sd = sdCohEpoch;
                
                coherence{day}{ep} = out;
                %end
            end
            save([animDir,animal,'coherence',region,'_',freqband,'_',daystr],'coherence')
            clear coherence
        end
        
    end
end