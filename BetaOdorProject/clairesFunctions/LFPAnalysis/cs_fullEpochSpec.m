%cs_FullEpochCoherence
clear
topDir = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
animals = {'CS41','CS42'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS

regions = {'TC'};

freqband = 'low'; %0-40Hz
%freqband = 'floor'; %0-15Hz

%---- Filters -----%
runfilter = 'strcmp($environment,''odorplace'') || strcmp($environment,''novelodor'') || strcmp($environment,''novelodor2'')';

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
params.tapers = [2 1]; %DO NOT GO LOWER THAN THIS- MESSES UP

for r = 1:length(regions)
    region = regions{r};
    
    
    %----- Tetrode Filter -----%
    
    switch region
        case 'CA1'
            %tetfilter = {'(strcmp($area, ''CA1'' ))','(strcmp($area, ''PFC'' ))'};
            tetfilter = 'strcmp($area, ''CA1'' ) && $numcells>0';
        case 'PFC'
            %tetfilter = {'(strcmp($area, ''PFC''))', '(strcmp($area, ''OB'' ))'};
            tetfilter = 'strcmp($area, ''PFC'' ) && $numcells>0';
        case 'OB'
            %tetfilter = {'(strcmp($area, ''CA1'' ))','(strcmp($area, ''OB''))'};
            tetfilter = 'strcmp($area, ''OB'' )';
        case 'TC'
            tetfilter = 'strcmp($area, ''TC'' )';
    end
    
    for a = 1:length(animals)
        animal = animals{a};
        %disp(['Doing ',animal]);
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        tetinfo = loaddatastruct(animDir,animal,'tetinfo');
        
        taskfiles = dir([animDir,animal,'task*']);
        epMatrix = [];
        for f = 1:length(taskfiles)
            load([animDir,taskfiles(f).name]);
            
            ep = evaluatefilter(task, runfilter);
            epMatrix = [epMatrix;ep];
        end
        
        days = unique(epMatrix(:,1));
        
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = epMatrix(epMatrix(:,1) == day,2);
            
            
            tets = evaluatefilter(tetinfo{day},tetfilter);
            tets = unique(tets(:,2));
            if isempty(tets)
                continue
            end
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);

                if isempty(tets)
                    continue
                end
                
                
                
                %if isempty(coherence{day}{ep}.Coh)
                allSpec = [];
                for t = 1:size(tets,1)
                    tet = tets(t);
                    tetstr = getTwoDigitNumber(tet);
    
                    %disp(['Doing ', num2str(tet1), '-',num2str(tet2)]);
                    
                eeg = loadeegstruct(animDir, animal, 'eeg', day, ep, tet);
                time = [eeg{day}{ep}{tet}.starttime:(1/eeg{day}{ep}{tet}.samprate):eeg{day}{ep}{tet}.endtime]';
                eeg = eeg{day}{ep}{tet}.data;
                
                
                %----- Do the specgram calc -----%
                disp(['Doing ',animal,' day ', num2str(day), ' epoch ', num2str(ep), ' tet ',num2str(tet)]);
                [S, Stime, Sfreq] = mtspecgramc(eeg,movingwin,params);
                
                %----- Find the mean and std for whole epoch and do zscore -----%
                load([animDir,'EEGSpec\',animal,'eeggndspec',freqband,daystr,'-',tetstr]) 
                meandayspec = eeggndspec{day}{1}{tet}.meandayspec; % Stored in first epoch
                stddayspec = eeggndspec{day}{1}{tet}.stddayspec;
                
                S = (S - meandayspec)./stddayspec;
                Stime = Stime + time(1);

                allSpec = cat(3,allSpec,S);
                end
                
                Spec = mean(allSpec,3);
                %out.index = [tet1,tet2];
                
                out.Spec = Spec;
                out.time = Stime;
                out.freq = Sfreq;

                spec{day}{ep} = out;
                %end
            end
            save([animDir,animal,'spec',freqband,region,daystr],'spec')
            clear spec
        end

    end
end