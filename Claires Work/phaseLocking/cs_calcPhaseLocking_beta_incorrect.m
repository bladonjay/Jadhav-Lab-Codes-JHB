function cs_calcPhaseLocking_beta_incorrect(a, day, region, windows)
trialstr = '_incorrect';
if length(windows) < day || isempty(windows{day})
    return
end
%% Params
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animal = animals{a};

animDir = [topDir, animal,'Expt\',animal,'_direct\'];
daystr = getTwoDigitNumber(day);


close all;
savedata = 1;
figopt = 0;

%% Delete all existing files

pldir = [animDir,'PhaseLocking\Incorrect\'];


oldfiles = dir([pldir, animal,'betaphaselock','_*','-',region,'_incorrect_',daystr,'.mat']);
cd (pldir)
if length(oldfiles) > 0
    delete (oldfiles.name)
end

%% Load Data
load([animDir, animal, 'tetinfo.mat'], 'tetinfo');
load([animDir, animal, 'cellinfo.mat'], 'cellinfo');
load([animDir, animal, 'odorTriggers',daystr,'.mat'],'odorTriggers');
load([animDir, animal, 'spikes', daystr,'.mat'], 'spikes'); % get spikes

epochs = find(~cellfun(@isempty,windows{day}));

for ep = 1:length(epochs)
    epoch = epochs(ep);
    epstr = getTwoDigitNumber(epoch);
    
    %% Get Time Windows
    
    betalist = windows{day}{epoch};
    
    
    if ~isstruct(betalist) || ~isempty(betalist.OB)
        if  isstruct(betalist)%there may be no high beta times for this region and epoch
            betalist = betalist.OB;
        end
       
        %betalist = betalist(~isnan(betalist(:,1)),:);
        
        %% Get Cells
        
        if strcmp(trialstr,'_prelearn') || strcmp(trialstr,'_postlearn')
            load([topDir,'AnalysesAcrossAnimals\npCells_novel_CA1.mat']);
            hpidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
            
            load([topDir,'AnalysesAcrossAnimals\npCells_novel_PFC.mat']);
            ctxidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
        else
            load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat']);
            hpidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
            
            load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat']);
            ctxidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
        end
        %                 cellfilter = 'isequal($area,''CA1'') && strcmp($tag, ''accepted'')' ;
        %                 CA1cells = evaluatefilter(cellinfo{day}{epoch},cellfilter);
        %                 cellfilter = 'isequal($area,''PFC'') && strcmp($tag, ''accepted'')' ;
        %                 PFCcells = evaluatefilter(cellinfo{day}{epoch},cellfilter);
        %
        %                 ctxidx = intersect(ctxidx, PFCcells, 'rows');
        %                 hpidx = intersect(hpidx, CA1cells, 'rows');
        
        
        % %
        hpnum = size(hpidx(:,1),1);
        ctxnum = size(ctxidx(:,1),1);
        
        
        
        %% -----load beta eeg using betatet-----%
        
%         tetfilter = ['((strcmp($area, ''',region,''')) && (strcmp($descrip2, ''betatet'')))'];
%         %tetfilter = ['((strcmp($area, ''OB'')) && (strcmp($descrip2, ''betatet'')))'];
%         
%         tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);

        tet = cs_getMostCellsTet(animal, day, ep, region);
        if isempty(tet)
            beta_phaselockPFC{day}{epoch} = [];
            beta_phaselockCA1{day}{epoch} = [];
            continue
        end
        
        tetstr = getTwoDigitNumber(tet);
        
        %tmpflist1 = sprintf('%s%sbeta%02d-%02d-%02d.mat', [animDir,'EEG/'],animal, day, epoch, reftet);
        load([animDir, 'EEG/', animal, 'beta',daystr,'-',epstr,'-',tetstr,'.mat'],'beta');
        t = geteegtimes(beta{day}{epoch}{tet});
        betaphase = beta{day}{epoch}{tet}.data(:,2);% beta phase
        
        %highBeta_correct = loaddatastruct(animDir, animal,'highBeta');
        %load([animDir, animal,'highBeta']);
        
        numtrials = size(betalist,1);
        
        %% Bootstrap for incorrect trials to match number with correct trials
        iterations = 1000;
        Idist = [];
        
        %% ------HP cells----%
        for cell = hpnum:-1:1
            cind = hpidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                goodspikes = s(isExcluded(s, betalist));
                sph_real = double(betaphase(lookup(goodspikes, t)))/10000;  
                
                kappa_dist= [];
                zrayl_dist = [];
                
                for i = 1:iterations
                    samp = datasample(1:size(betalist,1),numtrials);
                    trials = betalist(samp,:);
                    
                    sph = betaphase(lookup(s, t));
                    goodspikes = isExcluded(s, trials);
                    numgoodspikes = sum(goodspikes);
                    sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                    if length(sph)> 10 %don't bother with stats unless there are at least 10 spikes to use
                        
                        % Von Mises Distribution - From Circular Stats toolbox
                        [~, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                        [~, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                        
                    else
                        kappa=0;
                        zrayl=0;
                    end
                    
                    kappa_dist = [kappa_dist;kappa];
                    zrayl_dist = [zrayl_dist;zrayl];
                end
                
            else
                kappa_dist=0; 
                zrayl_dist=0;
                sph_real = [];

            end
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa_dist;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl_dist;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.sph = sph_real; %actual spike phases, not bootstrapped
            
        end
        
        
        
        %% ------CTX cells----%
        
        for cell = ctxnum:-1:1
            cind = ctxidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                goodspikes = s(isExcluded(s, betalist));
                sph_real = double(betaphase(lookup(goodspikes, t)))/10000;
                
                kappa_dist= [];
                zrayl_dist = [];
                
                for i = 1:iterations
                    samp = datasample(1:size(betalist,1),numtrials);
                    trials = betalist(samp,:);
                    
                    sph = betaphase(lookup(s, t));
                    goodspikes = isExcluded(s, trials);
                    numgoodspikes = sum(goodspikes);
                    sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                    if length(sph)> 10 %don't bother with stats unless there are at least 10 spikes to use
                        
                        % Von Mises Distribution - From Circular Stats toolbox
                        [~, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                        [~, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                        
                    else
                        kappa=0;
                        zrayl=0;
                        sph_real = [];
                    end
                    
                    kappa_dist = [kappa_dist;kappa];
                    zrayl_dist = [zrayl_dist;zrayl];
                end
                
            else
                kappa_dist=0; 
                zrayl_dist=0;
                sph_real = [];
            end
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa_dist;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl_dist;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.sph = sph_real;
        end
    end
end

%% Save
if savedata
    if exist('beta_phaselockPFC')
        save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockPFC');
    else
        beta_phaselockPFC{day} = [];
                save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockPFC');    
    end

    if exist('beta_phaselockCA1')
        save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockCA1');
    else
        beta_phaselockCA1{day} = [];
                save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockCA1');

    end
end



