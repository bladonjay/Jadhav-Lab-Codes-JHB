function cs_calcPhaseLocking_RR_incorrect(a, day, region, windows)
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


oldfiles = dir([pldir, animal,'respphaselock','_*','-',region,'_incorrect_',daystr,'.mat']);
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
    
    resplist = windows{day}{epoch};
    
    [~,~,iL, iR] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
    resplist = resplist([iL;iR],:);
    
    
%     if ~isstruct(resplist) || ~isempty(resplist.OB)
%         if  isstruct(resplist)%there may be no high resp times for this region and epoch
%             resplist = resplist.OB;
%         end
       
        %resplist = resplist(~isnan(resplist(:,1)),:);
        
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
        
        
        
        %% -----load resp eeg using resptet-----%
        
%         tetfilter = ['((strcmp($area, ''',region,''')) && (strcmp($descrip2, ''resptet'')))'];
%         %tetfilter = ['((strcmp($area, ''OB'')) && (strcmp($descrip2, ''resptet'')))'];
%         
%         tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);

        tet = cs_getMostCellsTet(animal, day, ep, region);
        if isempty(tet)
            resp_phaselockPFC{day}{epoch} = [];
            resp_phaselockCA1{day}{epoch} = [];
            continue
        end
        
        tetstr = getTwoDigitNumber(tet);
        
        %tmpflist1 = sprintf('%s%sresp%02d-%02d-%02d.mat', [animDir,'EEG/'],animal, day, epoch, reftet);
        load([animDir, 'EEG/', animal, 'resp',daystr,'-',epstr,'-',tetstr,'.mat'],'resp');
        t = geteegtimes(resp{day}{epoch}{tet});
        respphase = resp{day}{epoch}{tet}.data(:,2);% resp phase
        
        %highresp_correct = loaddatastruct(animDir, animal,'highresp');
        %load([animDir, animal,'highresp']);
        
        numtrials = size(resplist,1);
        
        %% Bootstrap for incorrect trials to match number with correct trials
        iterations = 1000;
        Idist = [];
        
        %% ------HP cells----%
        for cell = hpnum:-1:1
            cind = hpidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                goodspikes = s(isExcluded(s, resplist));
                sph_real = double(respphase(lookup(goodspikes, t)))/10000;  
                
                kappa_dist= [];
                zrayl_dist = [];
                
                for i = 1:iterations
                    samp = datasample(1:size(resplist,1),numtrials);
                    trials = resplist(samp,:);
                    
                    sph = respphase(lookup(s, t));
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
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa_dist;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl_dist;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.sph = sph_real; %actual spike phases, not bootstrapped
            
        end
        
        
        
        %% ------CTX cells----%
        
        for cell = ctxnum:-1:1
            cind = ctxidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                goodspikes = s(isExcluded(s, resplist));
                sph_real = double(respphase(lookup(goodspikes, t)))/10000;
                
                kappa_dist= [];
                zrayl_dist = [];
                
                for i = 1:iterations
                    samp = datasample(1:size(resplist,1),numtrials);
                    trials = resplist(samp,:);
                    
                    sph = respphase(lookup(s, t));
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
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa_dist;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl_dist;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.sph = sph_real;
        end
    %end
end

%% Save
if savedata
    if exist('resp_phaselockPFC')
        save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockPFC');
    else
        resp_phaselockPFC{day} = [];
                save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockPFC');    
    end

    if exist('resp_phaselockCA1')
        save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockCA1');
    else
        resp_phaselockCA1{day} = [];
                save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockCA1');

    end
end



