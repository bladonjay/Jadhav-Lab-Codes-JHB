function cs_calcPhaseLocking_beta(a, day, region, windows, trialstr)

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
if strcmp(trialstr,'_incorrect')
    pldir = [animDir,'PhaseLocking\Incorrect\'];
else
    pldir = [animDir,'PhaseLocking\'];
end

oldfiles = dir([pldir, animal,'betaphaselock','_*','-',region,trialstr,'_',daystr,'.mat']);
cd (pldir)
if length(oldfiles) > 0
    delete (oldfiles.name)
end

%% Load Data
load([animDir, animal, 'tetinfo.mat'], 'tetinfo');
load([animDir, animal, 'cellinfo.mat'], 'cellinfo');
load([animDir, animal, 'odorTriggers',daystr,'.mat'],'odorTriggers');
load([animDir, animal,'nosepokeWindow',daystr,'.mat'],'nosepokeWindow');
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
        betalist = nosepokeWindow{day}{epoch};
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
        
        %tetfilter = ['((strcmp($area, ''',region,''')) && (strcmp($descrip2, ''betatet'')))'];
        %tetfilter = ['((strcmp($area, ''OB'')) && (strcmp($descrip2, ''betatet'')))'];
        
        %tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);
        
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
        
        %% ------HP cells----%
        for cell = hpnum:-1:1
            cind = hpidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                %t = geteegtimes(beta{day}{epoch}{reftet});
                sph = betaphase(lookup(s, t));
                goodspikes = isExcluded(s, betalist);
                numgoodspikes = sum(goodspikes);
                sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                if length(sph)> 10 %don't bother with stats unless there are at least 10 spikes to use
                    % Rayleigh and Modulation: Originally in lorenlab Functions folder
                    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                    [m, ph] = modulation(sph);
                    phdeg = ph*(180/pi);
                    
                    % Von Mises Distribution - From Circular Stats toolbox
                    [betahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
%                     if kappa > 1
%                         disp('High Kappa Found')
%                     end
                    betahat_deg = betahat*(180/pi);
                    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                    % Make finer polar plot and overlay Von Mises Distribution Fit.
                    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                    % -------------------------------------------------------------
                    nbins = 50;
                    bins = -pi:(2*pi/nbins):pi;
                    count = histc(sph, bins);
                    
                    % Make Von Mises Fit
                    alpha = linspace(-pi, pi, 50)';
                    [pdf] = circ_vmpdf(alpha,betahat,kappa);
                else
                    stats=0;
                    m=0;
                    phdeg=0;
                    kappa=0;
                    betahat = 0;
                    betahat_deg=0;
                    prayl=NaN;
                    zrayl=0;
                    alpha=0;
                    pdf=0;
                end
                
            else
                sph = [];
                numgoodspikes = 0;
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                betahat = 0;
                betahat_deg=0;
                prayl=NaN;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.Nspikes = numgoodspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.betahat = betahat;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.betahat_deg = betahat_deg;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            
            
            
        end
        
        
        
        %% ------CTX cells----%
        
        for cell = ctxnum:-1:1
            cind = ctxidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                t = geteegtimes(beta{day}{epoch}{tet});
                sph = betaphase(lookup(s, t));
                goodspikes = isExcluded(s, betalist);
                numgoodspikes = sum(goodspikes);
                sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                if length(sph)> 10
                    % Rayleigh and Modulation: Originally in lorenlab Functions folder
                    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                    [m, ph] = modulation(sph);
                    phdeg = ph*(180/pi);
                    
                    % Von Mises Distribution - From Circular Stats toolbox
                    [betahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                    betahat_deg = betahat*(180/pi);
                    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                    % Make finer polar plot and overlay Von Mises Distribution Fit.
                    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                    % -------------------------------------------------------------
                    nbins = 50;
                    bins = -pi:(2*pi/nbins):pi;
                    count = histc(sph, bins);
                    
                    % Make Von Mises Fit
                    alpha = linspace(-pi, pi, 50)';
                    [pdf] = circ_vmpdf(alpha,betahat,kappa);
                else
                    stats=0;
                    m=0;
                    phdeg=0;
                    kappa=0;
                    betahat = 0;
                    betahat_deg=0;
                    prayl=NaN;
                    zrayl=0;
                    alpha=0;
                    pdf=0;
                end
                
            else
                sph = [];
                numgoodspikes = 0;
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                betahat = 0;
                betahat_deg=0;
                prayl=NaN;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.Nspikes = numgoodspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.betahat = betahat;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.betahat_deg = betahat_deg;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            
            if (figopt == 1)
                
            end
            %     close all
            
            
            %beta_phaselockPFC{day}{epoch} = [];
        end
    end
end

%% Save
if savedata
    if strcmp(trialstr,'_incorrect')
        if exist('beta_phaselockPFC')
            save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockPFC');
        end
        
        if exist('beta_phaselockCA1')
            save([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockCA1');
        end
    else
        
        if exist('beta_phaselockPFC')
            save([animDir,'PhaseLocking\',animal,'betaphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockPFC');
        end
        
        if exist('beta_phaselockCA1')
            save([animDir,'PhaseLocking\',animal,'betaphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'beta_phaselockCA1');
        end
    end
end



