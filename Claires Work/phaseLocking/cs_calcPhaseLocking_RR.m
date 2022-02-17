function cs_calcPhaseLocking_RR(a, day, region, windows, trialstr)

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

oldfiles = dir([pldir, animal,'respphaselock','_*','-',region,trialstr,'_',daystr,'.mat']);
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
    %correct trials only
    [cL, cR] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
    resplist = resplist([cL;cR],:);
    
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
        
        tet = cs_getMostCellsTet(animal,day,epoch,region);
        
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
        
        %% ------HP cells----%
        for cell = hpnum:-1:1
            cind = hpidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                %t = geteegtimes(resp{day}{epoch}{reftet});
                sph = respphase(lookup(s, t));
                goodspikes = isExcluded(s, resplist);
                numgoodspikes = sum(goodspikes);
                sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                if length(sph)> 10 %don't bother with stats unless there are at least 10 spikes to use
                    % Rayleigh and Modulation: Originally in lorenlab Functions folder
                    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                    [m, ph] = modulation(sph);
                    phdeg = ph*(180/pi);
                    
                    % Von Mises Distribution - From Circular Stats toolbox
                    [resphat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
%                     if kappa > 1
%                         disp('High Kappa Found')
%                     end
                    resphat_deg = resphat*(180/pi);
                    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                    % Make finer polar plot and overlay Von Mises Distribution Fit.
                    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                    % -------------------------------------------------------------
                    nbins = 50;
                    bins = -pi:(2*pi/nbins):pi;
                    count = histc(sph, bins);
                    
                    % Make Von Mises Fit
                    alpha = linspace(-pi, pi, 50)';
                    [pdf] = circ_vmpdf(alpha,resphat,kappa);
                else
                    stats=0;
                    m=0;
                    phdeg=0;
                    kappa=0;
                    resphat = 0;
                    resphat_deg=0;
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
                resphat = 0;
                resphat_deg=0;
                prayl=NaN;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.Nspikes = numgoodspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.resphat = resphat;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.resphat_deg = resphat_deg;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            resp_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            
            
            
        end
        
        
        
        %% ------CTX cells----%
        
        for cell = ctxnum:-1:1
            cind = ctxidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                t = geteegtimes(resp{day}{epoch}{tet});
                sph = respphase(lookup(s, t));
                goodspikes = isExcluded(s, resplist);
                numgoodspikes = sum(goodspikes);
                sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                if length(sph)> 10
                    % Rayleigh and Modulation: Originally in lorenlab Functions folder
                    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                    [m, ph] = modulation(sph);
                    phdeg = ph*(180/pi);
                    
                    % Von Mises Distribution - From Circular Stats toolbox
                    [resphat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                    resphat_deg = resphat*(180/pi);
                    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                    % Make finer polar plot and overlay Von Mises Distribution Fit.
                    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                    % -------------------------------------------------------------
                    nbins = 50;
                    bins = -pi:(2*pi/nbins):pi;
                    count = histc(sph, bins);
                    
                    % Make Von Mises Fit
                    alpha = linspace(-pi, pi, 50)';
                    [pdf] = circ_vmpdf(alpha,resphat,kappa);
                else
                    stats=0;
                    m=0;
                    phdeg=0;
                    kappa=0;
                    resphat = 0;
                    resphat_deg=0;
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
                resphat = 0;
                resphat_deg=0;
                prayl=NaN;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.Nspikes = numgoodspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.resphat = resphat;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.resphat_deg = resphat_deg;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            resp_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            
            if (figopt == 1)
                
            end
            %     close all
            
            
            %resp_phaselockPFC{day}{epoch} = [];
        end
%     end
end

%% Save
if savedata
    if strcmp(trialstr,'_incorrect')
        if exist('resp_phaselockPFC')
            save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockPFC');
        end
        
        if exist('resp_phaselockCA1')
            save([animDir,'PhaseLocking\Incorrect\',animal,'respphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockCA1');
        end
    else
        
        if exist('resp_phaselockPFC')
            save([animDir,'PhaseLocking\',animal,'respphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockPFC');
        end
        
        if exist('resp_phaselockCA1')
            save([animDir,'PhaseLocking\',animal,'respphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'resp_phaselockCA1');
        end
    end
end



