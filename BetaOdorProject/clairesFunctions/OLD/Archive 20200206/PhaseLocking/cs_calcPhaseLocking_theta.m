function cs_calcPhaseLocking_theta(a, day, region, windows, trialstr)

%% Params
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animal = animals{a};

animDir = [topDir, animal,'Expt\',animal,'_direct\'];

close all;
savedata = 1;
figopt = 0;

daystr = getTwoDigitNumber(day);

%% Load Data
load([animDir, animal, 'tetinfo.mat'], 'tetinfo');
load([animDir, animal, 'cellinfo.mat'], 'cellinfo');
load([animDir, animal, 'spikes', daystr,'.mat'], 'spikes'); % get spikes

epochs = find(~cellfun(@isempty,windows{day}));

% %% Get NP Cells
% load([topDir, 'AnalysesAcrossAnimals\npCells_CA1.mat'])
% hpcells = npCells(all(npCells(:,1:2) == [animnum, day],2),[3,4]);
% 
% load([topDir, 'AnalysesAcrossAnimals\npCells_PFC.mat'])
% ctxcells = npCells(all(npCells(:,1:2) == [animnum, day],2),[3,4]);

for ep = 1:length(epochs)
    epoch = epochs(ep);
    epstr = getTwoDigitNumber(epoch);
    
    %% Get Time Windows
    
    thetalist = windows{day}{epoch};
    if ~isstruct(thetalist) || (isstruct(thetalist) && ~isempty(thetalist.(region)))
        if  isstruct(thetalist)%there may be no high beta times for this region and epoch
            thetalist = thetalist.(region);
        end
        
        
        %% Get Cells
        load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat']);
        hpidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
        
        load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat']);
        ctxidx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
        
%         cellfilter = 'isequal($area,''CA1'') && strcmp($type, ''pyr'') && ($numspikes > 100)' ;
%         hpidx = evaluatefilter(cellinfo{day}{epoch},cellfilter);
%         
%         cellfilter = 'isequal($area,''PFC'') && strcmp($type, ''pyr'') && ($numspikes > 100)' ;
%         ctxidx = evaluatefilter(cellinfo{day}{epoch},cellfilter);
        
        %hpidx = hpidx(ismember(hpcells,hpidx,'rows'));
        
        
        hpnum = size(hpidx(:,1));
        ctxnum = size(ctxidx(:,1));
        
        
        
        %% -----load theta eeg using rrtet-----%
 
        tetfilter = '(strcmp($descrip3, ''rrtet''))';
        tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);
        if isempty(tet)
            rr_phaselockPFC{day}{epoch} = [];
            rr_phaselockCA1{day}{epoch} = [];
            continue
        end

        tetstr = getTwoDigitNumber(tet);
        
        %tmpflist1 = sprintf('%s%sbeta%02d-%02d-%02d.mat', [animDir,'EEG/'],animal, day, epoch, reftet);
        load([animDir, 'EEG/', animal, 'theta',daystr,'-',epstr,'-',tetstr,'.mat'],'theta');
        t = geteegtimes(theta{day}{epoch}{tet});
        thetaphase = theta{day}{epoch}{tet}.data(:,2);% beta phase
        
        %% ------HP cells----%
        for cell = hpnum:-1:1
            cind = hpidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
            s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
            %t = geteegtimes(beta{day}{epoch}{reftet});
            sph = thetaphase(lookup(s, t));
            goodspikes = isExcluded(s, thetalist);
            Nspikes = sum(goodspikes);
            sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
            if length(sph)>size(thetalist,1) %number of spikes has to be greater than number of trials
                % Rayleigh and Modulation: Originally in lorenlab Functions folder
                stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                [m, ph] = modulation(sph);
                phdeg = ph*(180/pi);
                
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                thetahat_deg = thetahat*(180/pi);
                [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                % Make finer polar plot and overlay Von Mises Distribution Fit.
                % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                % -------------------------------------------------------------
                nbins = 50;
                bins = -pi:(2*pi/nbins):pi;
                count = histc(sph, bins);
                
                % Make Von Mises Fit
                alpha = linspace(-pi, pi, 50)';
                [pdf] = circ_vmpdf(alpha,thetahat,kappa);
            else
                sph = [];
                Nspikes = 0;
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                thetahat = 0;
                thetahat_deg=0;
                prayl=0;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.Nspikes = Nspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.betahat = thetahat;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.betahat_deg = thetahat_deg;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            rr_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            
            else
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                betahat = 0;
                betahat_deg=0;
                prayl=0;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            if (figopt == 1)
                %CS - removed plotting code to clean it up. Run wenbo's
                %original code if plotting individual cells is necessary.
            end
            
            
        end
        
        
        
        %% ------CTX cells----%
        
        for cell = ctxnum:-1:1
            cind = ctxidx(cell,:);
            if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
            s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
            t = geteegtimes(theta{day}{epoch}{tet});
            sph = thetaphase(lookup(s, t));
            goodspikes = isExcluded(s, thetalist);
            Nspikes = sum(goodspikes);
            sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
            if length(sph)>size(thetalist,1) 
                % Rayleigh and Modulation: Originally in lorenlab Functions folder
                stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                [m, ph] = modulation(sph);
                phdeg = ph*(180/pi);
                
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                thetahat_deg = thetahat*(180/pi);
                [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                % Make finer polar plot and overlay Von Mises Distribution Fit.
                % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                % -------------------------------------------------------------
                nbins = 50;
                bins = -pi:(2*pi/nbins):pi;
                count = histc(sph, bins);
                
                % Make Von Mises Fit
                alpha = linspace(-pi, pi, 50)';
                [pdf] = circ_vmpdf(alpha,thetahat,kappa);
            else
                sph = [];
                Nspikes = 0;
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                thetahat = 0;
                thetahat_deg=0;
                prayl=0;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = tet;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.Nspikes = Nspikes;
            % output stats also, although you will have to redo this after combining epochochs
            % Rayleigh test
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.modln = m;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
            % Von Mises Fit
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.betahat = thetahat;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.betahat_deg = thetahat_deg;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
            rr_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
            else
                stats=0;
                m=0;
                phdeg=0;
                kappa=0;
                betahat = 0;
                betahat_deg=0;
                prayl=0;
                zrayl=0;
                alpha=0;
                pdf=0;
            end
            
            if (figopt == 1)
                
            end
            %     close all
            
            
        end
    else
        rr_phaselockPFC{day}{epoch} = [];
        rr_phaselockCA1{day}{epoch} = [];
    end
end

%% Save
if savedata
    if exist('rr_phaselockPFC')
    save([animDir,'PhaseLocking\',animal,'rrphaselock_PFC-',region,trialstr,'_',daystr,'.mat'],'rr_phaselockPFC');
    end
    
    if exist('rr_phaselockCA1')
    save([animDir,'PhaseLocking\',animal,'rrphaselock_CA1-',region,trialstr,'_',daystr,'.mat'],'rr_phaselockCA1');
    end
end



