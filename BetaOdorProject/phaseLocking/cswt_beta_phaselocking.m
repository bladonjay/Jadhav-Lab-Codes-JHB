function cswt_beta_phaselocking(animal, day, odorTriggers)

topDir = cs_setPaths();

animDir = [topDir, animal,'Expt\',animal,'_direct\'];

close all;
savedata = 1;
figopt = 0;

daystr = getTwoDigitNumber(day);

load([animDir, animal, 'tetinfo.mat']);
load([animDir, animal, 'cellinfo.mat']);
load([animDir, animal, 'spikes', daystr,'.mat']); % get spikes



epochs = find(~cellfun(@isempty,odorTriggers{1,day}));

for ep = 1:length(epochs)
    betalist = [];
    epoch = epochs(ep);
    %epstr = getTwoDigitNumber(epoch);

    trigs = odorTriggers{day}{epoch}.allTriggers;
    
    exclude_list = [];
    %[ctxidx, hpidx] = wt_findcellindex(animDir, animal, day, epoch, exclude_list); %(tet, cell) %may need to edit this to take more cells
    
    cellfilter = ['isequal($area,''CA1'') && strcmp($type, ''pyr'') && ($numspikes > 100)' ];
    hpidx = evaluatefilter(cellinfo{day}{epoch},cellfilter);
    
    cellfilter = ['isequal($area,''PFC'') && strcmp($type, ''pyr'') && ($numspikes > 100)' ];
    ctxidx = evaluatefilter(cellinfo{day}{epoch},cellfilter);
    
    hpnum = size(hpidx(:,1));
    ctxnum = size(ctxidx(:,1));


    betalist(:,1) = trigs;
    betalist(:,2) = trigs + 2;

%-----load beta eeg using CA1Ref tetrode-----%

    for ttt = 1:length(tetinfo{day}{epoch})
        if ~isempty(tetinfo{day}{epoch}{ttt})
            if isfield(tetinfo{day}{epoch}{ttt},'descrip')
               if strcmp(tetinfo{day}{epoch}{ttt}.descrip,'hpcRef')
                   reftet = ttt;
               end
            end
        end
    end
    %reftet = 10;               
    tmpflist1 = sprintf('%s%sbeta%02d-%02d-%02d.mat', [animDir,'EEG/'],animal, day, epoch, reftet);
    load(tmpflist1);
    t = geteegtimes(beta{day}{epoch}{reftet});
    tph = beta{day}{epoch}{reftet}.data(:,2);% beta phase
        
%%
%------HP cells----%
    for cell = hpnum:-1:1
        cind = hpidx(cell,:);
        s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
        %t = geteegtimes(beta{day}{epoch}{reftet});
        sph = tph(lookup(s, t));  
        goodspikes = isExcluded(s, betalist);
        numgoodspikes = sum(goodspikes);
        sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
        if length(sph)>1
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
            prayl=0;
            zrayl=0;
            alpha=0;
            pdf=0;
        end
        beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.index = cind;
        beta_phaselockCA1{day}{epoch}{cind(1)}{cind(2)}.tetindex = reftet;
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

        if (figopt == 1)
            %CS - removed plotting code to clean it up. Run wenbo's
            %original code if plotting individual cells is necessary.
        end
        

    end


    
%%
    %------CTX cells----%

    for cell = ctxnum:-1:1
        cind = ctxidx(cell,:);
        s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
        t = geteegtimes(beta{day}{epoch}{reftet});
        sph = tph(lookup(s, t));  
        goodspikes = isExcluded(s, betalist);
        numgoodspikes = sum(goodspikes);
        sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
        if length(sph)>1
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
            prayl=0;
            zrayl=0;
            alpha=0;
            pdf=0;
        end
        beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.index = cind;
        beta_phaselockPFC{day}{epoch}{cind(1)}{cind(2)}.tetindex = reftet;
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


    end
end

%%
if savedata
    save([animDir,animal,'betaphaselock_PFC',daystr,'.mat'],'beta_phaselockPFC');

    save([animDir,animal,'betaphaselock_CA1',daystr,'.mat'],'beta_phaselockCA1');
end



