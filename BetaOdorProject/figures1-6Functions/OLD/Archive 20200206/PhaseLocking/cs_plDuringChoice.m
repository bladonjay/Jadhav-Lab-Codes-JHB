%cs_plDuringChoice
clear
[topDir, figDir] = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
lfpregions = {'CA1','PFC','OB'};
cellregions = {'CA1','PFC'};
freq = 'beta';
trialtype = 'correct';
win = [.4 0]; %time to look before and after choice

savedata = 0;

for cr = 1:length(cellregions)
    
    cellregion = cellregions{cr};
    for lfp = 1:length(lfpregions)
        plcells = [];
        lfpregion = lfpregions{lfp};
        for a = 1:length(animals)
            animal = animals{a};
            
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
            dayepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
            days = unique(dayepochs(:,1));
            
            for d = 1:length(days)
                day = days(d);
                daystr = getTwoDigitNumber(day);
                disp(['Doing ', animal,' Day ', daystr]);
                
                %load data
                npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
                tetinfo = loaddatastruct(animDir, animal,'tetinfo');
                cellinfo = loaddatastruct(animDir, animal,'cellinfo');
                spikes = loaddatastruct(animDir, animal, 'spikes');
                
                epochs = dayepochs(dayepochs(:,1) == day,2);
                
                for epoch = epochs'
                    timewins = [npWins{day}{epoch}(:,2)-win(1) npWins{day}{epoch}(:,2)+win(2)];
                    epstr = getTwoDigitNumber(epoch);
                    
                    tetfilter = ['((strcmp($area, ''',lfpregion,''')) && (strcmp($descrip2, ''betatet'')))'];
                    betatet = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                    if isempty(betatet)
                        switch lfpregion
                            case 'CA1'
                                ref = 'hpcRef';
                            case 'PFC'
                                ref = 'pfcRef';
                        end
                        tetfilter = ['strcmp($descrip,''',ref,''') && strcmp($descrip2,''betatet'')'];
                        betatet = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                    end
                    
                    tetstr = getTwoDigitNumber(betatet);
                    
                    load([animDir, 'EEG/', animal, 'beta',daystr,'-',epstr,'-',tetstr,'.mat'],'beta');
                    t = geteegtimes(beta{day}{epoch}{betatet});
                    betaphase = beta{day}{epoch}{betatet}.data(:,2);% beta phase
                    
                    load([topDir,'AnalysesAcrossAnimals\npCells_',cellregion,'.mat']);
                    idx = npCells(ismember(npCells(:,[1,2]),[a, day], 'rows'),[3 4]);
                    cellnum = size(idx,1);
                    
                    for cell = cellnum:-1:1
                        cind = idx(cell,:);
                        if ~isempty(spikes{day}{epoch}{cind(1)}{cind(2)}.data)
                            s = spikes{day}{epoch}{cind(1)}{cind(2)}.data(:,1);
                            %t = geteegtimes(beta{day}{epoch}{reftet});
                            sph = betaphase(lookup(s, t));
                            goodspikes = isExcluded(s, timewins);
                            numgoodspikes = sum(goodspikes);
                            sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                            if length(sph)> 1
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
                                
                                if prayl < 0.05
                                    plcells = [plcells; a day epoch cind(1) cind(2)];
                                end
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
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.index = cind;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.tetindex = betatet;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.sph = sph;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.Nspikes = numgoodspikes;
                            % output stats also, although you will have to redo this after combining epochochs
                            % Rayleigh test
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.stats = stats;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.modln = m;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.phdeg = phdeg;
                            % Von Mises Fit
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.kappa = kappa;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.betahat = betahat;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.betahat_deg = betahat_deg;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.prayl = prayl;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.zrayl = zrayl;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.alpha = alpha;
                            betaphaselock{day}{epoch}{cind(1)}{cind(2)}.vmpdf = pdf;
                            
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
                        
                        
                    end
                end
                if savedata == 1
                save([animDir,'PhaseLocking\NPChoice\',animal,'betaphaselock_',cellregion,'-',lfpregion,'_',daystr,'.mat'],'betaphaselock');
                end
                clear betaphaselock
            end
        end
        save([topDir,'AnalysesAcrossAnimals\PhaseLocking\NPChoice\',cellregion,'-',lfpregion,'.mat'],'plcells');

    end
end

cs_selectivityPhaseLockOverlap_v4
