
%cs_phaseLocking
%created 3/16/2020 to streamline phaselocking pipeline

%phaselock files will contain simple "phaselock" variables. this file will
%contain data across all days for one animal, for one set of region pairs
%cell regions done independently
%calls a function that calculates phase locking for a single cell
%incorrect trials still saved independently
clear
%% Params
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'OB'};
freqs = {'beta','resp'};
trialtypes = {'correct','incorrect'};

%% get NP cells
load([topDir,'AnalysesAcrossAnimals\npInt_CA1']);
cells_CA1 = npInt;
load([topDir,'AnalysesAcrossAnimals\npInt_PFC']);
cells_PFC = npInt; clear npInt


for a = 1:length(animals)
    animal = animals{a};
    animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
    
    %% get animal data
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    
    load([animDir,animal,'highBeta']);
    windows = highBeta;
    load([animDir,animal,'highBeta_incorrect']);
    wins_incorr = highBeta; clear highBeta
    
    spikes = loaddatastruct(animDir, animal, 'spikes');
    
    days = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(days(:,1));
    
    for f = 1:length(freqs)
        freq = freqs{f};
        for cr = 1:length(cellregions)
            cellregion = cellregions{cr};
            
            for er = 1:length(eegregions)
                eegregion = eegregions{er};
                
                %                 if strcmp(cellregion,eegregion)
                %                     continue
                %                 end
                disp(['Doing ',animal, ' ',freq,' ',cellregion,'-',eegregion])
                for tt = 1:length(trialtypes)
                    trialtype = trialtypes{tt};
                    
                    switch trialtype
                        case 'correct'
                            trialstr = '';
                            
                        case 'incorrect'
                            trialstr = '_incorrect';
                    end
                    
                    for day = days'
                        
                        daystr = getTwoDigitNumber(day);
                        epochs = cs_getRunEpochs(animDir, animal, 'odorplace',day);
                        epochs = epochs(:,2);
                        
                        %get windows, lfp phases, and number of incorr trials
                        timelist = [];
                        lfptime = [];
                        phase = [];
                        numincorrtrials = 0;
                        for ep = 1:length(epochs)
                            epoch = epochs(ep);
                            epstr = getTwoDigitNumber(epoch);
                            
%                             if strcmp(trialtype,'correct')
%                                 list = windows{day}{epoch}.OB;
%                             elseif strcmp(trialtype,'incorrect')
%                                 list = wins_incorr{day}{epoch}.OB;
%                             end
                            %
                            [cl,cr,il,ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            wins_incorr = nosepokeWindow{day}{epoch}(sort([il;ir]),:);
                            if strcmp(trialtype,'correct')
                                list = nosepokeWindow{day}{epoch}(sort([cl;cr]),:);
                            elseif strcmp(trialtype,'incorrect')
                                list = nosepokeWindow{day}{epoch}(sort([il;ir]),:);
                            end
                            
                           % list = [list-1, list];
                            timelist = [timelist;list];
                            
                            numincorrtrials = numincorrtrials + size(wins_incorr,1);
                            %umincorrtrials = numincorrtrials + size(wins_incorr,1);
                            
                            tet = cs_getMostCellsTet(animal,day,epoch,eegregion);
                            tetstr = getTwoDigitNumber(tet);
                            
                            load([animDir, 'EEG/', animal, freq,daystr,'-',epstr,'-',tetstr,'.mat'],freq);
                            eval(['lfp = ',freq,';']);
                            
                            t = geteegtimes(lfp{day}{epoch}{tet});
                            lfptime = [lfptime;t'];
                            p = lfp{day}{epoch}{tet}.data(:,2);% beta phase
                            phase = [phase;p];
                        end
                        
                        %get cells
                        eval(['cells = cells_',cellregion,'(ismember(cells_',cellregion,'(:,[1,2]),[a, day], ''rows''),[3 4]);']);
                        
                        %do each cell separately
                        for c = 1:size(cells)
                            cell = cells(c,:);
                            allspikes = [];
                            
                            %combine all spikes over epochs
                            for epoch = epochs'
                                if ~isempty(spikes{day}{epoch}{cell(1)}{cell(2)}.data)
                                    s = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                                    allspikes = [allspikes;s];
                                end
                            end
                            
                            
                            
                            %get spike phase
                            sph_tmp = phase(lookup(allspikes, lfptime));
                            goodspikes = isExcluded(allspikes, timelist);
                            Nspikes = sum(goodspikes);
                            sph = double(sph_tmp(find(goodspikes))) / 10000;  % If no spikes, this will be empty
                            
                            
%                             phase_tmp = double(phase)/10000;
%                             time_tmp = lfptime(isExcluded(lfptime,[8287 8290]));                          sp = allspikes(goodspikes);
%                             phase_tmp = phase_tmp(isExcluded(lfptime,[8287 8290]));
%                             
%                             plot(time_tmp, phase_tmp);
%                             bad = find(sph == 1.5706);
%                             badtime = sp(bad);
%                             
%                             plotphase = double(phase)/10000;
%                             
                            %remove noise
                                %occasionally, long trains of consecutive
                                %"spikes" will have exactly the same phase.
                                %This is likely noise? 
                            
                            %calc stats
                            out = cs_calcPhaseLocking(sph);
                            
%                             %check result
%                             if ~isempty(sph)
%                                 binedges = -pi:(2*pi/20):pi;
%                                 count = histcounts(sph, binedges);
%                                 pct = (count./sum(count))*100;
%                                 newcount = [pct,pct];
%                                 
%                                 newbinedges = [binedges - pi, binedges(2:end) + pi];
%                                 newbins = newbinedges(2:end)-((binedges(2)-binedges(1))/2);
%                                 
%                                 figure, hold on
%                                 bar(newbins, newcount,1);
%                                 plot([0 0], [0 max(newcount+2)], 'k--');
%                                 axis([-2*pi, 2*pi, 0 max(newcount+2)])
%                                 xticks([-2*pi -pi 0 pi 2*pi])
%                                 xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
%                                 xlabel('Phase')
%                                 ylabel('% of spikes')
%                                 
%                                 
%                             end
%                             close all
                            
                            %store data in first epoch
                            phaselock{day}{1}{cell(1)}{cell(2)} = out;
                            phaselock{day}{1}{cell(1)}{cell(2)}.ind = cell;
                            phaselock{day}{1}{cell(1)}{cell(2)}.Nspikes = Nspikes;
                            
                            %do bootstrap for correct trials
                            if strcmp(trialtype,'correct')
                                [k_dist, z_dist] = cs_phaselockBootstrap(allspikes, lfptime, phase, timelist, numincorrtrials, 500);
                                phaselock{day}{1}{cell(1)}{cell(2)}.kappa_dist = k_dist;
                                phaselock{day}{1}{cell(1)}{cell(2)}.zrayl_dist = z_dist;
                            end
                            
                            
                        end
                        
                    end
                    
                    
                    
                    if strcmp(trialtype,'correct')
                        plDir = [animDir,'PhaseLocking\Interneurons\'];
                    elseif strcmp(trialtype,'incorrect')
                        plDir = [animDir, 'PhaseLocking\Interneurons\Incorrect\'];
                    end
                    
                    %if phaselocking variable does not exist, then there
                    %were no nosepoke cells on that day. in this case,
                    %delete any files that exist in the directory with
                    %matching name
                    if ~exist('phaselock','var')
                        oldfiles = dir([plDir, animal,'phaselock_',freq,trialstr,'_',cellregion,'-',eegregion,'.mat']);
                        if ~isempty(oldfiles)
                            delete (oldfiles.name)
                        end
                    else
                        %otherwise, save all days together
                        varfetch = cellfetch(phaselock,'Nspikes');
                        
                        save([plDir,animal,'phaselock_',freq,trialstr,'_',cellregion,'-',eegregion,'.mat'],'phaselock');
                        clear phaselock
                    end
                end
            end
        end
    end
end
cs_listPhaseLockedINs