
%cs_phaseLocking_all
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
eegregions = {'CA1','PFC','OB'};
freqs = {'beta','resp'};
trialtypes = {'correct','incorrect'};

%% get NP cells
load([topDir,'AnalysesAcrossAnimals\npCells_CA1']);
cells_CA1 = activeCells;
load([topDir,'AnalysesAcrossAnimals\npCells_PFC']);
cells_PFC = activeCells; clear npCells

rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);


for a = 1:length(animals)
    animal = animals{a};
    animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
    
    %% get animal data
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    
    
    spikes = loaddatastruct(animDir, animal, 'spikes');
    days = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(days(:,1));
    
    for f = 1:length(freqs)
        freq = freqs{f};
        for crg = 1:length(cellregions)
            cellregion = cellregions{crg};
            
            for er = 1:length(eegregions)
                eegregion = eegregions{er};
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

                            [cl,cr,il,ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            wins_incorr = nosepokeWindow{day}{epoch}(sort([il;ir]),:);
                            if strcmp(trialtype,'correct')
                                list = nosepokeWindow{day}{epoch}(sort([cl;cr]),:);
                            elseif strcmp(trialtype,'incorrect')
                                list = nosepokeWindow{day}{epoch}(sort([il;ir]),:);
                            end
                            
                            %list = [list-0.5, list+0.5];
                            timelist = [timelist;list];
                            
                            numincorrtrials = numincorrtrials + size(wins_incorr,1);
                            %umincorrtrials = numincorrtrials + size(wins_incorr,1);
                            
                            tet = cs_getMostCellsTet(animal,day,epoch,eegregion);
                            tetstr = getTwoDigitNumber(tet);
                            % here I insert a way to reconstitute the
                            % filtered data
                            eegfile=[animDir, 'EEG\', animal, freq,daystr,'-',epstr,'-',tetstr,'.mat'];
                            if ~exist(eegfile,'file')
                                try
                                    filterFile=sprintf('C:\\Users\\Jadhavlab\\Documents\\gitRepos\\LFP-analysis\\JadhavEEGFilter\\%sfilter.mat',...
                                        freqs{f});
                                catch
                                    fprintf('cant load the eeg file');
                                    return
                                end
                                [filtered]=jhb_LFPtetprocess(animDir,animal,day,epoch,tet,'f',filterFile,'band',freqs{f});
                            end
                            load([animDir, 'EEG\', animal, freq,daystr,'-',epstr,'-',tetstr,'.mat'],freq);
                            
                            eval(['lfp = ',freq,';']);
                            
                            t = geteegtimes(lfp{day}{epoch}{tet});
                            lfptime = [lfptime;t'];
                            p = lfp{day}{epoch}{tet}.data(:,2);% beta phase
                            phase = [phase;p];
                        end
                        
                        %get cells (only task Responsive cells!!!!
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
                            
                            %calc stats
                            out = cs_calcPhaseLocking(sph);
                            
                            %store data in first epoch
                            phaselock{day}{1}{cell(1)}{cell(2)} = out;
                            phaselock{day}{1}{cell(1)}{cell(2)}.ind = cell;
                            phaselock{day}{1}{cell(1)}{cell(2)}.Nspikes = Nspikes;
                            
                            %do bootstrap for correct trials
                            if strcmp(trialtype,'correct')
                                [k_dist, z_dist, mvl_dist] = cs_phaselockBootstrap_v2(allspikes, lfptime, phase, timelist, numincorrtrials, 500);
                                phaselock{day}{1}{cell(1)}{cell(2)}.kappa_dist = k_dist;
                                phaselock{day}{1}{cell(1)}{cell(2)}.zrayl_dist = z_dist;
                                phaselock{day}{1}{cell(1)}{cell(2)}.mvl_dist = mvl_dist;

                            end
                            
                            
                        end
                        
                    end
                    
                    
                    
                    if strcmp(trialtype,'correct')
                        plDir = [animDir,'PhaseLocking\'];
                        if ~exist(plDir,'dir')
                            mkdir(plDir);
                        end
                    elseif strcmp(trialtype,'incorrect')
                        plDir = [animDir, 'PhaseLocking\Incorrect\'];
                         if ~exist(plDir,'dir')
                            mkdir(plDir);
                        end
                    end
                    
                    %if phaselocking variable does not exist, then there
                    %were no nosepoke cells on that day. in this case,
                    %delete any files that exist in the directory with
                    %matching name
                    if ~exist('phaselock','var')
                        oldfiles = dir([plDir, animal,'phaselock_',freq,trialstr,'_',cellregion,'-',eegregion,'_all.mat']);
                        if ~isempty(oldfiles)
                            delete (oldfiles.name)
                        end
                    else
                        %otherwise, save all days together
                        varfetch = cellfetch(phaselock,'Nspikes');
                        
                        save([plDir,animal,'phaselock_',freq,trialstr,'_',cellregion,'-',eegregion,'_all.mat'],'phaselock');
                        clear phaselock
                    end
                end
            end
        end
    end
end
