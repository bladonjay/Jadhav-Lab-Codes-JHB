%find OB tetrode with largest difference between baseline and after-trigger
%RR (6-12 Hz, use theta), use for rrtet for that epoch.

topDir = cs_setPaths();


regions = {'OB'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS39'};
%win = [0.5 2]; %2 seconds after trig, 1 second before trig for baseline

for a = 1:length(animals)
    animal = animals{a};
    
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cd(dataDir);
    
    %trigfiles = dir([dataDir,animal,'nosepokeWindow*']);
    %trigfiles = {trigfiles.name};
    
    load([dataDir,animal,'tetinfo.mat']);
    
    dayepochmatrix = cs_getRunEpochs(dataDir, animal, 'odorplace');
    days = unique(dayepochmatrix(:,1));
    
    %     taskfiles = dir([dataDir, animal, 'task*']);
    %     taskfiles = {taskfiles.name};
    
    for r = 1:length(regions)
        region = regions{r};
        
        switch region
            case 'CA1'
                refstr = 'hpcRef';
            case 'PFC'
                refstr = 'pfcRef';
            case 'OB'
                refstr = 'obRef';
        end
        
        
        for d = 1:length(days)
            day = days(d);
            
            runeps = dayepochmatrix(dayepochmatrix(:,1) == day, 2);
            
            if isempty(runeps)
                continue
            else
                %runeps = runeps(:,2);
                
                %day = length(task);
                %daytrigs = nosepokeWindow{d};
                %runeps = find(~cellfun(@isempty,nosepokeWindow{day}));
                daystr = getTwoDigitNumber(day);
                load([dataDir, animal, 'nosepokeWindow',daystr,'.mat'])
                
                for e = 1:length(runeps)
                    epoch = runeps(e);
                    epstr = getTwoDigitNumber(epoch);
                    trigstart = nosepokeWindow{day}{epoch}(:,1);
                    trigend = nosepokeWindow{day}{epoch}(:,2);
                    
                    tetfilter = ['(isequal($area, ''',region,'''))'];
                    tets = evaluatefilter(tetinfo{d}{epoch}, tetfilter);
                    
                    if length(tets) > 1
                        
                        zscores = zeros(length(tets),1);
                        for t = 1:length(tets)
                            tet = tets(t);
                            tetstr = getTwoDigitNumber(tet);
                            
                            if isfield(tetinfo{d}{epoch}{tet}, 'descrip3')
                                tetinfo{d}{epoch}{tet} = rmfield(tetinfo{d}{epoch}{tet},'descrip3');
                            end
                            
                            %Use theta for RR- it is about 6-12 Hz, but increases
                            %during NP window
                            load([dataDir,'EEG\',animal,'theta',daystr,'-',epstr,'-',tetstr,'.mat'])
                            
                            thetamag = double(theta{d}{epoch}{tet}.data(:,3));
                            time = (theta{d}{epoch}{tet}.starttime:1/theta{d}{epoch}{tet}.samprate:theta{d}{epoch}{tet}.endtime)';
                            
                            %baseline is mean beta across epoch - might change
                            %this to 500ms before trig
                            baselinetimes = [trigstart - 0.5, trigstart];
                            baselinemags = thetamag(isExcluded(time,baselinetimes));
                            baseline = mean(baselinemags);
                            stdev = std(baselinemags);
                            
                            %find all beta periods during nosepoke, average across
                            %all times
                            betaaftertrig = mean(thetamag(any((time' > trigstart) & (time' <= (trigend)))));
                            zscores(t) = (betaaftertrig - baseline) / stdev;
                            
                            %                     betameans = zeros(length(trigs),2);
                            %                     for tr = 1:length(trigs)
                            %                         trig = trigs(tr);
                            %                         betaaftertrig = betamag((time' > trig) & (time <= trig+win(2)));
                            %                         betabaseline = betamag((time > trig-win(1)) & (time <= trig));
                            %
                            %                         betaaftertrig = mean(betaaftertrig);
                            %                         betabaseline = mean(betabaseline);
                            %
                            %                         betameans(tr,:) = [betabaseline, betaaftertrig];
                            %                     end
                            
                            %                     betamean = mean(betameans);
                            %                     diffs(t,1) = (betamean(2) - betamean(1));
                            
                        end
                        
                        
                        
                        [val,ind] = max(zscores);
                        rrtet = tets(ind);
                        
                        %if no tets had increase in beta, use ref. otherwise, save
                        %tet with highest zscore as betatet
                        if val <= 0
                            reftetfilter = ['(isequal($descrip,''',refstr,'''))'];
                            reftet = evaluatefilter(tetinfo{d}{epoch}, reftetfilter);
                            tetstr = getTwoDigitNumber(reftet);
                            tetinfo{d}{epoch}{reftet}.descrip3 = 'rrtet';
                            disp(['reftet ', tetstr, ' is ',region,' rrtet for ',animal, ' day ', daystr,' epoch ',epstr ]);
                            %
                        else
                            tetinfo{d}{epoch}{rrtet}.descrip3 = 'rrtet';
                            tetstr = getTwoDigitNumber(rrtet);
                            disp(['tet ', tetstr, ' is ',region,' rrtet for ',animal, ' day ', daystr,' epoch ',epstr ]);
                        end
                        
                    else
                        rrtet = tets;
                        tetinfo{d}{epoch}{rrtet}.descrip3 = 'rrtet';
                        tetstr = getTwoDigitNumber(rrtet);
                        disp(['tet ', tetstr, ' is ',region,' rrtet for ',animal, ' day ', daystr,' epoch ',epstr ]);
                    end
                    
                end
            end
            
        end
    end
    save([dataDir,animal,'tetinfo.mat'],'tetinfo');
    
end