%find tetrode with largest difference between baseline and after-trigger
%beta, use for betatet for that epoch. 
function cs_findBetaTets(animals,trialtype)
topDir = cs_setPaths();


regions = {'CA1','PFC','OB'};
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS41'};
%win = [0.5 2]; %2 seconds after trig, 1 second before trig for baseline

for a = 1:length(animals)
    animal = animals{a};
    
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cd(dataDir);
    
    %trigfiles = dir([dataDir,animal,'nosepokeWindow*']);
    %trigfiles = {trigfiles.name};
    
    load([dataDir,animal,'tetinfo.mat']);
    
    dayepochmatrix = cs_getRunEpochs(dataDir, animal, trialtype);
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
                refstr = 'OB';
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
                tets = evaluatefilter(tetinfo{day}{epoch}, tetfilter);
                
                zscores = zeros(length(tets),1);
                for t = 1:length(tets)
                    tet = tets(t);
                    tetstr = getTwoDigitNumber(tet);
                    
                    if isfield(tetinfo{day}{epoch}{tet}, 'descrip2')
                        tetinfo{day}{epoch}{tet} = rmfield(tetinfo{day}{epoch}{tet},'descrip2');
                    end
                    
                    load([dataDir,'EEG\',animal,'beta',daystr,'-',epstr,'-',tetstr,'.mat'])
                    
                    betamag = double(beta{day}{epoch}{tet}.data(:,3));
                    time = (beta{day}{epoch}{tet}.starttime:1/beta{day}{epoch}{tet}.samprate:beta{day}{epoch}{tet}.endtime)';
                    
                    %baseline is mean beta across epoch
                    baseline = mean(betamag);
                    stdev = std(betamag);
                    
                    %find all beta periods during nosepoke, average across
                    %all times
                    betaaftertrig = mean(betamag(any((time' > trigstart) & (time' <= (trigend)))));
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
                betatet = tets(ind);
                
                %if no tets had increase in beta, use 
                %tet with highest zscore as betatet
%                     if val <= 0
%                         reftetfilter = ['(isequal($descrip,''',refstr,'''))'];
%                         %reftetfilter = ['(isequal($area,''',refstr,'''))'];
%                         reftet = evaluatefilter(tetinfo{day}{epoch}, reftetfilter); 
%                         if isempty(reftet)
%                             try
%                                 reftetfilter = ['(isequal($area,''',refstr,'''))'];
%                                 reftet = evaluatefilter(tetinfo{day}{epoch}, reftetfilter); 
%                             catch
%                                 error('No tetrode found');
%                             end
%                         end
%                         tetstr = getTwoDigitNumber(reftet);
%                         tetinfo{day}{epoch}{reftet}.descrip2 = 'betatet';
%                         disp(['reftet ', tetstr, ' is ',region,' betatet for ',animal, ' day ', daystr,' epoch ',epstr ]);
%     %                   
%                     else
                        tetinfo{day}{epoch}{betatet}.descrip2 = 'betatet';
                        tetstr = getTwoDigitNumber(betatet);
                        disp(['tet ', tetstr, ' is ',region,' betatet for ',animal, ' day ', daystr,' epoch ',epstr ]);                  
%                     end
                
            end
            end

        end
    end
    save([dataDir,animal,'tetinfo.mat'],'tetinfo');

end