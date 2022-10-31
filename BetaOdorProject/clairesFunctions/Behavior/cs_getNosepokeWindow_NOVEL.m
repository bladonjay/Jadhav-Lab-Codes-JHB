

function cs_getNosepokeWindow_NOVEL(animals)

% Gives cell array with epochs for each day, and a structure for each epoch that
% includes all odor times, left odor times, right odor times, correct trial
% odor times, and incorrect trial odor times.

% Also creates nosepokeWindow files, which indicate when rat pulled out of
% nosepoke on each trial

%Also saves binary performance- 1s for correct trials, 0s for incorrect
%trials, for use with state space learning curve

%CS- edited 9/10/19 to get trig times for all task environments including
%novel

[topDir,~] = cs_setPaths();
envs = {'odorplace','novelodor','novelodor2','noodor'};

for a = 1:length(animals)
    animal = animals{a};
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    %cd(dataDir)
    
    taskfiles = dir([dataDir,animal,'task*']);
    
    %epochs = [];
    %     for tk = 1:length(taskfiles)
    %         load(taskfiles(tk).name);
    %
    %         epochfilter = '(strcmp($environment, ''novel1''))';
    %         eps = evaluatefilter(task,epochfilter);
    %         epochs = [epochs;eps];
    %     end
    
    numdays = length(taskfiles);
    
    for day = 1:numdays
        %day = d;
        daystring = getTwoDigitNumber(day);
        
        if strcmp(animal,'CS41') && day == 4
            continue
        end
        
        load(taskfiles(day).name);
        
        for en = 1:length(envs)
            env = envs{en};
            runEps = cs_getRunEpochs(dataDir, animal, env,day);
            runEps = runEps(:,2);
            %epochfilter = '(strcmp($environment, ''novel1''))';
            %         eps = evaluatefilter(task,epochfilter);
            %         epochs = [epochs;eps];
            
            load([dataDir, animal, 'DIO', daystring]);
            %numEpochs = length(runEps);
            
            for e = 1:length(runEps)
                ep = runEps(e);
                %epstring = getTwoDigitNumber(ep);
                
                %%  Get DIO times and states for reward wells, odor solenoids, nosepoke, and buzzer
                
                Din1times = dio{1,day}{1,ep}{1,1}.time(logical(dio{1,day}{1,ep}{1,1}.state));
                Din1hits = ones([length(Din1times),1]);
                leftwellnans = nan([length(Din1hits),3]);
                LeftWell = [Din1times, leftwellnans, Din1hits];
                clear('Din1times','Din1hits','leftwellnans');
                
                Din2times = dio{1,day}{1,ep}{1,2}.time(logical(dio{1,day}{1,ep}{1,2}.state));
                Din2hits = ones([length(Din2times),1])+1;
                rightwellnans = nan([length(Din2hits),3]);
                RightWell = [Din2times, rightwellnans, Din2hits];
                clear('Din2times','Din2hits','rightwellnans');
                
                Reward = sortrows([LeftWell; RightWell]);
                
                if strcmp(env, 'odorplace') || (strcmp(animal, 'CS42') && strcmp(env,'novelodor')) || (strcmp(animal,'CS41') && strcmp(env,'novelodor2')) ||...
                        (strcmp(animal, 'CS44') && strcmp(env,'novelodor')) || (strcmp(animal,'CS44') && strcmp(env,'novelodor2')) || ... 
                        (strcmp(animal,'CS44') && strcmp(env,'noodor'))
                    %
                    Dout6times = dio{1,day}{1,ep}{1,22}.time(logical(dio{1,day}{1,ep}{1,22}.state));
                    Dout6hits = ones([length(Dout6times),1]) + 5;
                    leftsolnans = nan([length(Dout6hits),1]);
                    LeftSolenoid = [Dout6times, leftsolnans, Dout6hits, leftsolnans, leftsolnans];
                    clear('Dout6times','Dout6hits','leftsolnans');
                    
                    Dout7times = dio{1,day}{1,ep}{1,23}.time(logical(dio{1,day}{1,ep}{1,23}.state));
                    Dout7hits = ones([length(Dout7times),1]) + 6;
                    rightsolnans = nan([length(Dout7hits),1]);
                    RightSolenoid = [Dout7times, rightsolnans, Dout7hits, rightsolnans, rightsolnans];
                    clear('Dout7times','Dout7hits','rightsolnans');
                else
                    
                    
                    Dout9times = dio{1,day}{1,ep}{1,25}.time(logical(dio{1,day}{1,ep}{1,25}.state));
                    Dout9hits = ones([length(Dout9times),1]) +5;
                    leftsolnans = nan([length(Dout9hits),1]);
                    LeftSolenoid = [Dout9times, leftsolnans, Dout9hits, leftsolnans, leftsolnans];
                    clear('Dout9times','Dout9hits','leftsolnans');
                    
                    
                    Dout10times = dio{1,day}{1,ep}{1,26}.time(logical(dio{1,day}{1,ep}{1,26}.state));
                    Dout10hits = ones([length(Dout10times),1]) +6;
                    rightsolnans = nan([length(Dout10hits),1]);
                    RightSolenoid = [Dout10times, rightsolnans, Dout10hits, rightsolnans, rightsolnans];
                    clear('Dout10times','Dout10hits','rightsolnans');
                end
                
                
                Odors = sortrows([LeftSolenoid; RightSolenoid]);
                
                
                %%%% Buzzer is now on Dout5 - use for all animals after CS31 %%%%
                
                if strcmp(animal, 'CS31') || strcmp(animal, 'CS42')
                    Dout9times = dio{1,day}{1,ep}{1,25}.time(logical(dio{1,day}{1,ep}{1,25}.state));
                    Dout9hits = ones([length(Dout9times),1]);
                    buzzernans = nan([length(Dout9hits),1]);
                    Buzzer = [Dout9times, buzzernans, buzzernans, Dout9hits, buzzernans];
                    clear('Dout9times','Dout9hits','buzzernans');
                    
                else
                    
                    Dout5times = dio{1,day}{1,ep}{1,21}.time(logical(dio{1,day}{1,ep}{1,21}.state));
                    Dout5hits = ones([length(Dout5times),1]);
                    buzzernans = nan([length(Dout5hits),1]);
                    Buzzer = [Dout5times, buzzernans, buzzernans, Dout5hits, buzzernans];
                    clear('Dout5times','Dout5hits','buzzernans');
                    
                end
                
                Din5times = dio{1,day}{1,ep}{1,5}.time;
                Din5states = double(dio{1,day}{1,ep}{1,5}.state);
                npNans = nan([length(Din5states),3]);
                NP = [Din5times, Din5states, npNans];
                clear('Din5times','Din5hits','NPnans');
                
                
                %%
                %create matrix with np state, which odor was dispensed, when buzzer
                %went off, and which reward well was trigered.
                
                attempts = sortrows([NP; Odors; Buzzer; Reward]);
                
                %for each buzzer, find the nosepoke time and the odor ID that was immediately before it.
                % use those np times for triggers, and also NP windows
                % separate into R/L odors, then find R/L trials and
                % correct/incorrect from there using existing method
                allTriggers = []; windows = [];
                rightTriggers = [];
                leftTriggers = [];
                correctTriggers = [];
                incorrectTriggers = [];
                binaryPerf = [];
                
                buzzerInds = find(~isnan(attempts(:,4)));
                for b = 1:length(buzzerInds)
                    
                    Bind = buzzerInds(b);
                    
                    %if a reward well is triggered before another buzzer, keep this
                    %buzzer time. Eliminates times when rat triggered buzzer twice
                    %before going to reward.
                    nextbuzzer = find(~isnan(attempts(Bind+1:end,4)),1,'first')+Bind;
                    if isempty(nextbuzzer)
                        nextbuzzer = size(attempts,1)+1;
                    end
                    nextreward = find(~isnan(attempts(Bind+1:end,5)),1,'first')+Bind;
                    
                    if nextreward < nextbuzzer
                        
                        %NPind = find(~isnan(attempts(1:Bind-1,2)),1,'last');
                        NPind = find(attempts(1:Bind-1,2) == 1,1,'last');
                        trigger = attempts(NPind,1);
                        NPoffind = find(attempts(NPind+1:end,2) == 0,1,'first')+NPind;
                        NPoff = attempts(NPoffind,1);
                        
                        window = [trigger, NPoff];
                        
                        %Statescript code for behaviors had errors in it for CS31,
                        %causing buzzer to go off sometimes if rat nosepoked twice
                        %in quick succession. In this case, take the nosepoke just
                        %before the last odor was dispensed.
                        if window(2)-window(1) < 0.5
                            odorind = find(~isnan(attempts(1:Bind-1,3)),1,'last');
                            NPind = find(attempts(1:odorind-1,2) == 1,1,'last');
                            trigger = attempts(NPind,1);
                            %                     NPoffind = find(attempts(NPind+1:end,2) == 0,1,'first')+NPind;
                            %                     NPoff = attempts(NPoffind,1);
                            
                            window = [trigger, NPoff];
                        end
                        
                        
                        %sometimes the buzzer still went off. not sure why.
                        %just exclude these trials (should be <10)
                        if window(2)-window(1) <0.5
                            trigger = [];
                            window = [];
                        end
                        
                        allTriggers = [allTriggers; trigger];
                        windows = [windows; window];
                        
                        
                        % Determine if it was a right or left trial, and if it was
                        % correct or incorrect
                        odorID = attempts(find(~isnan(attempts(1:Bind-1,3)),1,'last'),3);
                        rewardID = attempts(nextreward,5);
                        
                        if odorID == 6
                            leftTriggers = [leftTriggers; trigger];
                            
                            if rewardID == 1
                                correctTriggers = [correctTriggers; trigger];
                                binaryPerf = [binaryPerf;1];
                                
                            elseif rewardID == 2
                                incorrectTriggers = [incorrectTriggers; trigger];
                                binaryPerf = [binaryPerf;0];
                            else
                                warning(['check dio matrix, may be misaligned, index ', num2str(Bind)])
                            end
                            
                        elseif odorID == 7
                            rightTriggers = [rightTriggers; trigger];
                            
                            if rewardID == 2
                                correctTriggers = [correctTriggers; trigger];
                                binaryPerf = [binaryPerf;1];
                                
                            elseif rewardID == 1
                                incorrectTriggers = [incorrectTriggers; trigger];
                                binaryPerf = [binaryPerf;0];
                                
                            else
                                warning(['check dio matrix, may be misaligned, index ', num2str(Bind)])
                            end
                            
                        else
                            warning(['check dio matrix, may be misaligned, index ', num2str(Bind)])
                        end
                    end
                end
                
                %some sessions had both familiar and novel... for now, exclude
                %familiar trials
                allTriggers = sortrows([correctTriggers; incorrectTriggers]);
                triggers.allTriggers = allTriggers;
                windows = windows(ismember(windows(:,1), allTriggers),:);
                
                
                triggers.leftTriggers = leftTriggers;
                triggers.rightTriggers = rightTriggers;
                
                triggers.correctTriggers = correctTriggers;
                triggers.incorrectTriggers = incorrectTriggers;
                
                odorTriggers{1,day}{1,ep} = triggers;
                BinaryPerf{1,day}{1,ep} = binaryPerf;
                nosepokeWindow{1,day}{1,ep} = windows;
            end
        end
        save([dataDir,animal, 'odorTriggers', daystring],'odorTriggers');
        save([dataDir,animal, 'nosepokeWindow', daystring],'nosepokeWindow');
        
        clear odorTriggers
        clear nosepokeWindow
        
        if ~exist([dataDir,'BinaryPerf'])
            mkdir([dataDir,'BinaryPerf'])
        end
        %cd([dataDir,'BinaryPerf']);
        save([dataDir,'BinaryPerf\',animal,'BinaryPerf',daystring],'BinaryPerf');
        
        clear BinaryPerf
        
    end
    
end
