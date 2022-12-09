

function cs_getNosepokeWindow(animals)

% Gives cell array with epochs for each day, and a structure for each epoch that
% includes all odor times, left odor times, right odor times, correct trial
% odor times, and incorrect trial odor times. (odorTriggers files)

% Also creates nosepokeWindow files, which indicate when rat pulled out of
% nosepoke on each trial

%Also saves binary performance- 1s for correct trials, 0s for incorrect
%trials, for use with state space learning curve

[topDir,~] = cs_setPaths();

for a = 1:length(animals)
    animal = animals{a};
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    cd(dataDir)
    
    epochs = cs_getRunEpochs(dataDir, animal,'odorplace');
%     taskfiles = dir([animal,'task*']);
%     
%     epochs = [];
%     for tk = 1:length(taskfiles)
%         load(taskfiles(tk).name);
%         
%         epochfilter = '(strcmp($environment, ''odorplace''))';
%         eps = evaluatefilter(task,epochfilter);
%         epochs = [epochs;eps];
%     end
    
    days = unique(epochs(:,1));
    
    for d = 1:length(days)
        day = days(d);
        
        if (day<10)
            daystring = ['0',num2str(day)];
        else
            daystring = num2str(day);
        end
        
        runEps = epochs(epochs(:,1) == day, 2);
        numEpochs = length(runEps);
        
        load([animal,'DIO', daystring, '.mat']);
        % if the dio was accidentally saved as a 1 epoch matrix
        if length(dio{day})==1 && numEpochs==1 && runEps(1)>1     
            dio{1,day}{runEps(1)}=dio{1,day}{1}; % have to put the data in the correct epoch
        end
        
        
        for e = 1:length(runEps)
            ep = runEps(e); % right now there is only one epoch
            
            if (ep <10)
                epstring = ['0',num2str(ep)];
            else
                epstring = num2str(ep);
            end
            
            %%  Get DIO times and states for reward wells, odor solenoids, nosepoke, and buzzer
            % this is old, the dios were collapsed across epochs...
            %
            % code hangs here for CS41 because dio is not split by epochs
            % she collapsed epochs for this animal so there is a single
            % cell in that array e.g. size DIO is {1day} {1epochonly}
            % {32dios} instead of {1firstday}{3epochs}{32dios}
            %
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
            
            
            Dout6times = dio{1,day}{1,ep}{1,22}.time(logical(dio{1,day}{1,ep}{1,22}.state));
            Dout6hits = ones([length(Dout6times),1]) +5;
            leftsolnans = nan([length(Dout6hits),1]);
            LeftSolenoid = [Dout6times, leftsolnans, Dout6hits, leftsolnans, leftsolnans];
            clear('Dout6times','Dout6hits','leftsolnans');
            
            
            Dout7times = dio{1,day}{1,ep}{1,23}.time(logical(dio{1,day}{1,ep}{1,23}.state));
            Dout7hits = ones([length(Dout7times),1]) +6;
            rightsolnans = nan([length(Dout7hits),1]);
            RightSolenoid = [Dout7times, rightsolnans, Dout7hits, rightsolnans, rightsolnans];
            clear('Dout7times','Dout7hits','rightsolnans');
            
            
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
            
            triggers.allTriggers = allTriggers;
            
            triggers.leftTriggers = leftTriggers;
            triggers.rightTriggers = rightTriggers;
            
            triggers.correctTriggers = correctTriggers;
            triggers.incorrectTriggers = incorrectTriggers;
            
            %append existing file if it exists
            try
                odorTriggers = loaddatastruct(dataDir, animal,'odorTriggers',day);
            catch
                disp('');
            end
            
            odorTriggers{1,day}{1,ep} = triggers;
            BinaryPerf{1,day}{1,ep} = binaryPerf;
            nosepokeWindow{1,day}{1,ep} = windows;
        end
        save([dataDir,animal, 'odorTriggers', daystring],'odorTriggers');
        save([dataDir,animal, 'nosepokeWindow', daystring],'nosepokeWindow');
        
        clear odorTriggers
        clear nosepokeWindow
        
        if ~exist([dataDir,'BinaryPerf'])
            mkdir([dataDir,'BinaryPerf'])
        end
        cd([dataDir,'BinaryPerf']);
        save([animal,'BinaryPerf',daystring],'BinaryPerf');
        
        clear BinaryPerf
        
        cd(dataDir)
        
    end
    
end
