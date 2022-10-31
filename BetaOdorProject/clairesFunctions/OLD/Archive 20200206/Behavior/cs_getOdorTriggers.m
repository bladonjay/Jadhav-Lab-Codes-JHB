%% USE cs_getNosepokeWindow INSTEAD %%

function cs_getOdorTriggers(prefix)

% gives cell array with epochs for each day, and a structure for each epoch that
% includes all odor times, left odor times, right odor times, correct trial
% odor times, and incorrect trial odor times. 

%Enter dataDir WITH trailing slash

%(can also save binary performance)
[topDir,~] = cs_setPaths();

animDir = [topDir,prefix,'Expt\'];
dataDir =([animDir,prefix,'_direct\']);
cd(dataDir)

epochs = cs_getRunEpochs(dataDir, prefix, 'odorplace');
days = unique(epochs(:,1));

for d = 1:length(days)
    day = days(d);
    daystring = getTwoDigitNumber(day);
        
    load([prefix,'DIO', daystring, '.mat']);
    %numEpochs = length(dio{1,day});
    
    runEps = epochs(find(epochs(:,1) == day),2);
    
    for e = 1:length(runEps)
        epoch = runEps(e);
        
        epstring = getTwoDigitNumber(epoch);
        %if epoch <= numEpochs
         
        Din1times = dio{1,day}{1,epoch}{1,1}.time(find(dio{1,day}{1,epoch}{1,1}.state));  
        Din1hits = ones([length(Din1times),1]);
        leftwellzeros = zeros([length(Din1hits),1]);
        LeftWell = [Din1times, leftwellzeros, Din1hits, leftwellzeros];

        Din2times = dio{1,day}{1,epoch}{1,2}.time(find(dio{1,day}{1,epoch}{1,2}.state));
        Din2hits = ones([length(Din2times),1]);
        rightwellzeros = zeros([length(Din2hits),2]);
        RightWell = [Din2times, rightwellzeros, Din2hits];

        Dout6times = dio{1,day}{1,epoch}{1,22}.time(find(dio{1,day}{1,epoch}{1,22}.state));
        Dout6hits = ones([length(Dout6times),1]);
        leftsolzeros = zeros([length(Dout6hits),1]);
        LeftSolenoid = [Dout6times, Dout6hits, leftsolzeros];

        Dout7times = dio{1,day}{1,epoch}{1,23}.time(find(dio{1,day}{1,epoch}{1,23}.state));
        Dout7hits = ones([length(Dout7times),1]);
        rightsolzeros = zeros([length(Dout7hits),1]);
        RightSolenoid = [Dout7times, Dout7hits, rightsolzeros];

        if strcmp(prefix, 'CS31')
            Dout9times = dio{1,day}{1,epoch}{1,25}.time(find(dio{1,day}{1,epoch}{1,25}.state));
            Dout9hits = ones([length(Dout9times),1]);
            buzzerzeros = zeros([length(Dout9hits),1]);
            Buzzer = [Dout9times, buzzerzeros, Dout9hits];
            
        else
            
            % Buzzer is now on Dout5 - use for all animals after CS31
            
            Dout5times = dio{1,day}{1,epoch}{1,21}.time(find(dio{1,day}{1,epoch}{1,21}.state));
            Dout5hits = ones([length(Dout5times),1]);
            buzzerzeros = zeros([length(Dout5hits),1]);
            Buzzer = [Dout5times, buzzerzeros, Dout5hits];
        end
        
%%
        leftAttemps = sortrows([LeftSolenoid; Buzzer]);

        leftsolenoidInds = find(leftAttemps(:,2)); %finds all the times where the left solenoid was triggered
        leftRewardAvailInds = [];   
        for i = [1:length(leftsolenoidInds)] %each time left solenoid triggered, checks to see if buzzer went off immediately after. 
            currIndex = leftsolenoidInds(i);
            if currIndex ~= length(leftAttemps)
                buzzerCheck = leftAttemps(currIndex+1,3);
                    if buzzerCheck == 1
                    leftRewardAvailInds = [leftRewardAvailInds; currIndex]; %beginning of trials (after successful nose pokes), next reward well triggered is correct or incorrect trial
                    end
            end
        end

        leftTrials = leftAttemps(:,1:2);
        leftRewardAvail = leftTrials(leftRewardAvailInds,1:2);
        leftRewzeros = zeros([length(leftRewardAvailInds),2]);
        leftRewardAvail = [leftRewardAvail, leftRewzeros; LeftWell; RightWell];
        leftTrialData = sortrows(leftRewardAvail); %creates matrix with each time a successful nose poke occured, and when each reward well was triggered

        leftTrialStartInds = find(leftTrialData(:,2));
        leftCompleteTrialIndex = [];
        leftCorrectTrialIndex = [];
        leftIncorrectTrialIndex = [];
        for i = [1:length(leftTrialStartInds)]
            currIndex = leftTrialStartInds(i); %find the start of each trial
            if currIndex+1 <= length(leftTrialData)
                 leftChoice = leftTrialData(currIndex+1,3); 
                 rightChoice = leftTrialData(currIndex+1,4);
                
                if leftChoice == 1 || rightChoice == 1
                    leftCompleteTrialIndex = [leftCompleteTrialIndex; currIndex];
                end
                if leftChoice == 1
                    leftCorrectTrialIndex = [leftCorrectTrialIndex; currIndex];
                end
                if rightChoice == 1 %if right, incorrect choice
                    leftIncorrectTrialIndex = [leftIncorrectTrialIndex; currIndex];
                end
            end
        end

        leftCompleteTrials = leftTrialData(leftCompleteTrialIndex); %SAVE THIS - left trial times
        leftones = ones([length(leftCorrectTrialIndex) 1]);
        correctLeft = [leftTrialData(leftCorrectTrialIndex), leftones]; %SAVE - combine with right
        leftzeros = zeros([length(leftIncorrectTrialIndex) 1]);
        incorrectLeft = [leftTrialData(leftIncorrectTrialIndex), leftzeros]; %SAVE - combine with right

        leftCorrectIncorrect = sortrows([correctLeft; incorrectLeft]); %matrix of ones and zeros for correct/incorrect left trials only




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rightAttempts = sortrows([RightSolenoid; Buzzer]);

        rightsolenoidInds = find(rightAttempts(:,2)); 
        rightRewardAvailInds = [];   
        for i = 1:length(rightsolenoidInds)
            currIndex = rightsolenoidInds(i);
                if currIndex ~= length(rightAttempts)
                buzzerCheck = rightAttempts(currIndex+1,3);
                    if buzzerCheck == 1
                    rightRewardAvailInds = [rightRewardAvailInds; currIndex]; 
                    end
                end
        end

        rightTrials = rightAttempts(:,1:2);
        rightRewardAvail = rightTrials(rightRewardAvailInds,1:2);
        rightRewzeros = zeros([length(rightRewardAvailInds),2]);
        rightRewardAvail = [rightRewardAvail, rightRewzeros; RightWell; LeftWell];
        rightTrialData = sortrows(rightRewardAvail);

        rightTrialStartInds = find(rightTrialData(:,2));
        rightCompleteTrialIndex = [];
        rightCorrectTrialIndex = [];
        rightIncorrectTrialIndex = [];
        for i = [1:length(rightTrialStartInds)]
            currIndex = rightTrialStartInds(i);
            if currIndex+1 <= length(rightTrialData)
                leftChoice = rightTrialData(currIndex+1,3);
                rightChoice = rightTrialData(currIndex+1,4);
                if rightChoice == 1 || leftChoice == 1
                    rightCompleteTrialIndex = [rightCompleteTrialIndex; currIndex];
                end
                if rightChoice == 1
                    rightCorrectTrialIndex = [rightCorrectTrialIndex; currIndex];
                end
                if leftChoice == 1
                    rightIncorrectTrialIndex = [rightIncorrectTrialIndex; currIndex];
                end
            end
        end

        rightCompleteTrials = rightTrialData(rightCompleteTrialIndex); %SAVE THIS- Right trial times
        rightones = ones([length(rightCorrectTrialIndex) 1]);
        correctRight = [rightTrialData(rightCorrectTrialIndex), rightones];
        rightzeros = zeros([length(rightIncorrectTrialIndex) 1]);
        incorrectRight = [rightTrialData(rightIncorrectTrialIndex), rightzeros];

        rightCorrectIncorrect = sortrows([correctRight; incorrectRight]);

        LeftRightPerf = double(sortrows([leftCorrectIncorrect; rightCorrectIncorrect]));
        CorrectIncorrectBinary = LeftRightPerf(:,2);

        
        
        
        triggers.allTriggers = sortrows([leftCompleteTrials; rightCompleteTrials]);
        
        triggers.leftTriggers = leftCompleteTrials;
        triggers.rightTriggers = rightCompleteTrials;
        
        triggers.correctTriggers = sortrows([correctLeft(:,1); correctRight(:,1)]);
        triggers.incorrectTriggers = sortrows([incorrectLeft(:,1); incorrectRight(:,1)]);
        
        odorTriggers{1,day}{1,epoch} = triggers;
        
%         if ~exist([dataDir,'OdorTriggers'])
%             mkdir([dataDir,'OdorTriggers'])
%         end
        %cd([dataDir,'OdorTriggers']);
        save([dataDir,prefix, 'odorTriggers', daystring],'odorTriggers');
        
        if ~exist([dataDir,'BinaryPerf'])
            mkdir([dataDir,'BinaryPerf'])
        end
        cd([dataDir,'BinaryPerf']);
        save(['BinaryPerf',daystring,'-', epstring],'CorrectIncorrectBinary');
        
        
        clear CorrectIncorrectBinary
        cd(dataDir);
        %end
        
    end
    clear odorTriggers
end