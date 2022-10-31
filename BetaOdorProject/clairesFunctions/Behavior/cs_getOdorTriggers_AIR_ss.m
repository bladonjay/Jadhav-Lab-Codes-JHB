function cs_getOdorTriggers_AIR_ss(animal, days,epochs)
%topDir should be top Expt folder, WITH trailing slash

%This function was written for getting DIO triggers from statescript
%instead of DIO files. Written to get around DIO noise from Trodes 1.6.2
%with headstages with sensors, should be fixed in newer trodes versions.
%not recommended for normal use (use cs_getOdorTriggers.m instead, uses
%DIO.mat files).

%All days being run at the same time MUST have same number and order of run
%epochs (i.e. all days have run session epochs of 2,4,6)
topDir = cs_setPaths();
animDir = [topDir,animal,'Expt\'];
rawDir = [animDir,animal,'\'];
cd(rawDir)
dayfolders = dir() ;
    dayfolders = {dayfolders(3:end).name};

    for d = 1:length(days)
        runEps = {};
    day = days(d);
    daystring = getTwoDigitNumber(day);
    cd(dayfolders{day})
        ssfiles = dir('*.stateScriptLog');
        for i = 1:length(ssfiles)
            epochname = ssfiles(i).name;
            k = strfind(epochname,'None');
            if ~isempty(k)
                runEps{end+1} = epochname;
            end
        end
            
   
    for e = 1:length(runEps)
        filename = runEps{e}; %filename to use for SS file
        
        epstring = getTwoDigitNumber(epochs(e));
        
        %if ep <= numEpochs
        stateScriptStruct = parseStateScriptFile_16Ch(filename);

%         for i= 1:length(stateScriptStruct)
%            want_idx(i)= (~isempty(stateScriptStruct(i).inState)) ;
%         end
%         stateScriptStruct = stateScriptStruct(want_idx);
        rowstoremove = [];
        for s = 1:length(stateScriptStruct)

            ssinState = stateScriptStruct(s).inState;
            if isempty(ssinState)
                rowstoremove(end+1) = s;
            end
        end

        stateScriptStruct(rowstoremove) = [];

        time = cat(1,stateScriptStruct.time);
        time = time/1000; %convert to seconds from ms
        outStates = [time, cat(1,stateScriptStruct.outState)];
        inStates = [time, cat(1,stateScriptStruct.inState)];

        Din1hits = inStates((inStates(:,2)==1),[1 2]);
        Din1zeros = zeros(length(Din1hits),1);
        LeftWell = [Din1hits(:,1), Din1zeros, Din1hits(:,2), Din1zeros];

        Din2hits = inStates((inStates(:,3)==1),[1 3]);
        Din2zeros = zeros(length(Din2hits),2);
        RightWell = [Din2hits(:,1), Din2zeros, Din2hits(:,2)];

        Dout9hits = outStates((outStates(:,10)==1),[1 10]);
        Dout9zeros = zeros(length(Dout9hits),1);
        LeftSolenoid = [Dout9hits(:,1), Dout9hits(:,2), Dout9zeros];

        Dout10hits = outStates((outStates(:,11)==1),[1 11]);
        Dout10zeros = zeros(length(Dout10hits),1);
        RightSolenoid = [Dout10hits(:,1), Dout10hits(:,2), Dout10zeros];

        Dout5hits = outStates((outStates(:,6)==1),[1 6]);
        Dout5zeros = zeros(length(Dout5hits),1);
        Buzzer = [Dout5hits(:,1), Dout5zeros, Dout5hits(:,2)];

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
                    
                    if leftChoice == 1
                    leftCorrectTrialIndex = [leftCorrectTrialIndex; currIndex];
                    end
                    if rightChoice == 1 %if right, incorrect choice
                    leftIncorrectTrialIndex = [leftIncorrectTrialIndex; currIndex];
                    end
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
                    if rightChoice == 1
                    rightCorrectTrialIndex = [rightCorrectTrialIndex; currIndex];
                    end
                    if leftChoice == 1
                    rightIncorrectTrialIndex = [rightIncorrectTrialIndex; currIndex];
                
                    end
                
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
        
        odorTriggers{1,day}{1,epochs(e)} = triggers;
        
        
        
            if e<size(runEps,1)
                cd([animDir,'\',dayfolders{day}]);
            end
            
        
    end
        
        
        cd(animDir)
    end
    
    dataDir =([animDir,animal,'_direct\']);
        cd(dataDir)
%         if ~exist([dataDir,'OdorTriggers'])
%             mkdir([dataDir,'OdorTriggers'])
%         end
%         cd([dataDir,'OdorTriggers']);
        save([animal, 'airOdorTriggers', daystring],'odorTriggers');
        
%         if ~exist([dataDir,'BinaryPerf'])
%             mkdir([dataDir,'BinaryPerf'])
%         end
%         cd([dataDir,'BinaryPerf']);
        save([dataDir,'BinaryPerf\AirBinaryPerf',daystring,'-', epstring],'CorrectIncorrectBinary');
        
    clear odorTriggers
    clear CorrectIncorrectBinary
end