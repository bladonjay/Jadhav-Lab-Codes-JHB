% Claire Symanski- March 2 2017

% Extracts correct and incorrect trials from statescript logs from
% pretraining on specified days, and combines into matrix of 1s and 0s for
% use with state space code. 


function cs_stateSpacePerf_preTraining(prefix,days,epochtype)
% days should be cell array of day folder strings matching folder names
% from pretraining


%CS edit 12/28/2017 - days should now be in format of day number, i.e 1:10 

% rawdir = ['E:\Data\PreTraining\',prefix,'\'];
% datadir = ['E:\Data\OdorPlaceAssociation\',prefix,'_direct\Behavior\'];

rawdir = ['D:\OdorPlaceAssociation\',prefix,'Expt\PreTraining\'];
datadir = ['D:\OdorPlaceAssociation\',prefix,'Expt\PreTraining\'];
topDir = ['D:\OdorPlaceAssociation\',prefix,'Expt\'];

cd(rawdir)

files = dir();
filenames = {files(3:end).name};

BinaryPerfAll = [];
for d = 1:length(days)
    day = days{d};
    dayfolderind = strfind(filenames, day);
    dayfoldername = filenames{1,not(cellfun('isempty', dayfolderind))};
    dayfolder = [rawdir, dayfoldername,'\'];
    cd(dayfolder)

    ssLogs = dir('*stateScriptLog');
    
    for f = 1:length(ssLogs)
        
        ssLogName = ssLogs(f).name;
        
        k = strfind(ssLogName,epochtype);
        if ~isempty(k) 
                
        
        [BinaryPerf] = cs_getBinaryPerf_preTraining(ssLogName);
        
        BinaryPerfAll = [BinaryPerfAll; BinaryPerf];
        
        end
    end
        
end

%save([datadir,'binaryPerfPreTraining.mat'],'BinaryPerfAll');

getestprobcorrect_niceplot(BinaryPerfAll, 0.5, 0); 

cd(datadir)
% savename = [datadir,'PreTrainingBinaryPerf'];
%     save(savename,'BinaryPerfAll');
    