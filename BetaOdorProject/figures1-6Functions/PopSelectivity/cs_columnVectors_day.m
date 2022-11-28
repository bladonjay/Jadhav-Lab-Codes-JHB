function [leftfr,rightfr] = cs_columnVectors_day(animal, day, region, win, selectiveonly)

[topDir,~] = cs_setPaths();
winstr = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
winsize = win(1)+win(2);
leftfr = []; rightfr = [];
cellinds = [];
totalcells = 0;

animDir = [topDir, animal, 'Expt\',animal,'_direct\'];

load([animDir,animal,'cellinfo.mat'])

%cellfilter = ['isequal($area,''',region,'''))'];

if selectiveonly == 1
    cellfilter = ['((isequal($area,''',region, ...
        ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
    
    selstr = 'selective';
else
    cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region,''')))'];
    
    selstr = '';
end

cells = evaluatefilter(cellinfo,cellfilter);

noeps = cells(:,[1 3 4]);
cells = unique(noeps,'rows');


daystr = getTwoDigitNumber(day);

daycells = cells(cells(:,1) == day,:);

if size(daycells,1) < 3
    columnVectors = [];
    return
end
load([animDir,animal,'spikes',daystr,'.mat'])
load([animDir,animal,'odorTriggers',daystr,'.mat'])
runeps = find(~cellfun(@isempty,odorTriggers{day}));


%Gather trig times
trigs_all= []; correct_left_all = []; correct_right_all = [];
for ep = 1:length(runeps)
    epoch = runeps(ep);
    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
    trigs = odorTriggers{day}{epoch}.allTriggers;
    
    correct_left_all = [correct_left_all; trigs(correct_left)];
    correct_right_all = [correct_right_all; trigs(correct_right)];
    trigs_all = [trigs_all ; trigs];
    
end

columnVectors = [];

leftfr = [];
    rightfr = [];
for c = 1:size(daycells,1)
    
    leftspikes = [];
    rightspikes = [];
    
    
    cell = daycells(c,:);
    
    runspikes = [];
    
    for ep = 1:length(runeps)
        epoch = runeps(ep);
        
        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
            runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
        end
        
        
    end
    
    %Find spikes in trigger windows
    for t = 1:length(correct_left_all)
        trig = correct_left_all(t);
        trigwin = [trig-win(1), trig+win(2)];
        winspikes =runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
        
        leftspikes = [leftspikes;length(winspikes)];
    end
    
    for t = 1:length(correct_right_all)
        trig = correct_right_all(t);
        trigwin = [trig-win(1), trig+win(2)];
        winspikes = runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
        
        rightspikes = [rightspikes;length(winspikes)];
    end
    
    
    %                 sumspikes = sum(sum([leftspikes; rightspikes]));
    %
    %                 if sumspikes >= length(trigs_all) %only take the cell if its total spikes during all trig windows is greater than total number of trials. (avg one spike per trial)
    leftfr = [leftfr, leftspikes/winsize];
    rightfr = [rightfr, rightspikes/winsize];
    
    %v = [leftfr;rightfr];
    
    %columnVectors  = [columnVectors, v];
    
    %                     totalcells = totalcells+1;
    %                     cellinds{totalcells,1} = [a, cell];
    %
    %                     lefttrials{totalcells,1}= leftspikes;
    %                     righttrials{totalcells,1} = rightspikes;
    %                 end
    
end

