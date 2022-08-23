function [leftfr,rightfr, lefttrials, righttrials, cellinds, np_left, np_right,np_all] = sj_columnVectors_day_npCells(animal, animno, day, region, win, binsize)

%[topDir,~] = cs_setPaths();
topDir = '/Users/Shantanu/data25/OLF_CS/Data/';
figdir = '/Users/Shantanu/data25/OLF_CS/Data/sj_Figures';

winstr = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
winsize = win(1)+win(2);
leftfr = []; rightfr = [];
cellinds = [];
totalcells = 0;

animDir = [topDir,animal,'_direct/'];
%load([animDir,animal,'cellinfo.mat'])

%cellfilter = ['isequal($area,''',region,'''))'];
% if selectiveonly == 1
%     cellfilter = ['((isequal($area,''',region, ...
%         ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
%     
%     selstr = 'selective';
% else
%     cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region,''')))'];
%     
%     selstr = '';
% end

% Old daycells - by evaluating filter
% -----------------------------------
% cells = evaluatefilter(cellinfo,cellfilter);
% noeps = cells(:,[1 3 4]);
% cells = unique(noeps,'rows');
% daystr = getTwoDigitNumber(day);
% daycells = cells(cells(:,1) == day,:);
% %daynp = nosepokeWindow{day};

% if size(daycells,1) < 3
%     columnVectors = [];
%     return
% end


% Get daycells by looking in npCells structure
% -----------------------------------
load([topDir, 'AnalysesAcrossAnimals/npCells_PyrInt_',region,'.mat'])
daycells_anim = npCells(ismember(npCells(:,[1,2]),[animno, day],'rows'),:); % npCells have [animno day tet cell] - already no epochs

% Remove the animno column from daycells
% -------------------------------------
daycells = daycells_anim(:,[2 3 4]);

daystr = getTwoDigitNumber(day);
load([animDir,animal,'spikes',daystr,'.mat']);
load([animDir,animal,'odorTriggers',daystr,'.mat']);
load([animDir,animal,'nosepokeWindow',daystr,'.mat']);
runeps = find(~cellfun(@isempty,odorTriggers{day}))


%Gather trig times
trigs_all= []; correct_left_all = []; correct_right_all = []; np_all=[];np_left=[];np_right=[];
for ep = 1:length(runeps)
    epoch = runeps(ep);
    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
    trigs = odorTriggers{day}{epoch}.allTriggers;
    
    correct_left_all = [correct_left_all; trigs(correct_left)];
    correct_right_all = [correct_right_all; trigs(correct_right)];
    trigs_all = [trigs_all ; trigs];
    np_all = [np_all; nosepokeWindow{day}{epoch}];
    np_left = [np_left; nosepokeWindow{day}{epoch}(correct_left,:)];
    np_right = [np_right; nosepokeWindow{day}{epoch}(correct_right,:)];
    
end


%Initialize
columnVectors = [];

leftfr = []; rightfr = [];  
cellinds = []; totalcells = 0;
nonrespleft_flag=0; nonrespright_flag=0;
nonrespleftidx=[]; nonresprightidx=[]; % Keep tab of trials with nonresp
%np=[]; leftnp=[]; rightnp=[];

left_popspks_tr=[]; %Send Total Popln Spikes in trial, for Left & Right
right_popspks_tr=[];

for c = 1:size(daycells,1)
    
    leftspikes = []; rightspikes = [];   % Binned spikes
    leftspikes_count = []; rightspikes_count = []; % Total spikes
    
    
    cell = daycells(c,:);
    
    runspikes = [];
    
    for ep = 1:length(runeps)
        epoch = runeps(ep);
        
        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
            runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
        end
        
    end
    
    
    
    if ~isempty(runspikes)
        
        %Find spikes in trigger windows
        for t = 1:length(correct_left_all)
            trig = correct_left_all(t);
            trigwin = [trig-win(1), trig+win(2)];
            winspikes =runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
            bins = (trig-win(1):binsize:trig+win(2));
            if isempty(winspikes)
                %keyboard; break;
                nonrespleft_flag = 1; winspikes=0;
                %nonrespleftidx=[ nonrespleftidx;t];
            end
            binspikecount = histcounts(winspikes,bins);
            leftspikes = [leftspikes;binspikecount];    % Binned spikes in trial
            
            % leftspikes_count = [leftspikes_count;length(winspikes)]; 
            % Total spikes - this is bkcnd+resp, take only resp below
            
            % For firing rate, ignore the background window
            trigwin_resp = [trig, trig+win(2)];
            winspikes_resp =runspikes(runspikes > trigwin_resp(1) & runspikes <= trigwin_resp(2));
            leftspikes_count = [leftspikes_count;length(winspikes_resp)]; % Total spikes   
         
            % Nosepoke
            %currnp_left = np_left_all(t); % Same as trigger time                 
        end
        
        for t = 1:length(correct_right_all)
            trig = correct_right_all(t);
            trigwin = [trig-win(1), trig+win(2)];
            winspikes = runspikes(runspikes > trigwin(1) & runspikes <= trigwin(2));
            if isempty(winspikes)
                %keyboard; break;
                nonrespright_flag = 1; winspikes=0;
                %nonresprightidx=[ nonresprightidx;t];
            end
            bins = (trig-win(1):binsize:trig+win(2));
            binspikecount = histcounts(winspikes,bins);
            rightspikes = [rightspikes;binspikecount];  % Binned spikes
            
            %rightspikes_count = [rightspikes_count;length(winspikes)];  % Total spikes           
            % Total spikes - this is bkcnd+resp, take only resp below
            
            % For firing rate, ignore the background window
            trigwin_resp = [trig, trig+win(2)];
            winspikes_resp =runspikes(runspikes > trigwin_resp(1) & runspikes <= trigwin_resp(2));
            rightspikes_count = [rightspikes_count;length(winspikes_resp)];  % Total spikes  
            
            
             % Nosepoke
            %currnp_right = np_right_all(t);   
            
        end
        
        
        %  sumspikes = sum(sum([leftspikes; rightspikes]));
        %  if sumspikes >= length(trigs_all) %only take the cell if its total spikes during all trig windows is greater than total number of trials. (avg one spike per trial)
        
        sumspikes = sum(sum([leftspikes_count;rightspikes_count]));
        
        if sumspikes >= length(trigs_all) %only take the cell if its total spikes during all trig windows is greater than total number of trials. (avg one spike per trial)
            
            leftfr = [leftfr, leftspikes_count/winsize];
            rightfr = [rightfr, rightspikes_count/winsize];
            %leftnp = [leftnp; currnp_left];
            %rightnp = [rightnp; currnp_right];
            
            totalcells = totalcells+1;
            %cellinds{totalcells,1} = [a, cell];
            cellinds{totalcells,1} = [cell];
            
            lefttrials{totalcells,1}= leftspikes;
            righttrials{totalcells,1} = rightspikes;
            
            left_popspks_tr{totalcells,1} = leftspikes_count;   % Send Nspks for each cell for each trial: Sum across cells will give Popspks/tr
            right_popspks_tr{totalcells,1} = rightspikes_count;
            
            
        end
        
        %v = [leftfr;rightfr];
        
        %columnVectors  = [columnVectors, v];
        
        %                     totalcells = totalcells+1;
        %                     cellinds{totalcells,1} = [a, cell];
        %
        %                     lefttrials{totalcells,1}= leftspikes;
        %                     righttrials{totalcells,1} = rightspikes;
        %                 end
        
        
    end   % runspikes not empty
    
end


% Check if size of trials in firing rate and nosepoke match. Otherwise
% omit, non-resp trials for which no cell fired (unlikely)

if size(np_left,1)~=size(leftfr,1)
    keyboard;
    %nonresp_l = unique(nonrespleftidx); % This does not make sense
    %np_left(nonresp_l,:)=[];
end

if size(np_right,1)~=size(rightfr,1)
    keyboard;
    %nonresp_r = unique(nonresprightidx);
    %np_right(nonresp_r,:)=[];
end
