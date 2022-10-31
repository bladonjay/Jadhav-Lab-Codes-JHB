function out = DFAcssj_openfieldrate_leftrightall(index, excludeperiods, spikes, linpos, pos, runTrialBounds, varargin)

% Equivalent to sj_openfieldrate and DFsj_openfieldrate - WITH TIMEFILTER
% DO NOT NEED POS for this version UNLIKE DFsj_openfieldrate since
% speed is supposed to have been filtered out by timefilter
%
% Does NOT separate by trajectory like DFAsj_twodoccupancy or twodoccupancy
%
% Shantanu: Adding Abs Vel from "pos" field. Accurate calculation of
% velocity, and independent of minvel inputed/used in sj_getbehavestate
% 13 Jun 2011

%[out] = openfieldratemap(index, excludeperiods,spikes,behavestate,pos, binsize, std, minabsvel)
%[out] = openfieldratemap(index, excludeperiods,spikes,behavestate,pos, binsize)
%
%Calculates the 2d occupancy normalized firing rate for the cell
%
%spikes - the spike data structure 
%behavestate - the state data structure
%pos - the position data structure
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%
%The output is a structure with n matrices. The matrices are: 
% 1) occupancy, 
% 2) bin vector x, 
% 3) bin vector y, 
% 4) bin spike count, 
% 5)occ normailized firing per bin, and 
% 6) smoothed occ normalized firing. 
%

% parse the options
appendindex = 0;
binsize = 1; % cm square 
std = 2; 
threshocc = 0.02; % Threshold occupancy in seconds

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'threshocc'
                std = varargin{option+1};
            case 'trialtype'
                trialtype = varargin{option+1}; %CS added 8/14/2017 - all trials = 2, left = 1, right = 0 (default is all trials, -1)
	    otherwise
                error(['Option ',varargin{option},' unknown.']);
	end   
    else
        error('Options must be strings, followed by the variable');
    end
end


warning('OFF','MATLAB:divideByZero');

disp(['Doing Day ', num2str(index(1)), ' epoch ', num2str(index(2)), ' tetrode ', num2str(index(3)), ' cell ', num2str(index(4))]);

% Spikes
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
% Position Data
posdata = pos{index(1)}{index(2)}.data;
runTrialBounds = runTrialBounds{index(1)}{index(2)}.data;

statematrix = linpos{index(1)}{index(2)}.statematrix;
lindist = statematrix.lindist;
statevector = statematrix.traj;

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
% that already exists based on behavestate
% statematrix.time and posdata(:,1) are equivalent
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time, 



% spikes.data fields: time x y dir not_used amplitude(highest variance channel) posindex x-sm y-sm dir-sm
if size(spikes,2)>7
    posx = 8; posy = 9;
else
    posx = 2; posy = 3;
end
posx = 2; posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);
    % Find any posidxs that are outside the range of linpos (statevector)
    outliers = find(spikes(:,4) > length(statevector));
    spikes(outliers,:)=[];
    spikes(:,5) = statevector(spikes(:,4)); %add the traj state for each spike
    %spikes(:,6) = absvel(spikes(:,4)); %add the absvel for each spike
else
    spikes = [0 0 -1];
end


timestep = posdata(2,1) - posdata(1,1);


    goodspikes = [];
    goodpos = [];
    if trialtype ~= 2
         leftTrials=runTrialBounds((runTrialBounds(:, 3) == 1 & runTrialBounds(:, 4) == 1) | (runTrialBounds(:, 3) == 0 & runTrialBounds(:, 4) == 0), :);
         rightTrials=runTrialBounds((runTrialBounds(:, 3) == 0 & runTrialBounds(:, 4) == 1) | (runTrialBounds(:, 3) == 1 & runTrialBounds(:, 4) == 0), :);
        
         if trialtype == 0
             trialbounds = rightTrials;
         end
         if trialtype == 1
             trialbounds = leftTrials;
         end
%          specifiedtrials = (runTrialBounds(:,3) == trialtype); %gets only left or right trials
%         trialbounds = runTrialBounds(specifiedtrials,:);
        for trial = 1:length(trialbounds(:,1))
            goodspikes = [goodspikes; spikes(( spikes(:,1) >= trialbounds(trial,1) & spikes(:,1) <= trialbounds(trial,2)),:)];
            goodpos = [goodpos; posdata(( posdata(:,1) >= trialbounds(trial,1) & posdata(:,1) <= trialbounds(trial,2)),:)];
        end
    else 
        goodspikes = spikes;
        goodpos = posdata;
    end

    %index(3),index(4), size(spikes)


    % Non sharp-wave spikes
    % Shantanu - I don't use Column 5
    %goodspikeind = (spikes(:,5) ~= -1);
    %goodspikes = spikes(goodspikeind,:);

    %make a cell array, where each cell contains data for one trajectory.
    %inside each cell [binlocation occupancy spikecount firingrate]

    tmpposition = (goodpos(:,[2 3]));
    tmpspikes = (goodspikes(:,[2 3]));

    if ~isempty(tmpposition)
        minx = floor(min(tmpposition(:,1))) - 1;
        maxx = ceil(max(tmpposition(:,1))) + 1;
        binx = (minx:binsize:maxx);
        miny = floor(min(tmpposition(:,2))) - 1;
        maxy = ceil(max(tmpposition(:,2))) + 1;
        biny = (miny:binsize:maxy);
        %[out.occupancy out.xticks out.yticks] = hist2(tmpposition, binx, biny);
        % Shantanu - I am using hist2 in ~/Src/Matlab/sj_Utility
        [output.occupancy, output.xticks, output.yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
        %[out.occupancy out.xticks out.yticks] = hist2d(tmpposition, min([length(binx),length(biny)]));

        if ~isempty(tmpspikes)

            % From calcopenfieldoccupancy / DFAsj_twodoccupancy

            %1) Get spikerate from spikes and occupancy
            [output.spikes, BX, BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
            nonzero = find(output.occupancy ~= 0);
            output.spikerate = zeros(size(output.spikes));
            output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );

            %2) Smooth spikerate and occupancy
            g = gaussian2(std,(6*std));
            output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
            %smoothedoccupancy = [];
            %smoothedoccupancy = zeros(size(output.spikes));
            smoothedoccupancy = filter2(g, output.occupancy);

            % 3) Turn occupancy to seconds and set spikerate wherever occupancy
            %is < threshold occupancy in seconds to 0

            output.occupancy = timestep*output.occupancy;
            output.smoothedoccupancy = timestep*smoothedoccupancy;

            %zero = find(smoothedoccupancy == 0);
            zero = find(output.smoothedoccupancy <= threshocc);
            %output.smoothedspikerate(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
            output.smoothedspikerate(zero) = -1; 

        else

            output.spikes = [];
            output.spikerate = [];
            output.smoothedspikerate = [];
            output.smoothedoccupancy = [];

        end
    end


if appendindex == 0
    out = output;
elseif appendindex == 1
    output.index = index;
    out = output;
end

warning('ON','MATLAB:divideByZero');


