%CS- modified from DFSsj_HPexpt_placefield2_savefields

%Linfields:
%1st column is distance in cm
%2nd column is occupancy
%3rd column is spike count
%4th column is occ-normalized firing rate
%5th column is smoothed occ-normalized firing rate <-- Use this one
%6th column is smoothed occupancy
%7th column is smoothed spike count

clear;
runscript = 1;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options
savelinfields = 1; % To save trajdata in [prefix-linfields-day] file for each day
% Also mapdata in [prefix-mapfields-day] file for day
% Save mapdata as opnfield rate, as well as separate trajectories
topDir = cs_setPaths();
% savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
% val=1; savefile = [savedir 'HPa_allfields']; % HPa fields
%val=2; savefile = [savedir 'HPb_allfields']; % HPb fields
%val=3; savefile = [savedir 'HPc_allfields']; % HPc fields
%val=4; savefile = [savedir 'Ndl_allfields']; % Ndl fields


minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;

plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
plotanimidx = []; % To pick animals for plotting


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
%     switch val
%         case 1
%             animals = {'HPa'}
%         case 2
%             animals = {'HPb'}
%         case 3
%             animals = {'HPc'}
%         case 4
%             animals = {'Ndl'}
%     end
    animals = {'CS31','CS33','CS34','CS35','CS44'};
    %animals = {'CS44'};
    %region = 'CA1';
    %Filter creation
    %--------------------------------------------------------
    
    % Epoch filter
    % -------------
    %dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    %epochfilter = 'isequal($type, ''run'')';
    runepochfilter = 'isequal($environment, ''odorplace'')';
    
    
    % Cell filter
    % -----------
    
    % All cells - no condition
    % -----
    % cellfilter = '( strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') ||  strcmp($area, ''PFC'') ) ';
    
    % All cells with Nspk consition
    % ------
    cellfilter = '(strcmp($area, ''CA1'') || strcmp($area, ''PFC'')) && (strcmp($type, ''pyr'')) && ($numspikes > 100)';
    % cellfilter = '((strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7))  ||  (strcmp($area, ''PFC'') && ($numspikes > 100) ) ';
    
    %placecellfilter = '( strcmp($tag, ''PFC'') && ($numspikes > 100))';
    %placecellfilter = '( strcmp($tag, ''CA1Pyr'') || strcmp($tag, ''iCA1Pyr'') || strcmp($tag, ''PFC'') )';
    %placecellfilter = '( strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''iCA1Pyr'') || strcmp($tag2, ''PFC'') ) && ($numspikes > 100)';
    
    
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
    
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Abs linear velocity(thrs 5)/ velocity(thrs 3) time filter
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} }
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',[],'tetfilter',riptetfilter,'minthresh',2} }
    
    % Linear velocity time filter implemented in getlinstate
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} };
    
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    %timefilter = { {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Note - For includestates=1:6, state will never be -1
    
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    psf = createfilter('animal',animals,'epochs',runepochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator', iterator);
    % psf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ----------------------
    psf = setfilterfunction(psf, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 3);
    %psfold = setfilterfunction(psf, 'DFAsj_filtercalclinfields', {'spikes', 'linpos', 'pos'}, 'binsize', 2, 'minabsvel', 3);
    pmf = setfilterfunction(psf, 'DFAsj_openfieldrate_tf', {'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 3);
    %pof = setfilterfunction(psf, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 2);
    
    %traj IDs: outbound L = 1, inbound L = 2, outbound R = 3, outbound R =4
    
    disp('Finished filter creation');
    
    % Run analysis
    % ------------
    psf = runfilter(psf);  % Place Field Stability
    pmf = runfilter(pmf);  % Place Field Map
    %pof = runfilter(pof);  % Place Field Map - Separate trajectories
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata savelinfields
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

% if ~exist('savedata')
%     return
% end

% ----------------------------------

gatherdata = 1; savegatherdata = 1;

% Get trajdata and days and epochs
trajdata = []; index=[];
allanimindex=[]; allmaps=[]; alltrajs=[];

for an = 1:length(psf)
    for i=1:length(psf(an).output{1})
        index{an}(i,:)=psf(an).output{1}(i).index; %cellindex
        alltrajdata{an}{i}=psf(an).output{1}(i).trajdata;
        allmapdata{an}{i}=pmf(an).output{1}(i);
        %allmapdata_sep{an}{i}=pof(an).output{1}(i);
        
        % Only indexes
        animindex=[an psf(an).output{1}(i).index]; % Put animal index in front
        allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
        allmaps{i} = pmf(an).output{1}(i);
        alltrajs{i} = psf(an).output{1}(i).trajdata;
    end
end


% For each index, save linfields and mapfields
% --------------------------------------------

for an = 1:length(psf)
    prefix = psf(an).animal{1};
    animdirect = psf(an).animal{2};
    
    curranimidxs = index{an}; % All indexes for cells from current animal. day-ep-tet-cell
    uniquedays = unique(curranimidxs(:,1));
    
    for d = 1:length(uniquedays)
        day = uniquedays(d);
        dayidxs = find(curranimidxs(:,1)==day);
        daystr = getTwoDigitNumber(day);
        curranimdayidxs = curranimidxs(dayidxs,:);
        linfields = []; mapfields = []; % Reset linfields and mapfields
        
        for c = 1:length(curranimdayidxs)
            curridx = curranimdayidxs(c,:);
            currmapidx = find( index{an}(:,1)==curridx(1) & index{an}(:,2)==curridx(2) & index{an}(:,3)==curridx(3) & index{an}(:,4)==curridx(4));
            linfields{curridx(1)}{curridx(2)}{curridx(3)}{curridx(4)}=alltrajdata{an}{currmapidx};
            mapfields{curridx(1)}{curridx(2)}{curridx(3)}{curridx(4)}=allmapdata{an}{currmapidx}; % Save everything .smoothedspikerate is what you want
        end
        
        % Save for current day
         if savelinfields==1
            
            save([animdirect,prefix,'linfields',daystr],'linfields');
            save([animdirect,prefix,'mapfields',daystr],'mapfields');
            
         end  
    end % end day
end % end an


        
