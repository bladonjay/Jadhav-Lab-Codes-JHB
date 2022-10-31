clear
%% Params
animal = 'CS44';
%topDir = 'F:\Data\OdorPlaceAssociation\';
topDir = 'D:\OdorPlaceAssociation\';
animDir = [topDir, animal, 'Expt\', animal, '_direct\'];
day = 5;
epoch = 2;
region = 'PFC';

%% Gather Data
daystr = getTwoDigitNumber(day);

load([animDir,animal,'cellinfo.mat']);
load([animDir,animal,'odorTriggers',daystr,'.mat']);
load([animDir,animal,'spikes',daystr,'.mat']);

lefttrigs = odorTriggers{day}{epoch}.leftTriggers;
righttrigs = odorTriggers{day}{epoch}.rightTriggers;
numtrigs = length(lefttrigs)+length(righttrigs);

trajlabel = [ones(1,length(lefttrigs)), repmat(2,1,length(righttrigs))];

[c_left, c_right, i_left, i_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});

c_left = [c_left, ones(length(c_left),1)];
i_left = [i_left, zeros(length(i_left),1)];
ci_left = sortrows([c_left; i_left],1); ci_left = ci_left(:,2);

c_right = [c_right, ones(length(c_right),1)];
i_right = [i_right, zeros(length(i_right),1)];
ci_right = sortrows([c_right; i_right],1); ci_right = ci_right(:,2);

trajcorrect = [ci_left; ci_right]';


epochdata = spikes{day}{epoch};


%cellfilter = ['isequal($area,''', region,''') && strcmp($type, ''pyr'') && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective''))'];

cellfilter = ['isequal($area,''', region,''') && strcmp($type, ''pyr'') && $numspikes > 0'];
cells = evaluatefilter(cellinfo{day}{epoch},cellfilter);

%% Get matrix of FR data

frmatrix = zeros(size(cells,1), numtrigs);
for c = 1:size(cells, 1)
    cind = cells(c,:);
    
    spikes = epochdata{cind(1)}{cind(2)}.data(:,1);
    
    for lt = 1:length(lefttrigs)
        trig = lefttrigs(lt);
        win = [trig, trig+1];
        %spikesinwin = spikes(spikes >= win(1) & spikes < win(2));
        spikesinwin = spikes(isExcluded(spikes, win));
        winfr = sum(spikesinwin)/(win(2)-win(1));
        frmatrix(c,lt) = winfr;
    end
    
    for rt = 1:length(righttrigs)
        trig = righttrigs(rt);
        win = [trig, trig+1];
        %spikesinwin = spikes(spikes>= win(1) & spikes < win(2));
        spikesinwin = spikes(isExcluded(spikes, win));
        winfr = sum(spikesinwin)/(win(2)-win(1));
        frmatrix(c,length(lefttrigs)+rt) = winfr;
    end

end
frmatrix = frmatrix(any(frmatrix,2),:); %remove rows with all zeros

%% Do MDS
opts = statset('MaxIter',1000);
dissimilarities = pdist(frmatrix','cityblock');
Y = mdscale(dissimilarities,3,'criterion','sammon','Options',opts);
signals = Y(1:length(frmatrix(1,:)),:)';
% signals_incorr = Y(1+length(frmatrix(1,:)):end,:)';

%% Plotting

figure('color','w')
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultaxesfontsize',14);

labelsall = [trajlabel;trajcorrect]';

%correct left
inds = ismember(labelsall,[1 1],'rows');
plot3(signals(1,inds),signals(2,inds),signals(3,inds),...
    'o','MarkerEdgeColor','w',...
    'MarkerFaceColor','b',...
    'MarkerSize',10)
hold on

%incorrect left
inds = ismember(labelsall,[1 0],'rows');
plot3(signals(1,inds),signals(2,inds),signals(3,inds),...
    'o','MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',10)
hold on

%correct right
inds = ismember(labelsall,[2 1],'rows');
plot3(signals(1,inds),signals(2,inds),signals(3,inds),...
    'o','MarkerEdgeColor','w',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
hold on


%incorrect right
inds = ismember(labelsall,[2 0],'rows');
plot3(signals(1,inds),signals(2,inds),signals(3,inds),...
    'o','MarkerEdgeColor','r',...
    'MarkerFaceColor','w',...
    'MarkerSize',8)

legend({'CorrectLeft','IncorrectLeft','CorrectRight','IncorrectRight'})

% for traj_seg = find(trajlabel == 1) %left trials
%     if trajcorrect(traj_seg) == 1 %correct left
%         plot3(signals(1,traj_seg),signals(2,traj_seg),signals(3,traj_seg),...
%             'o','MarkerEdgeColor','w',...
%                 'MarkerFaceColor','b',...
%                 'MarkerSize',10)
%         hold on
%     else %incorrect left
%         plot3(signals(1,traj_seg),signals(2,traj_seg),signals(3,traj_seg),...
%             'o','MarkerEdgeColor','b',...
%                 'MarkerFaceColor','w',...
%                 'MarkerSize',8)
%         hold on
%     end
% end
% 
% for traj_seg = find(trajlabel == 2) %right trials
%     if trajcorrect(traj_seg) == 1 %correct right
%         plot3(signals(1,traj_seg),signals(2,traj_seg),signals(3,traj_seg),...
%              'o','MarkerEdgeColor','w',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10)
%         hold on
%     else %incorrect right
%         plot3(signals(1,traj_seg),signals(2,traj_seg),signals(3,traj_seg),'rx')
%         plot3(signals(1,traj_seg),signals(2,traj_seg),signals(3,traj_seg),...
%             'o','MarkerEdgeColor','r',...
%                 'MarkerFaceColor','w',...
%                 'MarkerSize',8)
%         hold on
%     end
% end




