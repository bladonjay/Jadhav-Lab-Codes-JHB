% Run_PCA_SJ_Copy

%% this runs shantanus pca code and figures but with my datastruct instead
% of claires data.  All you need are the odor periods and the spikes of all
% your cells

%{
Crucial algorithm:

1. for each session, get a 3d matrix.  m trials, n cells and o timebins
    - a design matrix of trial type needs to accompany
    - to get the 3d matrix, grab a 2d mat first for each cell
2. for each session
    for each timebin
        - pca that timebin across trials, then grab the average left and the
        average right first three.  measure euclidean distance between those 2
        points in 3d space  here you can also grab a shuffled distance
3. this should yield a m session by n timebin matrix of real and shuffled
distances
4. 

%}
%% hardcoded parameters
win = [-.5 .25]; % Implies from -0.2 to 1 not sure how though
step=0.01; % 10 msec bins
step_ms=step*1000;
binsize=0.1; % 100 msec windows
binctrs=win(1):step:win(2);
binedges=[win(1)-step/2:binsize:win(end)-step/2;...
    win(1)+step/2:binsize:win(end)+step/2];

selectiveonly=0; % only selective units?
PYRonly=0; % only pyramidal cells?
regions={'PFC','CA1'}; % both regions


for i=1:length(SuperRat)
    % for each session
    % first pull the zero, for me it iwll be momemt from odor end
    events=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    % now only pull reasonably timed trials
    okdurs=diff(events,1,2)>.5 & diff(events,1,2)<2;
    design=array2table([SuperRat(i).trialdata.leftright10 SuperRat(i).trialdata.CorrIncorr10],...
        'VariableNames',{'LeftRight10','CorrIncorr10'});
    events=events(okdurs,:); design=design(okdurs,:);
    
    
    for reg=1:2
        % grab our cells
        cellinds=contains({SuperRat(i).units.area},regions{reg}) &...
            cellfun(@(a) a(3)<pcrit, {SuperRat(i).units.taskResponsive});

        if PYRonly
            cellinds=cellinds & contains({SuperRat(i).units.type},'pyr');
        end
        if selectiveonly
            cellinds=cellinds & cellfun(@(a) a(2)<pcrit, {SuperRat(i).units.OdorSelective});
        end
        mycells=cellfun(@(a) a.ts, {SuperRat.units(cellinds)});
        % 
        cellrates={};
        for cl=1:length(mycells)
            % lock to event end, use bin centers
            for bn=1:length(binctrs)
                cellrates{cl}(bn,:)=event_spikes(mycells(cl).ts,...
                events(design.CorrIncorr20==1,2)+binctrs(bn)-binsize/2,...
                0,binsize);
            end
        end
        cellmat=cat(cellrates,3);
        % save out 3d matrix, with accompanying metadata
        SuperRat(i).PCA.ratemat=cellmat;
        SuperRat(i).PCA.design=design;
    end
end


%% and now to turn into PCA space
% to run first rat, plot an example 3d pc space with r,l, and incorrect r
% and incorrect l

% will want an n timebin by m condition by o pc matrix

i=1;
% design tells you correct, incorrect, right, left
design=SuperRat(i).PCA.design;
[~,~,design.class]=unique(design,'rows');
ratemat=SuperRat(i).PCA.ratemat;
% now for each timebin, pca


% third dimension is window
for tbin=1:size(ratemat,3)
    [coeff,scores,~,~,ev]=pca(ratemat(:,:,tbin));
    % take top explained variance, probably top 3
    % will have to iterate for accumarray
    winpc(tbin,:,1)=accumarray(scores(:,1),design.class,1,@mean,nan);
    winpc(tbin,:,2)=accumarray(scores(:,1),design.class,1,@mean,nan);
    winpc(tbin,:,3)=accumarray(scores(:,1),design.class,1,@mean,nan);
end
figure;
subplot(3,1,1)
% first plot first two pc against time
plot3(winpc(:,1),winpc); % should 
    
% now colorcode all three
subplot(3,1,2);

% now plot euclidean distance between r and l for correct and incorrect
subplot(3,1,3);


            


