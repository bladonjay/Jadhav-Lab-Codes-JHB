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


%% for single session to show raw pcs
win = [-.2 .8]; % Implies from -0.2 to 1 not sure how though
step=0.025; % 10 msec bins
step_ms=step*1000;
binsize=0.15; % 100 msec windows
binctrs=win(1):step:win(2);
binedges=[win(1)-step/2:binsize:win(end)-step/2;...
    win(1)+step/2:binsize:win(end)+step/2];

selectiveonly=0; % only selective units?
PYRonly=0; % only pyramidal cells?
regions={'PFC','CA1'}; % both regions
pcrit=.05;
nboots=400;

i=4; % session 4 is cs33 day 1
events=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
% now only pull reasonably timed trials
okdurs=diff(events,1,2)>.5 & diff(events,1,2)<2;
design=array2table([SuperRat(i).trialdata.leftright10 SuperRat(i).trialdata.CorrIncorr10],...
    'VariableNames',{'LeftRight10','CorrIncorr10'});
events=events(okdurs,:); design=design(okdurs,:);
load("redBlackBlue.mat");
colors1=redBlackBlue(length(redBlackBlue)/2:-1:1,:); cInds=round(rescale(1:length(binctrs),1,length(colors1)));
colors2=redBlackBlue(length(redBlackBlue)/2:end,:);
figure;
for reg=1:2
    % grab our cells
    cellinds=find(contains({SuperRat(i).units.area},regions{reg})); % &...
    %cellfun(@(a) a(3)<pcrit, {SuperRat(i).units.taskResponsive}));

    if PYRonly==1
        cellinds=cellinds & contains({SuperRat(i).units.type},'pyr');
    end
    if selectiveonly==1
        cellinds=cellinds & cellfun(@(a) a(2)<pcrit, {SuperRat(i).units.OdorSelective});
    end
    %
    ratemat={};
    if length(cellinds)>3
        for cl=1:length(cellinds)
            % for each time the matrix is ncells by m trials
            for bn=1:length(binctrs)
                % for start of sniff period
                [~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                    events(design.CorrIncorr10==1,1)+binctrs(bn)-binsize/2,...
                    0,binsize);
                % for end of sniff period
                %[~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                %    events(design.CorrIncorr10==1,2)+binctrs(bn)-binsize/2,...
                %   0,binsize);

            end
        end
        subplot(2,2,reg);
        % some testing here
        oktrials=design.LeftRight10(design.CorrIncorr10==1);
        rlcorr=[]; pcdistBOOT=[]; ctrldist=[]; pcdist=[]; pcraw=[];
        % for each timebin
        for tm=1:length(ratemat)
            % okmat is only rows(cells) that have spikes, and all trials
            % (cols)
            okmat=ratemat{tm}(mean(ratemat{tm}>0,2)>.1,:); % cant be more than 80% empty
            
            % pca takes rows as observations(trials) and cols as
            % variables(cells)
            [~,coeff]=pca(zscore(okmat)');
            if size(coeff,2)>=3
                d=dist([mean(coeff(oktrials==1,1:3)); mean(coeff(oktrials==0,1:3))]');
                % distance between mean values
                pcdist(tm)=d(1,2); % save out one distance
                pcraw(tm,:)=[mean(coeff(oktrials==1,1:3)) mean(coeff(oktrials==0,1:3))];

                % now a shuffle
                for bt=1:nboots
                    oktrialsbt=oktrials(randperm(length(oktrials)));
                    d=dist([mean(coeff(oktrialsbt==1,1:3)); mean(coeff(oktrialsbt==0,1:3))]');
                    pcdistBOOT(bt,tm)=d(1,2);
                end

                % for second control, take half the odor trials and find the
                % dist from one half to another half
                trialmix1=find(oktrials==1);
                trialmix2=find(oktrials==0);
                % even odd split for first odor
                d1=dist([mean(coeff(trialmix1(1:2:end),1:3)); mean(coeff(trialmix1(2:2:end),1:3))]');
                d2=dist([mean(coeff(trialmix2(1:2:end),1:3)); mean(coeff(trialmix2(2:2:end),1:3))]');
                ctrldist(tm)=mean([d1(1,2) d2(1,2)]); % mean dist across rl but for eo split

                if tm>1

                    p=plot3(pcraw(tm-1:tm,1),pcraw(tm-1:tm,2),pcraw(tm-1:tm,3),'.-',...
                        'color',colors1(cInds(tm),:),'LineWidth',2);
                    hold on;
                    p(2)=plot3(pcraw(tm-1:tm,4),pcraw(tm-1:tm,5),pcraw(tm-1:tm,6),'.-',...
                        'color',colors2(cInds(tm),:),'LineWidth',2);

                end
            end

        end
        if ~isempty(pcdist)
            legend(p,{'left trials','right trials'}); xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

            subplot(2,2,reg+2);
            patch([binctrs fliplr(binctrs)],...
                [mean(pcdistBOOT)+std(pcdistBOOT,1)*2 fliplr(mean(pcdistBOOT)-std(pcdistBOOT,1)*2)],...
                [.7 .7 .7],'LineStyle','none','FaceAlpha',.5); hold on;
            p=plot(binctrs,mean(pcdistBOOT),'color',[.3 .3 .3],'LineWidth',2);
            p(2)=plot(binctrs,ctrldist,'r'); p(3)=plot(binctrs,pcdist,'b');
            xlabel('Time from odor offset'); ylabel('Ensemble PC Distance');
            legend(p,{'Null','Even-Odd Control','Right-Left Distance'})



        end
    end
end
%%
%
%
%  This runs across sessions, for this bit, he runs it first locked to
%  start(panel D) and then he uses a precision accuracy... not sure how
%  this is calculated, but something like .95% of trials are correctly
%  classified? as correct?
%
%

%% hardcoded parameters
win = [-1 1]; % Implies from -0.2 to 1 not sure how though
step=0.025; % 10 msec bins
step_ms=step*1000;
binsize=0.12; % 120 msec windows
binctrs=win(1):step:win(2);
binedges=[win(1)-step/2:binsize:win(end)-step/2;...
    win(1)+step/2:binsize:win(end)+step/2];

selectiveonly=0; % only selective units?
PYRonly=0; % only pyramidal cells?
regions={'PFC','CA1'}; % both regions
pcrit=.05;
nboots=400;
PCmat=repmat({[]},length(SuperRat),2);

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
        cellinds=find(contains({SuperRat(i).units.area},regions{reg})); % &...
            %cellfun(@(a) a(3)<pcrit, {SuperRat(i).units.taskResponsive}));

        if PYRonly==1
            cellinds=cellinds & contains({SuperRat(i).units.type},'pyr');
        end
        if selectiveonly==1
            cellinds=cellinds & cellfun(@(a) a(2)<pcrit, {SuperRat(i).units.OdorSelective});
        end
        % 
        ratemat={};
        if length(cellinds)>5
            for cl=1:length(cellinds)
                % for each time the matrix is ncells by m trials
                for bn=1:length(binctrs)
                    % locked to start
                    [~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                        events(design.CorrIncorr10==1,1)+binctrs(bn)-binsize/2,...
                        0,binsize);
                    % locked to end
                    %[~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                    %    events(design.CorrIncorr10==1,2)+binctrs(bn)-binsize/2,...
                    %    0,binsize);
                end
            end

            % some testing here
            oktrials=design.LeftRight10(design.CorrIncorr10==1);
            clear pcdistBOOT pcdist pcraw ctrldist;
            for tm=1:length(ratemat)
                % i want to plot the correlation between the two vectors
                okmat=ratemat{tm}(mean(ratemat{tm}>0,2)>.2,:); % cant be more than 80% empty
                [~,coeff]=pca(ratemat{tm}');
                d=dist([mean(coeff(oktrials==1,1:3)); mean(coeff(oktrials==0,1:3))]');
                % distance between mean values
                pcdist(tm)=d(1,2); % save out one distance
                pcraw(tm,:)=[mean(coeff(oktrials==1,1:3)) mean(coeff(oktrials==0,1:3))];

                % now a shuffle
                for bt=1:nboots
                    oktrialsbt=oktrials(randperm(length(oktrials)));
                    d=dist([mean(coeff(oktrialsbt==1,1:3)); mean(coeff(oktrialsbt==0,1:3))]');
                    pcdistBOOT(bt,tm)=d(1,2);
                end

                % for second control, take half the odor trials and find the
                % dist from one half to another half
                trialmix1=find(oktrials==1);
                trialmix2=find(oktrials==0);
                % even odd split for first odor
                d1=dist([mean(coeff(trialmix1(1:2:end),1:3)); mean(coeff(trialmix1(2:2:end),1:3))]');
                d2=dist([mean(coeff(trialmix2(1:2:end),1:3)); mean(coeff(trialmix2(2:2:end),1:3))]');
                ctrldist(tm)=mean([d1(1,2) d2(1,2)]); % mean dist across rl but for eo split
            end
            % now save out into large matrix
            % matrix is n sessions by 2 region cellmat
            % in each cell is a table: pcdist, ctrldist, nulldist and rows
            % are the timesteps
            PCmat{i,reg}=[pcdist' ctrldist' mean(pcdistBOOT)',...
                mean(pcdistBOOT)'-std(pcdistBOOT,1)'*2 mean(pcdistBOOT)'+std(pcdistBOOT,1)'*2];
        end
    end
    fprintf('Sess %d done \n',i);
end
%% now across session averages
% basically plot means across animals
% 'VariableNames',{'PCdist','ctrldist','nulldist','nullbot','nullupper'});
figure;
regcolors=[0 162/256 181/256 ;242/256 100/256 86/256];
killrows=cellfun(@(a) isempty(a), PCmat);
PCmat(sum(killrows,2)>0,:)=[];
for reg=1:2
    subplot(2,2,reg);
    % first patch
    patch([binctrs fliplr(binctrs)], [mean(cell2mat(cellfun(@(a) a(:,4), PCmat(:,reg), 'UniformOutput', false)'),2)',...
        fliplr(mean(cell2mat(cellfun(@(a) a(:,5), PCmat(:,reg), 'UniformOutput', false)'),2)')],[.7 .7 .7],...
        'LineStyle','none','FaceAlpha',.5);
        hold on;
    % gather the real data
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,3), PCmat(:,reg), 'UniformOutput', false)'),2),'.k-','LineWidth',2);
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,2), PCmat(:,reg), 'UniformOutput', false)'),2),'.r-','LineWidth',2);
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,1), PCmat(:,reg), 'UniformOutput', false)'),2),'.-','LineWidth',2,'color',regcolors(reg,:));
    xlim([-.2 .8])
end
% and using session variance instead of trial variance (i dont know which
% one shantanu used)
figure;
regcolors=[0 162/256 181/256 ;242/256 100/256 86/256];
killrows=cellfun(@(a) isempty(a), PCmat);
PCmat(sum(killrows,2)>0,:)=[];
Vfactor=norminv(.99);
for reg=1:2
    subplot(2,2,reg);
    % first patch
    nullmat=cell2mat(cellfun(@(a) a(:,3), PCmat(:,reg), 'UniformOutput', false)');
    patch([binctrs fliplr(binctrs)], [mean(nullmat,2)+std(nullmat,1,2)*Vfactor;...
        fliplr(mean(nullmat,2)-std(nullmat,1,2)*Vfactor)],[.7 .7 .7],...
        'LineStyle','none','FaceAlpha',.5);
        hold on;
    % gather the real data
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,3), PCmat(:,reg), 'UniformOutput', false)'),2),'.k-','LineWidth',2);
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,2), PCmat(:,reg), 'UniformOutput', false)'),2),'.r-','LineWidth',2);
    plot(binctrs,mean(cell2mat(cellfun(@(a) a(:,1), PCmat(:,reg), 'UniformOutput', false)'),2),'.-','LineWidth',2,'color',regcolors(reg,:));
    xlim([-.2 .8])
end


% this looks very poor, i dont think this is replicable.
            

%% now raw difference in mean rate vectors

win = [-.2 .8]; % Implies from -0.2 to 1 not sure how though
step=0.025; % 10 msec bins
step_ms=step*1000;
binsize=0.1; % 100 msec windows
binctrs=win(1):step:win(2);
binedges=[win(1)-step/2:binsize:win(end)-step/2;...
    win(1)+step/2:binsize:win(end)+step/2];

selectiveonly=0; % only selective units?
PYRonly=0; % only pyramidal cells?
regions={'PFC','CA1'}; % both regions
pcrit=.05;
nboots=400;

i=11; % session 6
events=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
% now only pull reasonably timed trials
okdurs=diff(events,1,2)>.5 & diff(events,1,2)<2;
design=array2table([SuperRat(i).trialdata.leftright10 SuperRat(i).trialdata.CorrIncorr10],...
    'VariableNames',{'LeftRight10','CorrIncorr10'});
events=events(okdurs,:); design=design(okdurs,:);
load("redBlackBlue.mat");
colors1=redBlackBlue(length(redBlackBlue)/2:-1:1,:); cInds=round(rescale(1:length(binctrs),1,length(colors1)));
colors2=redBlackBlue(length(redBlackBlue)/2:end,:);
figure;
for reg=1:2
    % grab our cells
    cellinds=find(contains({SuperRat(i).units.area},regions{reg})); % &...
    %cellfun(@(a) a(3)<pcrit, {SuperRat(i).units.taskResponsive}));

    if PYRonly==1
        cellinds=cellinds & contains({SuperRat(i).units.type},'pyr');
    end
    if selectiveonly==1
        cellinds=cellinds & cellfun(@(a) a(2)<pcrit, {SuperRat(i).units.OdorSelective});
    end
    %
    ratemat={};
    if length(cellinds)>5
        for cl=1:length(cellinds)
            % for each time the matrix is ncells by m trials
            for bn=1:length(binctrs)
                % for start of sniff period
                [~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                    events(design.CorrIncorr10==1,1)+binctrs(bn)-binsize/2,...
                    0,binsize);
                % for end of sniff period
                %[~,ratemat{bn}(cl,:)]=event_spikes(SuperRat(i).units(cellinds(cl)).ts(:,1),...
                %    events(design.CorrIncorr10==1,2)+binctrs(bn)-binsize/2,...
                %   0,binsize);
                
            end
        end
        subplot(2,2,reg);
        % some testing here
        oktrials=design.LeftRight10(design.CorrIncorr10==1);
        rlcorr=[]; clear pcdistBOOT ctrldist pcdist pcraw;
        % for each timebin
        for tm=1:length(ratemat)
            % i want to plot the correlation between the two vectors
            corrmat=zscore(ratemat{tm},1,2);
            pcdist(tm)=mean(abs(mean(corrmat(:,oktrials==1),2)-mean(corrmat(:,oktrials==0),2)));


            % now a shuffle
            for bt=1:nboots
                oktrialsbt=oktrials(randperm(length(oktrials)));
                pcdistBOOT(bt,tm)=mean(abs(mean(corrmat(:,oktrialsbt==1),2)-mean(corrmat(:,oktrialsbt==0),2)));
            end

            % for second control, take half the odor trials and find the
            % dist from one half to another half
            trialmix1=find(oktrials==1);
            trialmix2=find(oktrials==0);
            % even odd split for first odor
            d1=mean(abs(mean(corrmat(:,trialmix1(1:2:end)),2)-mean(corrmat(:,trialmix1(2:2:end)),2)));
            d2=mean(abs(mean(corrmat(:,trialmix2(1:2:end)),2)-mean(corrmat(:,trialmix2(2:2:end)),2)));
            ctrldist(tm)=mean([d1 d2]); % mean dist across rl but for eo split
            
        end

        subplot(2,2,reg+2);
        patch([binctrs fliplr(binctrs)],...
            [mean(pcdistBOOT)+std(pcdistBOOT,1)*2 fliplr(mean(pcdistBOOT)-std(pcdistBOOT,1)*2)],...
            [.7 .7 .7],'LineStyle','none','FaceAlpha',.5); hold on;
        p=plot(binctrs,mean(pcdistBOOT),'color',[.3 .3 .3],'LineWidth',2);
        p(2)=plot(binctrs,ctrldist,'r'); p(3)=plot(binctrs,pcdist,'b');
       xlabel('Time from odor offset'); ylabel('Ensemble PC Distance');
       legend(p,{'Null','Even-Odd Control','Right-Left Distance'})
        
        
        
    end
end
