function [LinPlaceField,FiresDuringRun,PFexist,RunRates,FieldProps,SmoothOccup] = CalcLinPlaceField(LinCoords,Spikes,varargin)
%function [outputArg1,outputArg2] = CalcLinPlaceField(inputArg1,inputArg2)

% INPUTS:
%   LinCoords: n by 3 vector of ts, linearized position, and velocity
%   Spikes: a vector of spike timestamps
% OPTIONAL:
%   calctheta: empty if you dont want theta data, otherwise its an n by 2
%       vector if ts and LFP phases
%   verbose: show raster and bootstrap data?
%   nBoots: number of bootstrap perms on pf to do
% OUTPUTS:
%   LinField
%   LinOccup
%   PFexist
%   FiresDuringRun
%   RunRates
%   FieldProps
%   PhaseStats

%% InputParser and Preallocated Outputs

p=inputParser;
addOptional(p,'calctheta',false);
addOptional(p,'verbose',false);
addOptional(p,'nBoots',500);
parse(p,varargin{:});
calctheta=p.Results.calctheta;
verbose=p.Results.verbose;
nBoots=p.Results.nBoots;

PFexist=0;
FiresDuringRun=0;
RunRates={}; % this will fill in before silent cells get spit out
LinPlaceField=zeros(1,99); % initialize the rates
SmoothOccup=nan(1,99);
FieldProps=struct('PFmax',nan,'PFmaxpos',nan,'info',nan,...
    'sparsity',nan,'Zpfmax',nan,'PFsize',nan,...
    'PFmaxP',nan,'infoP',nan,'sparsityP',nan);

%% grab all the spikes WITHIN the run epochs
breaks=find(diff(LinCoords(:,1))>1); % only scrub down time that last more than a seconds
runepochs=[[LinCoords(1); LinCoords(breaks+1,1)] [LinCoords(breaks,1); LinCoords(end,1)]];

[Espikes,Erates,~,~,trspikes]=event_spikes(Spikes(:,1),runepochs(:,1),0,diff(runepochs,1,2));
if verbose % show the position raster
    trpos=cellfun(@(a) interp1(LinCoords(:,1), LinCoords(:,2),a,'linear'), trspikes,'UniformOutput',false);
    trpos=trpos(~cellfun(@(a) isempty(a), trpos));
    figure;
    for i=1:length(trpos)
        xvals=linearize(repmat([0;1;nan],1,(length(trpos{i}))))+i;
        plot(linearize(repmat(trpos{i}',3,1)),xvals,'k');
        hold on;
    end
end
RunRates=Erates; % run the splitter score later

if length(Espikes)<20
    %fprintf('no spikes detected ses %d unit %d tr %d (%.2f secs) \n',...
    %ses, j, tr, sum(keepinds)/30);
    return;% if not enough spikes, return to invoking fx
end
% temptraj- 1 is ts, 2 is pos, espikes are epoched spike ts
spikepos=interp1(LinCoords(:,1),LinCoords(:,2),Espikes,'nearest'); % now interp to position
% now you can calculate precession and phase
% preference
if any(calctheta)
    spikephase=interp1(calctheta(:,1),calctheta(:,2),Espikes,'nearest');
    thetastats.spikeprefP(tr)=circ_rtest(spikephase);
    thetastats.spikeprefV(tr)=circ_r(spikephase);
    thetastats.spikeprefrho(tr)=circ_mean(spikephase);
    [thetastats.slope(tr),thetastats.phasestart(tr),~,thetastats.precessionpP]=...
        corrC2Lin_Kempter2012(spikepos,spikephase);
end

% get occupancy map
[bincts,bins]=histcounts(LinCoords(:,2),1:100); % in bins
occupancy=bincts/30; % convert to real time (seconds)
SmoothOccup=SmoothMat2(occupancy,[5 0],2); % smooth over 2 pixels

fullcurves={}; cellsort=[];

if verbose, figure; end
% now for each unit capture a mean tuning curve for each run
% skaggs is run on presmoothed spike and occcupancy
if nBoots>0
    if~isempty(Espikes)
        Bmax=nan(1,nBoots); Binfo=nan(1,nBoots); Bsparsity=nan(1,nBoots);
        for bt=1:nBoots
            EBspikes=trspikes;
            for i=1:length(trspikes)
                % circ shift the spikes in this epoch by a min of 1/2 second
                if ~isempty(trspikes{i})
                    EBspikes{i}=PermuteSpikeTimes(trspikes{i},runepochs(i,:)',0.2);
                end
            end
            Bspikepos=interp1(LinCoords(:,1),LinCoords(:,2),cell2mat(EBspikes'),'nearest');
            bnspikes=histcounts(Bspikepos,bins); % get number of spikes per position
            smoothspikes=SmoothMat2(bnspikes,[5 0],2); % two pixel kernel
            % calculate null info
            [Binfo(bt),Bsparsity(bt)]=Skaggs_basic(SmoothOccup./max(SmoothOccup),...
                smoothspikes,nanmean(smoothspikes));
            % calculate null peak rate
            bLinPlaceField=SmoothMat2(bnspikes./occupancy,[5 0],2);
            if verbose, plot(bLinPlaceField); hold on; end
            Bmax(bt)=max(bLinPlaceField);
            %waitbar(bt/runboots,ha,'working on it ...');
        end
    elseif isempty(Espikes)
        Binfo=nan; Bsparsity=nan; Bmax=nan;
    end
end

[nspikes,bins]=histcounts(spikepos,bins); % get number of spikes per position
smoothspikes=SmoothMat2(nspikes,[5 0],2);
% now save this tuning curve
% smooth the rate map AFTER dividing by raw
% occupancy
LinPlaceField=SmoothMat2(nspikes./occupancy,[5 0],2);

if verbose
    plot(LinPlaceField,'LineWidth',3);
end

%%%%%%%%%%
%calculate all the fieldprops
%%%%%%%%%%
% now qualify this as a pf or not
[FieldProps.PFmax,FieldProps.PFmaxpos]=max(LinPlaceField); % has to be above ??
FieldProps.Zpfmax=max(nanzscore(LinPlaceField)); % has to be above ??

% run skaggs on individually smoothed occupancy and
% spikes
[FieldProps.info,FieldProps.sparsity]=Skaggs_basic(SmoothOccup./max(SmoothOccup),...
    smoothspikes,nanmean(smoothspikes));

% and now pf size e.g. size of top 75%% of max rate
% this will underestimate, because it finds the
% first that drops below the high rate 'flood
% fill' method
pfstart=find(LinPlaceField(1:FieldProps.PFmaxpos) < FieldProps.PFmax*.25,1,'last');
if isempty(pfstart), pfstart=1; end
pfend=FieldProps.PFmaxpos+find(LinPlaceField(FieldProps.PFmaxpos:end)...
    < FieldProps.PFmax*.25,1,'first');
if isempty(pfend), pfend=100; end

% try to calculate pf size, otherwise it doesnt exist
FieldProps.PFsize=pfend-pfstart;

% now add to the session struct



if nBoots>0
    FieldProps.PFmaxP=1-normcdf(FieldProps.PFmax,nanmean(Bmax),nanstd(Bmax)); % peak statistically high?
    FieldProps.infoP=1-normcdf(FieldProps.info,nanmean(Binfo),nanstd(Binfo)); % information statistically high?
    FieldProps.sparsityP=1-normcdf(FieldProps.sparsity,nanmean(Bsparsity),nanstd(Bsparsity)); % Sparsity?
    PFexist=FieldProps.PFmaxP<.05 && FieldProps.PFmax>2 && FieldProps.PFsize<75; % is it a place cell
else
    PFexist=FieldProps.PFmax>2 && FieldProps.Zpfmax>2 && FieldProps.PFsize<75; % is it a place cell
end
FiresDuringRun=FieldProps.PFmax>1; % is this cell active? e.g. does it pass 1 hz at anywhere on track



if verbose % defunct because i folded pfmax etc into fieldprops
    drawnow;
    title(sprintf('Info: %.2f P=%.5f Peak: %.2f, p=%.5f',FieldProps.info,...
        FieldProps.infoP,FieldProps.PFmax,FieldProps.PFmaxP));
    answer = questdlg('Stop verbose?', 'Verbose Menu', 'Yes, quiet','no, more figs','Yes, quiet');
    verbose=contains(answer,'no');
end
end

