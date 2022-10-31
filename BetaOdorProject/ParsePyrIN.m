% ParsePyrIN

% this script contains code to parse pyramidal cells and interneurons based
% on the burst index in the autocorrelogram and the overall firing rate



%% this generates the ac and calculates the stats off that
% this adds the short ISI's for 1 and 2 msec, it adds the numspikes, the
% mean rate, the burst probability, the burst index, and the theta index

%% heres where I make the designation between pyrams and ins for each region

% this is exactly how claire ran this analysis except that she had a method
% for parsing cells whose rate was between 9 and 9.5 by using waveforms...
% she never passed me waveform data so i am just letting those cells
% through

for ses= 1:length(SuperRat)
    for j=1:length(SuperRat(ses).units)
        if contains(SuperRat(ses).units(j).area,'CA1')
            % I think we need to positively select PYRs and leave ins
            % to whatever is left
            if SuperRat(ses).units(j).meanrate<=9.5
                SuperRat(ses).units(j).type='pyr';
            else
                SuperRat(ses).units(j).type='in';
            end
            
        elseif contains(SuperRat(ses).units(j).area,'PFC')
            if SuperRat(ses).units(j).meanrate<9.5
                SuperRat(ses).units(j).type='pyr';
            else
                SuperRat(ses).units(j).type='in';
            end
            
        end        
    end
end



























%%
% dont use this code, its over the top...
% it gives you a deep dive into burstiness, acgrams etc.
%{
verbose=0;

for ses = 1:length(SuperRat)
    tic
    tracking=SuperRat(ses).tracking.data(:,1);
    elapsed=diff(tracking);
    sesslength=sum(elapsed(elapsed<1)); % sum elapsed time thats less than a second
    for j=1:length(SuperRat(ses).units)
        %if isempty(SuperRat(ses).units(j).tag)
        %    SuperRat(ses).units(j).tag='accepted'; % probably not best practice, but whatever
        %end
        
        spikets=SuperRat(ses).units(j).ts(:,1);

        % first count isi violations (1ms then 2 ms)
        SuperRat(ses).units(j).shortISI=sum(diff(spikets)<=0.001);
        SuperRat(ses).units(j).shortISI(2)=sum(diff(spikets)<=0.002);
        
        % mean rate and number of spikes
        SuperRat(ses).units(j).numspikes=length(spikets);
        SuperRat(ses).units(j).meanrate=length(spikets)/sesslength;
        
        % now for the AC stuff;
        timebin=1/1000; % 1 ms per bin
        tseries=accumarray(ceil(spikets/timebin),1); % acgram it at 1 ms steps
        acgram=xcorr(tseries,tseries,50);
        acgram(51)=nan;
        
        if sum(acgram(1:50))>30
            burstindex=max(acgram(40:50))/nanmean(acgram(1:10));
        else
            burstindex=nan;
        end
        isi=diff(spikets); % get isis
        closeisi=[isi<.009;0]; % get short pre isis
        closeisi(:,2)=circshift(closeisi,1); % get short post isis
        % burst prob is % spikes that have a short pre or post isi
        burstprob=mean(sum(closeisi,2)>0);
        SuperRat(ses).units(j).propbursts=burstprob;
        SuperRat(ses).units(j).burstIndex=burstindex;
        
        % now get theta in acgram
        timebin=1/100; acgramdist=70;
        step=1000*timebin;
        x=(-step*acgramdist):step:(step*acgramdist);
        tseries=accumarray(ceil(spikets/timebin),1);
        acgram2=xcorr(tseries,tseries,acgramdist);
        % nan out the middle point
        acgram2(acgramdist-1:acgramdist+1)=nan;
        acgram2=acgram2/max(acgram2);
        SuperRat(ses).units(j).ACgram=acgram2;
        % mean out the center spikes (+10, 0, & -10 msec centers) which will bin
        % out to +15 to -15
        acgram2(acgramdist-1:acgramdist+1)=nanmean(acgram2(acgramdist-4:acgramdist+4));
        
        royermodel= @(a, b, c, f, tau1, tau2, x)...
            (a*(cos(2*.001*pi*f*x)+1)+b).*exp(-abs(x)/tau1)+... % theta with the slow decay
            c*exp(-(x.^2)/tau2^2); % plus fast decay
        
        FreqBand=[6 12]; % theta!
        
        [royermx] = fit(x', acgram2, royermodel, ...
            'StartPoint', [1, .1,   2,   nanmean(FreqBand),     500,  25], ...
            'Lower',      [0,   0, -10,   FreqBand(1),          75,   0], ...
            'Upper',     [10,  10,  10,  FreqBand(2),           2000, 75],...
            'MaxFunEvals',10^10,'TolFun',10^-10, 'Robust', 'LAR');
        royercf=coeffvalues(royermx);
        SuperRat(ses).units(j).ACtheta=royercf(2)/royercf(1); % b/a is the % of slow that is theta
        % the other way is to find the inflections in the power spectrum in the
        % ac, and to also calculate the peak power from 6-12 and divide that by
        % avearge across other freqs (andy did this)
        if verbose
            figure;
            slow=(royercf(1)*(cos(2*.001*pi*royercf(4)*x)+1)+royercf(2)).*exp(-abs(x)/royercf(5));
            fast=royercf(3)*exp(-(x.^2)/royercf(6)^2);
            hfig=figure;
            b=bar([-700:10:700]',acgram2,'FaceColor','flat','EdgeColor','flat');
            b.CData=[.7 .7 .7];
            hold on; plot(x,slow); plot(x,fast); plot(x,slow+fast);
            drawnow;
            answer = questdlg('Stop verbose?', 'Verbose Menu', 'Yes, quiet','no, more figs','Yes, quiet');
            verbose=contains(answer,'no');
        end
    end
    fprintf('Session %d done in %.2f seconds \n',ses,toc);
end
%}

%% now plot them to get a good impression of pyr vs IN
for ses=1:length(SuperRat)
    for i=1:length(SuperRat(ses).units)
        if isempty(SuperRat(ses).units(i).tag)
            SuperRat(ses).units(i).tag='mua';
        end
    end
end
AllCells=[];
for ses=1:length(SuperRat)
    units=SuperRat(ses).units;
    SesCells=[[units.meanrate]' [units.propbursts]' [units.burstIndex]'...
        cellfun(@(a) contains(a,'CA1'),{units.area})'...
        cellfun(@(a) contains(a,'accepted'),{SuperRat(ses).units.tag})'];
    AllCells=[AllCells; SesCells];
end



figure; 
subplot(1,2,1); scatter(log2(AllCells(AllCells(:,4)==1 & AllCells(:,5)==1,1)),...
    log2((AllCells(AllCells(:,4)==1 & AllCells(:,5)==1,3))),15,'filled');
hold on;
scatter(log2(AllCells(AllCells(:,4)==1 & AllCells(:,5)==0,1)),...
    log2((AllCells(AllCells(:,4)==1 & AllCells(:,5)==0,3))),15,'filled');
xlabel('Log_2 Firing Rate'); ylabel('Log_2 Burst Index');
legend('Accepted','MUA');
set(gca,'YLim',[0 10]); title('HPC units');


subplot(1,2,2); scatter(log2(AllCells(AllCells(:,4)==0 & AllCells(:,5)==1,1)),...
    log2((AllCells(AllCells(:,4)==0 & AllCells(:,5)==1,3))),15,'filled');
hold on;
scatter(log2(AllCells(AllCells(:,4)==0 & AllCells(:,5)==0,1)),...
    log2((AllCells(AllCells(:,4)==0 & AllCells(:,5)==0,3))),15,'filled');
xlabel('Log_2 Firing Rate'); ylabel('Log_2 Burst Index');
legend('Accepted','MUA');
set(gca,'YLim',[0 10]); title('PFC units');

% for hpc units it looks like INS have a rate above 3 Hz and a burst index 
% of below 3
% for PFC it looks like INS have a firing rate above 7 Hz



%% now the above plot but labeled correctly


keepfields={'meanrate','propbursts','burstIndex','area','tag','type'};
AllCells=[];
for ses=1:length(SuperRat)
    units=SuperRat(ses).units;
    allnames=fieldnames(units); extranames=~ismember(allnames,keepfields);
    units=rmfield(units,allnames(extranames==1));
    AllCells=[AllCells units];
end
HPCinds=cellfun(@(a) contains(a,'CA1'), {AllCells.area});
alltypes=cellfun(@(a) contains(a,'pyr'), {AllCells.type});
alltypes=alltypes+2*cellfun(@(a) contains(a,'in'), {AllCells.type});
% this plots log2 mean rate against log2 burst index
figure; sp=subplot(1,2,1); 
scatter(log2([AllCells(HPCinds & alltypes==1).meanrate]),...
    log2([AllCells(HPCinds & alltypes==1).burstIndex]),15,'filled');
hold on; scatter(log2([AllCells(HPCinds & alltypes==2).meanrate]),...
    log2([AllCells(HPCinds & alltypes==2).burstIndex]),15,'filled');
hold on; scatter(log2([AllCells(HPCinds & alltypes==0).meanrate]),...
    log2([AllCells(HPCinds & alltypes==0).burstIndex]),15,'filled');
xlabel('Log_2 Firing Rate'); ylabel('Log_2 Burst Index');

PFCinds=cellfun(@(a) contains(a,'PFC'), {AllCells.area});
sp(2)=subplot(1,2,2); 
scatter(log2([AllCells(PFCinds & alltypes==1).meanrate]),...
    log2([AllCells(PFCinds & alltypes==1).burstIndex]),15,'filled');
hold on; scatter(log2([AllCells(PFCinds & alltypes==2).meanrate]),...
    log2([AllCells(PFCinds & alltypes==2).burstIndex]),15,'filled');
hold on; scatter(log2([AllCells(PFCinds & alltypes==0).meanrate]),...
    log2([AllCells(PFCinds & alltypes==0).burstIndex]),15,'filled');
xlabel('Log_2 Firing Rate'); ylabel('Log_2 Burst Index');
linkaxes(sp);
%% now lets look at 'accepted' vs 'mua'


% this requires these cells to be fully processed.  This means spatial
% properties and object response properties

% first lets pull out our cells with a lower mean rate
% 1. mean rate 2. pburst 3. burst index 4. 1ms isi 5. 2ms isi 6. actheta
% 7. is ca1, 8. is accepted, 9. numspikes 10 odor responsive 11 odor
% selective 12. change at odor onset

AllCells=[];
for ses=1:length(SuperRat)
    units=SuperRat(ses).units;
    SesCells=[[units.meanrate]' [units.propbursts]' [units.burstIndex]'...
        cell2mat(cellfun(@(a) a, {units.shortISI}, 'UniformOutput', false)')...
        [units.ACtheta]' cellfun(@(a) contains(a,'CA1'),{units.area})'...
        cellfun(@(a) contains(a,'accepted'),{units.tag})' [units.numspikes]'...
        cellfun(@(a) a(4),{units.OdorResponsive})'...
        cellfun(@(a) a(3),{units.OdorSelective})'...
        cellfun(@(a) a(2),{units.OdorResponsive})'...
        ];
    AllCells=[AllCells; SesCells];
end

% start with CA1
AllHC=AllCells(AllCells(:,7)==1 & sum(isnan(AllCells),2)==0,:);

% mean rate
figure;
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==1,1)),50);
sp=subplot(7,2,1); bar(a(2:end)-(a(2)-a(1))/2,b); title('accepted'); xlabel('mean rate (Log_2)'); ylabel('hist bin counts');
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==0,1)),a);
sp(2)=subplot(7,2,2); bar(a(2:end)-(a(2)-a(1))/2,b); title('MUA'); xlabel('mean rate(Log_2)'); ylabel('hist bin counts');
linkaxes(sp,'x'); 


% burst prob
[b,a]=histcounts(AllHC(AllHC(:,8)==1,2),50);
sp=subplot(7,2,3); bar(a(2:end)-(a(2)-a(1))/2,b); title('accepted''s'); xlabel('Burst Prob'); ylabel('hist bin counts');
[b,a]=histcounts(AllHC(AllHC(:,8)==0,2),a);
sp(2)=subplot(7,2,4); bar(a(2:end)-(a(2)-a(1))/2,b); title('MUA'); xlabel('Burst Prob'); ylabel('hist bin counts');
linkaxes(sp,'x'); 

% 1ms IS %%
[b,a]=histcounts(AllHC(AllHC(:,8)==1,4)./AllHC(AllHC(:,8)==1,9),50);
sp=subplot(7,2,5); bar((a(2:end)-(a(2)-a(1))/2)*100,b); title('accepted'); xlabel('1ms ISI(%)'); ylabel('hist bin counts');
[b,a]=histcounts(AllHC(AllHC(:,8)==0,4)/AllHC(AllHC(:,8)==0,9),a);
sp(2)=subplot(7,2,6); bar((a(2:end)-(a(2)-a(1))/2)*100,b); title('MUA'); xlabel('1ms ISI(%)'); ylabel('hist bin counts');
linkaxes(sp,'x'); xlim([0 2]);


% 2ms ISI %%
[b,a]=histcounts(AllHC(AllHC(:,8)==1,5)./AllHC(AllHC(:,8)==1,9),50);
sp=subplot(7,2,7); bar((a(2:end)-(a(2)-a(1))/2)*100,b); title('accepted'); xlabel('2ms ISI(%)'); ylabel('hist bin counts');
[b,a]=histcounts(AllHC(AllHC(:,8)==0,5)/AllHC(AllHC(:,8)==0,9),a);
sp(2)=subplot(7,2,8); bar((a(2:end)-(a(2)-a(1))/2)*100,b); title('MUA'); xlabel('2ms ISI (%)'); ylabel('hist bin counts');
linkaxes(sp,'x'); xlim([0 2]);

% total spikes
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==1,9)),50);
sp=subplot(7,2,9); bar(a(2:end)-(a(2)-a(1))/2,b); title('accepted'); xlabel('total spikes (Log_2)'); ylabel('hist bin counts');
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==0,9)),a);
sp(2)=subplot(7,2,10); bar(a(2:end)-(a(2)-a(1))/2,b); title('MUA'); xlabel('total spikes(Log_2)'); ylabel('hist bin counts');
linkaxes(sp,'x');

% theta modulation index
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==1,6)),50);
sp=subplot(7,2,11); bar(a(2:end)-(a(2)-a(1))/2,b); title('accepted'); xlabel('Theta Mod index (Log_2)'); ylabel('hist bin counts');
[b,a]=histcounts(log2(AllHC(AllHC(:,8)==0,6)),a);
sp(2)=subplot(7,2,12); bar(a(2:end)-(a(2)-a(1))/2,b); title('MUA'); xlabel('Theta Mod index (Log_2)'); ylabel('hist bin counts');
linkaxes(sp,'x'); 

sgtitle('CA1 units');

%% odor responsive
figure; subplot(2,1,1);
% this is [acc- acc+; mua- mua+]'
bary=[nanmean(AllHC(AllHC(:,8)==1,10)<.05 & AllHC(AllHC(:,8)==1,12)<0)*100,...
    nanmean(AllHC(AllHC(:,8)==1,10)<.05 & AllHC(AllHC(:,8)==1,12)>0)*100;...
    nanmean(AllHC(AllHC(:,8)==0,10)<.05 & AllHC(AllHC(:,8)==0,12)<0)*100,...
    nanmean(AllHC(AllHC(:,8)==0,10)<.05 & AllHC(AllHC(:,8)==0,12)>0)*100];

hb=bar([1 2],bary,'stacked');    
set(gca,'XTickLabel',{sprintf('accepted (%d)',sum(AllHC(:,8)==1)),...
    sprintf('MUA (%d)',sum(AllHC(:,8)==0))});
xtips1 = hb(1).XEndPoints; ytips1 = hb(2).YEndPoints;
labels1 = string([sum(AllHC(AllHC(:,8)==1,10)<.05) sum(AllHC(AllHC(:,8)==0,10)<.05)]);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
ylabel('% of units that are odor responsive'); legend('Lowers at odor','elevates at odor');

% odor selectivity
subplot(2,1,2);
hb=bar([nanmean(AllHC(AllHC(:,8)==1,11))*100 nanmean(AllHC(AllHC(:,8)==0,11))*100]);    
set(gca,'XTickLabel',{sprintf('accepted (%d)',sum(AllHC(:,8)==1)),...
    sprintf('MUA (%d)',sum(AllHC(:,8)==0))});
xtips1 = hb(1).XEndPoints; ytips1 = hb(1).YEndPoints;
labels1 = string([sum(AllHC(AllHC(:,8)==1,11)) sum(AllHC(AllHC(:,8)==0,11))]);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel('% of units that are odor selective');

sgtitle('CA1 Units');
%%

% rate against burstiness
figure; 
subplot(2,2,1);
scatter(log2(AllHC(AllHC(:,8)==1,1)),AllHC(AllHC(:,8)==1,2),15,'filled');
hold on;
scatter(log2(AllHC(AllHC(:,8)==0,1)),AllHC(AllHC(:,8)==0,2),15,'filled');
xlabel('mean rate (log_2)'); ylabel('burst ratio'); legend('accepted','MUA');

% rate vs isis
subplot(2,2,2); scatter(log2(AllHC(AllHC(:,8)==1,1)),AllHC(AllHC(:,8)==1,4)./AllHC(AllHC(:,8)==1,9),15,'filled');
hold on;
scatter(log2(AllHC(AllHC(:,8)==0,1)),AllHC(AllHC(:,8)==0,4)./AllHC(AllHC(:,8)==0,9),15,'filled');
xlabel('mean rate (log_2)'); ylabel('ISI violations');

% rate vs theta mod
subplot(2,2,3); scatter(log2(AllHC(AllHC(:,8)==1,1)),log2(AllHC(AllHC(:,8)==1,6)),15,'filled');
hold on;
scatter(log2(AllHC(AllHC(:,8)==0,1)),log2(AllHC(AllHC(:,8)==0,6)),15,'filled');
xlabel('mean rate (log_2)'); ylabel('theta modulation'); 


% ac theta vs object responsiveness
subplot(2,2,4); scatter(log2(AllHC(AllHC(:,8)==1,10)).*(((AllHC(AllHC(:,8)==1,12)>0)-.5)*2),log2(AllHC(AllHC(:,8)==1,6)),15,'filled');
hold on;
scatter(log2(AllHC(AllHC(:,8)==0,10)).*(((AllHC(AllHC(:,8)==0,12)>0)-.5)*2),log2(AllHC(AllHC(:,8)==0,6)),15,'filled');
xlabel({'<- depress at odor                elevate at odor ->';'object responsivity (log_2P)'}); ylabel('theta modulation Index');

sgtitle('CA1 units');



