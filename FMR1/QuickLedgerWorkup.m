%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 9);

% Specify sheet and range
opts.Sheet = "XF201-204";
opts.DataRange = "A122:I239";

% Specify column names and types
opts.VariableNames = ["Date", "TimeStart", "TimeEnd", "TaskParadigm", "AnimalPos112", "AnimalPos2ab", "Filename", "Notes","rates"];
opts.VariableTypes = ["datetime", "double", "double", "categorical", "double", "double", "string", "double","double"];

% Specify variable properties
opts = setvaropts(opts, "Filename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["TaskParadigm", "Filename"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Date", "InputFormat", "");

% Import the datas

dataloc=uigetfile();
dataloc="E:/Dropbox (Personal)/Brandeis/Social Task Recording Ledger 2021-1-29.xlsx"
SocialTaskLedger = readtable(dataloc, opts, "UseExcel", false);



%% parsing social task performance 



pitor=1;
vitor=1;
ipairnames={};
pairstat={}; % all pairings
individstat={}; % all individual animals
individ={'201','202','203','204'};
pairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
paircats=[1 1 2 3 1 1];
figure;
% so question one, how do each of our four guys stack up to eachother
for i=1:4
    icontains=cellfun(@(a) a==str2double(individ{i}), table2cell(SocialTaskLedger(:,[5 6])),'UniformOutput',false);
    individstat{i}=SocialTaskLedger(sum(cell2mat(icontains),2)>=1,:);
    plot(SmoothMat2(table2array(individstat{i}(:,9)),[10,10],2)); hold on;
    fprintf('%s on average = %.2f \n',individ{i},nanmean(table2array(individstat{i}(:,9))));
end
legend(individ);
% now lets see how each pair did
figure;
% so question one, how do each of our four guys stack up to eachother
for i=1:6
    icontains=cellfun(@(a) a==str2double(individ{pairs(i,1)}), table2cell(SocialTaskLedger(:,[5 6])),'UniformOutput',false);
    icontains=[icontains cellfun(@(a) a==str2double(individ{pairs(i,2)}), table2cell(SocialTaskLedger(:,[5 6])),'UniformOutput',false)];
    
    pairstat{i}=SocialTaskLedger(sum(cell2mat(icontains),2)>1,:);
    plot(SmoothMat2(table2array(pairstat{i}(:,9)),[10,10],2)); hold on;
    ipairnames{i}=num2str([individ{pairs(i,1)}, individ{pairs(i,2)}]);
    fprintf('%s on average = %.2f \n',ipairnames{i},...
        nanmean(table2array(pairstat{i}(:,9))));
end
legend(ipairnames);
% now we can compare
ctrls=cell2mat(cellfun(@(a) a.rates, pairstat(paircats==2), 'UniformOutput', false));

mixes=cell2mat(cellfun(@(a) a.rates, pairstat(paircats==1), 'UniformOutput', false));
fxes=cell2mat(cellfun(@(a) a.rates, pairstat(paircats==3), 'UniformOutput', false));
signrank(mean(mixes,2,'omitnan'),ctrls)

signrank(mean(mixes,2,'omitnan'),fxes)
figure; subplot(1,2,1);
plot([ctrls fxes]'); colororder(gca,parula(18));
hold on; boxplot([ctrls fxes]); set(gca,'XTickLabel',{'CTRL','FX'});
ylabel('Rewards per minute');
%errorbar(1,nanmean(ctrls), nanstd(ctrls),'kx','CapSize',0);
%errorbar(2,nanmean(fxes), nanstd(fxes),'kx','CapSize',0); xlim([.6 2.4]);
subplot(1,2,2);
plot(SmoothMat2([ctrls mean(mixes,2,'omitnan') fxes],[1,10],1),'*-');
legend('Controls','Mixes','FX');
ylabel('Rewards per minute');





%%
for i=1:height(SocialTaskLedger)
    if ~isnan(SocialTaskLedger.TimeStart(i))
    mypair(1)=find(contains(individ,num2str(SocialTaskLedger.AnimalPos112(i))));
    mypair(2)=find(contains(individ,num2str(SocialTaskLedger.AnimalPos2ab(i))));
    % now number this pair
    pairmatch=find(sum(pairs==mypair,2)==2);
    if isempty(pairmatch), pairmatch=find(sum(fliplr(pairs)==mypair,2)==2); end
    SocialTaskLedger.Pair(i)=pairmatch;
    
    for j=1:4
        SocialTaskLedger.(['has' num2str(j)])(i)=contains(individ(j),num2str(SocialTaskLedger.AnimalPos112(i))) ||...
            contains(individ(j),num2str(SocialTaskLedger.AnimalPos2ab(i)));
    end
    SocialTaskLedger.trial(i)=sum(SocialTaskLedger.Date(1:i)==SocialTaskLedger.Date(i));
    end
end
    
SocialTaskLedger(isnan(SocialTaskLedger.TimeStart),:)=[];
    
% need to do a signrank test for all comparisons...
pairp=[];
for i=1:length(pairstat)-1
    for j=i+1:length(pairstat)
        [a,b]=signrank(pairstat{i}.rates,pairstat{j}.rates);
        fprintf('%s (%.2f) vs %s (%.2f) pairwise p=%.4f \n',...
            [individ{pairs(i,1)} '-' individ{pairs(i,2)}],...
            nanmean(pairstat{i}.rates),...
            [individ{pairs(j,1)} '-' individ{pairs(j,2)}],...
            nanmean(pairstat{j}.rates), a);
        pairp(i,j)=a;
    end
end