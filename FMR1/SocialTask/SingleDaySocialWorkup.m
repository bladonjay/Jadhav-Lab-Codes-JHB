%%Motivation Workup

% the idea here is to get an impression for how well a set of animals is
% doing on a single day


myfolder=uigetdir();

allsess=dir(myfolder);
allsess=allsess(cellfun(@(a) contains(a,'stateScriptLog'), {allsess.name}));
% sort by recording time
[~,index] = sortrows({allsess.date}.'); allsess = allsess(index); clear index
% import options
[opts]=setStateScriptOpts();

% set the rat run schedule;
runRats=[nan nan; 1 3; 2 1; 3 4; 4 1; 2 3; 4 2; 2 1; 1 3; 3 1; 3 4;...
    nan nan; 2 3; 4 2];

for i=1:length(allsess); allsess(i).runRats=runRats(i,:); end
allsess(cellfun(@(a) any(isnan(a)), {allsess.runRats}))=[];
% parse the names

for i=1:length(allsess)
    DataFile = readtable(fullfile(allsess(i).folder,allsess(i).name), opts);
    DataRaw = table2cell(DataFile);
    % clean out empty rows
    numIdx = cellfun(@(x) ~isempty(x{1}), DataRaw(:,1)); DataTips = DataRaw(numIdx,1);
    
    % convert to char and split out
    DataTips2=cellfun(@(a) char(a{1}), DataTips, 'UniformOutput',false);
    ledger=cellfun(@(a) a(1:find(a==' ',1,'first')), DataTips2,'UniformOutput',false);
    ledger(:,2)=cellfun(@(a) a(find(a==' ',1,'first')+1:end), DataTips2,'UniformOutput',false);
    myevents={};
    debounce=1; % 1 second debounce... basically thats well too fast for an animal to transition arms
    [myevents{1},myevents{2}] = parseSocialEvents(ledger,debounce);
    % first turn each into a table
    for rt=1:2
    allsess(i).ratsamples{rt}=table(myevents{rt}(:,1),myevents{rt}(:,2),...
            myevents{rt}(:,3),myevents{rt}(:,4),myevents{rt}(:,5),...
            myevents{rt}(:,6), myevents{rt}(:,7),'VariableNames',...
            {'start','end','thiswell','lastwell','Reward','match','Goal Well'});
    end
    % now get where the partner has been
    for tr=1:height(allsess(i).ratsamples{1})
        hislastpoke=find(allsess(i).ratsamples{2}.start<allsess(i).ratsamples{1}.start(tr),1,'last');
        if ~isempty(hislastpoke)
            allsess(i).ratsamples{1}.hiswell(tr)=allsess(i).ratsamples{2}.thiswell(hislastpoke);
        else
            allsess(i).ratsamples{1}.hiswell(tr)=nan;
        end
    end
    
    for tr=1:height(allsess(i).ratsamples{2})
        hislastpoke=find(allsess(i).ratsamples{1}.start<allsess(i).ratsamples{2}.start(tr),1,'last');
        if ~isempty(hislastpoke)
            allsess(i).ratsamples{2}.hiswell(tr)=allsess(i).ratsamples{1}.thiswell(hislastpoke);
        else
            allsess(i).ratsamples{2}.hiswell(tr)=nan;
        end
    end
    % now get efficiency
    for rt=1:2
        myevents=allsess(i).ratsamples{rt}; % 1 then 2
        hisevents=allsess(i).ratsamples{3-rt}; % 2 then 1
        armtrans=myevents(diff(myevents.thiswell)~=0,:);
        allsess(i).RunCounts(rt)=height(armtrans);
        wins=armtrans.thiswell==armtrans.hiswell;
        allsess(i).Efficiency(rt)=nanmean(wins);
    end
    allsess(i).sessDur=max(cellfun(@(a) a.end(end), allsess(i).ratsamples))-...
        min(cellfun(@(a) a.start(1), allsess(i).ratsamples));
end

%%
% now lets aggregate per rat how they did
motivationAll=[];
for i=1:4
    tally=1;
    for ses=1:length(allsess)
        ratpos=find(allsess(ses).runRats==i);
        if any(ratpos)
            motivationAll(tally,i)=allsess(ses).RunCounts(ratpos)/allsess(ses).sessDur;
            tally=tally+1;
        end
    end
end
figure;
boxplot(motivationAll);
figure;
bar(nanmean(motivationAll));
hold on;
errorbar(1:4,nanmean(motivationAll),SEM(motivationAll,1),'.');
figure;
[a,b,c]=anovan(motivationAll(:),linearize(repmat(1:4,size(motivationAll,1),1)))
[a,b]=multcompare(c);
