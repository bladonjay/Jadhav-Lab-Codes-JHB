function [rat1events, rat2events]=parseSocialEvents(dataRaw,debounce)
% function [rat1events, rat2events, allevents]=parseSocialEvents(ledger)
% splits the stream of events by rat and tracks what they and their partner
% were doing for each sample/arm transition

% rows are as follows;
%{
1. timestamp in
2. timestamp out
3. well
4. last well
5. was it rewarded?
6. was it a match
7. if there was a reward arm assigned, this is its assignment
%}

% debounce is the lockout period for a second 'sample' once the first
% sample has happened.  Sometimes the IRbeams are finicky and will trip
% multiple times in a very short window.  This prevents that decently well

cellData = table2cell(dataRaw); % convert to cell array
numIdx = cellfun(@(x) ~isempty(x{1}), cellData(:,1)); % now look down first column
rawEvents = cellData(numIdx,1);

% convert to char and split out
DataTips2=cellfun(@(a) char(a{1}), rawEvents, 'UniformOutput',false);
DataAll=cellfun(@(a) a(1:find(a==' ',1,'first')), DataTips2,'UniformOutput',false);
DataAll(:,2)=cellfun(@(a) a(find(a==' ',1,'first')+1:end), DataTips2,'UniformOutput',false);

ledger=DataAll;


if ~exist('debounce','var')
    debounce=0;
end
% preallocate
rat1events=zeros(1,7); counter1=2;
rat2events=zeros(1,7); counter2=2;
% take only real events
eventlist=string(ledger);
oktouse=cellfun(@(a) any(a), isstrprop(ledger(:,2),'alpha'));
eventlist=eventlist(oktouse,:);
% go and assign each events to a sample with a duration
for i=1:length(eventlist)
    %  okay so very first event must be a poke
    % rows are 1 ts, 2 rat ID, rat 1 or rat 2
    if contains(eventlist(i,2),'poke','IgnoreCase',true) && ...
            ~contains(eventlist(i,2),'position','IgnoreCase',true)% the space is important
        charevent=char(eventlist(i,2));
        % now choose which rat it is
        if isstrprop(charevent(end),'digit') % if its a number, rat 1
            if ~contains(eventlist(i,2),'un','IgnoreCase',true)% only take poke or unpoke events... we know what follows
                rat1events(counter1,1) = str2double(eventlist(i,1))/1000; % recorded in msec
                rat1events(counter1,3)=str2double(charevent(end)); % well number
                rat1events(counter1,[4 7])=rat1events(counter1-1,[3 7]); % last well
                counter1=counter1+1; % advance event ledger
            else
                rat1events(counter1-1,2)=str2double(eventlist(i,1))/1000;
            end
        else
            if ~contains(eventlist(i,2),'un','IgnoreCase',true)% only take poke or unpoke events... we know what follows
                rat2events(counter2,1) = str2double(eventlist(i,1))/1000;
                rat2events(counter2,3)=double(charevent(end))-64;
                rat2events(counter2,[4 7])=rat2events(counter2-1,[3 7]);
                counter2=counter2+1;
            else
                rat2events(counter2-1,2)=str2double(eventlist(i,1))/1000;
            end
        end
    elseif contains(eventlist(i,2),'Matched','IgnoreCase',true) || ...
        contains(eventlist(i,2),'Mathced','IgnoreCase',true)% eitgher match and no reward or match and reward
        rat1events(end,6)=1; % register the previous sample as a match
        rat2events(end,6)=1;
        if contains(eventlist(i+1,2),'count') % if we add to the reward ct, this was rewarded
            rat1events(end,5)=1; % register the previous sample as rewarded
            rat2events(end,5)=1;
            rewardStr=char(eventlist(i+2,2));
            rat1events(end,7)=str2double(rewardStr(end));
            rat2events(end,7)=str2double(rewardStr(end));
        end % if it didnt, this wasnt rewarded
    end
end

% back fix to assign correct column well to change after correct sample not
% before
winids=find(rat1events(:,6)==1);
rat1events(winids,7)=rat1events(winids-1,7);
winids=find(rat2events(:,6)==1);
rat2events(winids,7)=rat2events(winids-1,7);

rat1events=rat1events(2:end,:); % remove the first event (that we used for last sample)
rat2events=rat2events(2:end,:);
%%%%%%%%%%%
% working with debouncing
%{
returnevents=rat1events(rat1events(:,3)==rat1events(:,4),:);
nonreturns=rat1events(rat1events(:,3)~=rat1events(:,4),:);
returnlags=diff(linearize(returnevents(:,1:2)'));
returnlags=returnlags(2:2:end);
nonreturnlags=diff(linearize(nonreturns(:,1:2)'));
nonreturnlags=nonreturnlags(2:2:end);
figure;
[y,x]=histcounts(log(returnlags),[-3:.2:4]);
[y2,x]=histcounts(log(nonreturnlags),[-3:.2:4]);
newx=mean([x(1:end-1); x(2:end)]);
figure; bar(newx,y,1,'FaceAlpha',.5,'EdgeColor','none');
hold on;
bar(newx,y2,1,'FaceAlpha',.5,'EdgeColor','none');
xlabel('seconds 10^x')
legend('return lags','arm transition lags');
%}
%%%%%%%%%%
% first, i notice that some samples are stuck, so we need to destick
rat1events(rat1events(:,2)==0,:)=[];
rat2events(rat2events(:,2)==0,:)=[];



% basically merge all consecutive events at the same port that are within a
% second of eachother
if debounce>0
    i=1;
    while i<length(rat1events)
        if rat1events(i,3)==rat1events(i+1,3) &&...
                (rat1events(i+1,1)-rat1events(i,2))<debounce
            
            rat1events(i,2)=rat1events(i+1,2); % take new end
            rat1events(i,5)=rat1events(i,5)+rat1events(i+1,5); % take if rewareded
            rat1events(i,6)=rat1events(i,6)+rat1events(i+1,6); % take if matchec
            rat1events(i+1,:)=[]; % delete next row           
        else
            i=i+1;
        end
    end
    i=1;
    while i<length(rat2events)
        if rat2events(i,3)==rat2events(i+1,3) &&...
                (rat2events(i+1,1)-rat2events(i,2))<debounce
            
            rat2events(i,2)=rat2events(i+1,2); % take new end
            rat2events(i,5)=rat2events(i,5)+rat2events(i+1,5); % take if rewareded
            rat2events(i,6)=rat2events(i,6)+rat2events(i+1,6); % take if matchec
            rat2events(i+1,:)=[]; % delete next row            
        else
            i=i+1;
        end
    end
end

end
