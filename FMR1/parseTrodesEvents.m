function [myevents,eventlist] = parseTrodesEvents(ledger)
%function [myevents,eventlist] = parseTrodesEvents(ledger)
% parses a n by 2 cell mat of characters coming from the trodes ledger
myevents=[]; counter=1;
eventlist=string(ledger);
oktouse=cellfun(@(a) any(a), isstrprop(ledger(:,2),'alpha'));
eventlist=eventlist(oktouse,:);

% so in the end we want a n by 4 vector, 

% 1 start ts    (break time) 
% 2 end ts      (unbreak time)
% 3 port        (1, 2, or 3)
% 4 rewarded    (1 or 0(


for i=1:length(eventlist)
    % if its a port 1 poke
    if contains(eventlist(i,2),'portin1 up','IgnoreCase',true)
        myevents(counter,1)=eventlist(i,1);
        myevents(counter,3)=1;
        % now find the up
        nextup=find(contains(eventlist(i:end,2),'portin1 up','IgnoreCase',true),1,'first')+i;
        myevents(counter,2)=eventlist(nextup,1);
        % now find whether its rewarded
        if any(contains(eventlist(i:nextup,2),'reward','IgnoreCase',true) &...
                ~contains(eventlist(i:nextup,2),'not','IgnoreCase',true))
            myevents(counter,4)=1;
        else
            myevents(counter,4)=0;
        end
        counter=counter+1;
    % if its a port 2 poke
    elseif contains(lower(eventlist(i,2)),'portin2 up','IgnoreCase',true)
        myevents(counter,1)=eventlist(i,1);
        myevents(counter,3)=2;
        % now find the up
        nextup=find(contains(eventlist(i:end,2),'portin2 up','IgnoreCase',true),1,'first')+i;
        myevents(counter,2)=eventlist(nextup,1);
        % now find reward
        if any(contains(eventlist(i:nextup,2),'reward','IgnoreCase',true) &...
                ~contains(eventlist(i:nextup,2),'not','IgnoreCase',true))
            myevents(counter,4)=1;
        else
            myevents(counter,4)=0;
        end
        counter=counter+1;
    % if its a port 3 poke
    elseif contains(lower(eventlist(i,2)),'portin3 up','IgnoreCase',true)
        myevents(counter,1)=eventlist(i,1);
        myevents(counter,3)=3;
        % now find the up
        nextup=find(contains(eventlist(i:end,2),'portin3 up','IgnoreCase',true),1,'first')+i;
        myevents(counter,2)=eventlist(nextup,1);
        % now find reward
        if any(contains(eventlist(i:nextup,2),'reward','IgnoreCase',true) &...
                ~contains(eventlist(i:nextup,2),'not','IgnoreCase',true))
            myevents(counter,4)=1;
        else
            myevents(counter,4)=0;
        end
        counter=counter+1;
    end
end
% and convert to time
myevents(:,1:2)=myevents(:,1:2)/1000;


end

