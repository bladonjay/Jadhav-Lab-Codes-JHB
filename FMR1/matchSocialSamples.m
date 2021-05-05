function [allevents,eachevents] = matchSocialSamples(myevents,hisevents)
% function [fulltable,eachrat] = matchSocialSamples(rat1events,rat2events)
%
%
%
% output format
% should probably be a table in this form

% rat1 timestamp
% rat1 entry/exit
% rat1 well ID

% rat2 timestamp
% rat2 entry/exit
% rat2 well ID

% whether it was rewarded
% whether it was a match


% the final thing i want is to take each animals arm transitions, see
% whether they are leaving from where the other naimal is, see whether the
% other animal is currently at that other well

% first concatenate rats
myevents.ratnum(:)=1; hisevents.ratnum(:)=2;
allevents=[myevents; hisevents];

tempevents1=allevents(:,[1 3:8]);
tempevents1.in1out0(:)=1;
tempevents2=allevents(:, [2 3:8]);
tempevents2.in1out0(:)=0;
tempevents2.Properties.VariableNames{1} = 'start';
allevents=[tempevents1; tempevents2];
allevents = sortrows(allevents,'start','ascend');
allevents.start=allevents.start-allevents.start(1);

eachevents{1}=allevents(allevents.ratnum==1,:);
eachevents{2}=allevents(allevents.ratnum==2,:);


end

