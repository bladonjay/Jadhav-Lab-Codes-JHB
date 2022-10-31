%cs_getTCtrace

%smooth thermocouple trace can be obtained by subtracting the signals from
%the two channels from one another. Do this for all days, and save in EEG
%file. 

animal = 'CS42';
topDir = cs_setPaths();

animDir = [topDir, animal,'Expt\',animal,'_direct\'];
tetinfo = loaddatastruct(animDir,animal,'tetinfo');
filt = 'strcmp($area,''TC'')';
tet = evaluatefilter(tetinfo,filt);
tet = unique(tet(:,3));

tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
filt = '~isempty($area)';
dayeps = evaluatefilter(tetinfo,filt);
dayeps = unique(dayeps(:,[1 2]),'rows');

%dayeps(find(dayeps(:,1)==4),:) = [];

for ep = 22:size(dayeps,1)
    day = dayeps(ep,1);
    ep = dayeps(ep,2);

    disp(['Doing day ',daystr,' ep ',epstr]);
daystr = getTwoDigitNumber(day);
epstr = getTwoDigitNumber(ep);
tetstr = getTwoDigitNumber(tet);

load([animDir,'EEG_TC\',animal,'eeg',daystr,'-',epstr,'-',tetstr,'ch01']);
ch1 = eeg{day}{ep}{tet}.data(:,1);
load([animDir,'EEG_TC\',animal,'eeg',daystr,'-',epstr,'-',tetstr,'ch02']);
ch2 = eeg{day}{ep}{tet}.data(:,1);

tcdata = ch1-ch2;

eeg{day}{ep}{tet}.nTrodeChannel = '1-2';
eeg{day}{ep}{tet}.data = tcdata;

%figure, plot(tcdata(1:1000))

save([animDir,'EEG\',animal,'eeg',daystr,'-',epstr,'-',tetstr],'eeg');
end
