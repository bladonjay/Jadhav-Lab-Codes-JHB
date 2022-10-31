%Plot raw EEG with markers for odor trigger times

function cs_raweegwithtrigs(prefix, day, epoch, tet)

datadirect = ['E:\Data\OdorPlaceAssociation\' prefix '_direct\'];
eegdirect = ['E:\Data\OdorPlaceAssociation\' prefix '_direct\EEG\'];

if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

if (epoch<10)
    epstring = ['0',num2str(epoch)];
else
    epstring = num2str(epoch);
end

if (tet < 10)
    tetstring = ['0',num2str(tet)];
else
    tetstring = num2str(tet);
end

eegFile = load([eegdirect, prefix, 'eegref', daystring, '-', epstring, '-', tetstring, '.mat']);

eegdata = eegFile.eegref{1,day}{1,epoch}{1,tet}.data;
starttime = eegFile.eegref{1,day}{1,epoch}{1,tet}.starttime;
endtime = eegFile.eegref{1,day}{1,epoch}{1,tet}.endtime;
samprate = eegFile.eegref{1,day}{1,epoch}{1,tet}.samprate;




% eegFile = load([eegdirect, prefix, 'eeg', daystring, '-', epstring, '-', tetstring, '.mat']);
% 
% eegdata = eegFile.eeg{1,day}{1,epoch}{1,tet}.data;
% starttime = eegFile.eeg{1,day}{1,epoch}{1,tet}.starttime;
% endtime = eegFile.eeg{1,day}{1,epoch}{1,tet}.endtime;
% samprate = eegFile.eeg{1,day}{1,epoch}{1,tet}.samprate;

eegtimes = [starttime:(1/samprate):endtime]';

figure,
plot(eegtimes, eegdata)
hold on

triggerFile = load([datadirect, 'OdorTriggers', daystring, '-', epstring, '.mat']);
OdorTriggers_all = triggerFile.OdorTriggers_all;




for i = 1:length(OdorTriggers_all)
plot([OdorTriggers_all(i) OdorTriggers_all(i)], [-2000 2000], 'k-', 'LineWidth', 2)
end

end