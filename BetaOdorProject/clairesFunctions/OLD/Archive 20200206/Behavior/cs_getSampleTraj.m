function cs_getSampleTraj(animal, day, epoch)

[topDir, figDir] = cs_setPaths();

speedthresh = 10;

animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
daystr = getTwoDigitNumber(day);
load([animDir,animal,'pos',daystr,'.mat'])
load([animDir,animal,'odorTriggers',daystr,'.mat'])
triggers = odorTriggers{day}{epoch}.allTriggers;
postime = pos{day}{epoch}.data(:,1);
posdata = pos{day}{epoch}.data(:,[2 3]);
for t = 1:length(triggers)
    trigger = triggers(t);
    times = [trigger, (trigger+5)];
    position = posdata(find(postime > times(1) & postime <times(2)),:);
    figure, plot(position(:,1), position(:,2),'k-');
    pause
    close
end
