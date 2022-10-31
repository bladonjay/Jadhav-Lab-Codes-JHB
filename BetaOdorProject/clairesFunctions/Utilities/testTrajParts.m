function testTrajParts(animal, day, epoch)
%run inside animal_direct folder

%plots each stem part individually for each trail so you can check to make
%sure all make sense


daystr = getTwoDigitNumber(day);


load([animal, 'runTrajBounds', daystr,'.mat']);
load([animal, 'pos', daystr,'.mat']);


stemintime = runTrajBounds{1,day}{1,epoch}.data(:,9);
stemouttime = runTrajBounds{1,day}{1,epoch}.data(:,10);
posdata = pos{1,day}{1,epoch}.data;


figure,
set(gcf,'Position',[50 100 600 400])
hold on

xvals = posdata(:,2);
yvals = posdata(:,3);
postime = posdata(:,1);

for i = 1:length(stemintime)
stemin = (lookup(stemintime(i), postime));

stemout = (lookup(stemouttime(i), postime));

stemX = xvals(stemin:stemout);
stemY = yvals(stemin:stemout);

plot(stemX,stemY,'k-');

disp(i)
pause

end

