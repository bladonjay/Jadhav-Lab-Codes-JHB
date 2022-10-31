function cs_getTriggeredResp(animal, day, epoch, trial)
%good one: CS33, day 3 ep 2 trial 23
[topDir,figDir] = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];
load([animDir, animal, 'tetinfo.mat']);

tets(1) = cs_getMostCellsTet(animal,day,epoch,'CA1');
tets(2) = cs_getMostCellsTet(animal,day,epoch,'PFC');
tets(3) = cs_getMostCellsTet(animal,day,epoch,'OB');
tets(4) = cs_getMostCellsTet(animal,day,epoch,'TC');

close all

colors = {rgb('MediumIndigo'), rgb('DarkAquamarine'), rgb('LightSkyBlue'), rgb('Tomato')};

daystr = getTwoDigitNumber(day);
epochstr = getTwoDigitNumber(epoch);

odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
trigs = odorTriggers{day}{epoch}.allTriggers;

figure,
for t = 1:3
    tet = tets(t);
    
tetstr = getTwoDigitNumber(tet);
load([animDir, 'EEG\',animal,'resp',daystr,'-',epochstr,'-',tetstr,'.mat']);
resp = resp{1, day}{1, epoch}{1, tet}  ;

trigger = trigs(trial);

window = [0.5 1];
start = trigger-window(1);
stop = trigger+window(2);

timerange = resp.timerange;
resptime = [timerange(1):1/1500:timerange(2)];
resptime = resptime';

respamp = resp.data(:,1);
respamp = double(respamp);

resptimes = (resptime >= start & resptime<= stop);
respinwin = respamp(resptimes);

%stop = stop-trigger;
%start = start-trigger;

windowtimes = [-window(1):1/1500:window(2)];
windowtimes = windowtimes(1:end-1);



eval(['p',num2str(t),'= plot(windowtimes,respinwin,''Color'',colors{t},''LineWidth'',2);']);
eval(['p',num2str(t),'.Color(4) = 0.75;']);
hold on
plot([0 0], [-350 350], 'r--', 'LineWidth', 2);
axis([-window(1) window(2) -350 350]);
%plot([0 0], [-min(betainwin)-(min(betainwin)*0.1) max(betainwin)+(max(betainwin)*0.1)], 'r--','LineWidth',2);

set(gcf, 'Position', [300 500 800 350]);
                set(gca,'fontsize',20);
box off

xlabel('Time (seconds)');
ylabel('Resp rhythm amplitude (\muV)');
end

t = 4;
tet = tets(4);
tetstr = getTwoDigitNumber(tet);
load([animDir, 'EEG\',animal,'eeg',daystr,'-',epochstr,'-',tetstr,'.mat']);
resp = eeg{1, day}{1, epoch}{1, tet}  ;
trigger = trigs(trial);

window = [0.5 1];
start = trigger-window(1);
stop = trigger+window(2);

timerange = resp.timerange;
resptime = [timerange(1):1/1500:timerange(2)];
resptime = resptime';

respamp = resp.data(:,1);
respamp = double(respamp);

resptimes = (resptime >= start & resptime<= stop);
respinwin = respamp(resptimes);

%stop = stop-trigger;
%start = start-trigger;

windowtimes = [-window(1):1/1500:window(2)];
windowtimes = windowtimes(1:end-1);



eval(['p',num2str(t),'= plot(windowtimes,respinwin,''Color'',colors{t},''LineWidth'',2);']);
eval(['p',num2str(t),'.Color(4) = 0.75;']);
hold on
plot([0 0], [-350 350], 'r--', 'LineWidth', 2);
axis([-window(1) window(2) -350 350]);
%plot([0 0], [-min(betainwin)-(min(betainwin)*0.1) max(betainwin)+(max(betainwin)*0.1)], 'r--','LineWidth',2);

set(gcf, 'Position', [300 500 800 350]);
                set(gca,'fontsize',20);
box off

xlabel('Time (seconds)');
ylabel('Resp rhythm amplitude (\muV)');


% figtitle = ['triggeredBeta'];
%     newfigfile = [figDir,'NicePPTFigures\',figtitle];
%     saveas(gcf,newfigfile,'fig');
%     print('-dpdf', newfigfile);
%     print('-djpeg', newfigfile);
%     %pause
%     
% axis([0.4 0.6 -300 300])
% set(gcf, 'Position', [250 500 800 350]);
% figtitle = ['triggeredBeta_zoom'];
%     newfigfile = [figDir,'NicePPTFigures\',figtitle];
%     saveas(gcf,newfigfile,'fig');
%     print('-dpdf', newfigfile);
%     print('-djpeg', newfigfile);
% 
% axis([-0.5 -0.3 -300 300])
% figtitle = ['triggeredBeta_zoom_beforeodor'];
%     newfigfile = [figDir,'NicePPTFigures\',figtitle];
%     print('-dpdf', newfigfile);
%     print('-djpeg', newfigfile);

