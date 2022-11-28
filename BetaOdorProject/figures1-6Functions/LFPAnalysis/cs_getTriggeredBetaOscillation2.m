function cs_getTriggeredBetaOscillation2(animal, day, epoch, trial)

topDir = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];
load([animDir, animal, 'tetinfo.mat']);

tetfilter = 'strcmp($descrip2, ''betatet'') & strcmp($area,''CA1'')';
tets(1) = evaluatefilter(tetinfo{day}{epoch},tetfilter);

tetfilter = 'strcmp($descrip2, ''betatet'') & strcmp($area,''PFC'')';
tets(2) = evaluatefilter(tetinfo{day}{epoch},tetfilter);

tetfilter = 'strcmp($descrip2, ''betatet'') & strcmp($area,''OB'')';
tets(3) = evaluatefilter(tetinfo{day}{epoch},tetfilter);


% tets = [ca1, pfc, ob];



close all
[topDir,figDir]= cs_setPaths();

colors = {rgb('MediumIndigo'), rgb('DarkAquamarine'), rgb('LightSkyBlue')};
dataDir = [topDir,animal,'Expt\',animal,'_direct\'];

daystr = getTwoDigitNumber(day);
epochstr = getTwoDigitNumber(epoch);

figure,
for t = 1:length(tets)
    tet = tets(t);
    
tetstr = getTwoDigitNumber(tet);

load([dataDir, 'EEG\',animal,'beta',daystr,'-',epochstr,'-',tetstr,'.mat']);
beta = beta{1, day}{1, epoch}{1, tet}  ;
load([dataDir, animal, 'odorTriggers', daystr, '.mat'])
odorTriggers = odorTriggers{1, day}{1, epoch}.correctTriggers  ;

trigger = odorTriggers(trial);

window = [0.2 1.2];
start = trigger-window(1);
stop = trigger+window(2);

timerange = beta.timerange;
betatime = [timerange(1):1/1500:timerange(2)];
betatime = betatime';

betaamp = beta.data(:,1);
betaamp = double(betaamp);

betatimes = (betatime >= start & betatime<= stop);
betainwin = betaamp(betatimes);

%stop = stop-trigger;
%start = start-trigger;

windowtimes = [-window(1):1/1500:window(2)];
windowtimes = windowtimes(1:end-1);



eval(['p',num2str(t),'= plot(windowtimes,betainwin,''Color'',colors{t},''LineWidth'',2);']);
eval(['p',num2str(t),'.Color(4) = 0.75;']);
hold on
plot([0 0], [-350 350], 'r--', 'LineWidth', 2);
axis([-window(1) window(2) -350 350]);
%plot([0 0], [-min(betainwin)-(min(betainwin)*0.1) max(betainwin)+(max(betainwin)*0.1)], 'r--','LineWidth',2);

set(gcf, 'Position', [300 500 800 350]);
                set(gca,'fontsize',20);
box off

xlabel('Time (seconds)');
ylabel('Beta amplitude (\muV)');
end

figtitle = ['triggeredBeta'];
    newfigfile = [figDir,'NicePPTFigures\',figtitle];
    saveas(gcf,newfigfile,'fig');
    print('-dpdf', newfigfile);
    print('-djpeg', newfigfile);
    
axis([0.4 0.6 -300 300])
set(gcf, 'Position', [250 500 800 350]);
figtitle = ['triggeredBeta_zoom'];
    newfigfile = [figDir,'NicePPTFigures\',figtitle];
    saveas(gcf,newfigfile,'fig');
    print('-dpdf', newfigfile);
    print('-djpeg', newfigfile);

