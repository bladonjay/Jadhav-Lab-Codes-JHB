%look across trial. plot each spike as color according to phase.

triggers = odorTriggers{1}{2}.allTriggers;
 
s = spikes{1, 1}{1, 2}{1, 20}{1, 8}.data(:,1);
 
tetfilter = ['strcmp($area,''',region,''') & strcmp($descrip2,''betatet'')']; 
tet = evaluatefilter(tetinfo{1}{2}, tetfilter);

%s = spikes{1}{2}{1}{2}.data(:,1);
t = geteegtimes(beta{1}{2}{29});
betaphase = double(beta{1}{2}{29}.data(:,2));% beta phase
sph = betaphase(lookup(s, t));        
c = linspace(-pi,pi);
figure,  hold on
allt = [];
allph = [];

wins = betaWindows{1}{2};
for tr = 1:length(triggers)
    trig = triggers(tr);
    win = [wins(tr,1), trig+1.5];
    if ~any(isnan(win))
    goodspikes = s(find(isExcluded(s, win)));
    ph = sph(find(isExcluded(s,win)));
    ph = double(ph)/10000;
    above = find(ph>pi);
    if ~isempty(above)
        ph(above) = pi;
    end
   
    below = find(ph<-pi);
     if ~isempty(below)
         ph(below) = -pi;
     end
    
    spiketimes = goodspikes - win(1);
    allt = [allt;spiketimes];
    allph = [allph; ph];
%     colors = c(lookup(ph,c));
%     scatter(spiketimes, zeros(length(spiketimes),1)+tr, [], colors);
    end
end

[sortedT, i] = sort(allt);
sortedph = allph(i);


% if mod(length(i),2)
%     sortedT = sortedT(1:end-1);
%     i = i(1:end-1);
% end

earlyphase = allph(i(1:length(i)/2));

latephase = allph(i(length(i)/2+1:end));
figure
scatter(sortedT,sortedph,'r.');

[R,p] = corrcoef(allt, allph);

% earlyind = find(allt <= 0.2);
% earlyphase = allph(earlyind);
% mnearly = mean(earlyphase);
% 
% 
% lateind = find(allt > 0.5);
% latephase = allph(lateind);
% mnlate = mean(latephase);
%[~,p] = ttest2(earlyphase, latephase)

