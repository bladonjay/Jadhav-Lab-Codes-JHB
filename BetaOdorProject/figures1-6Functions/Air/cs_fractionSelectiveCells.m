%tabulate number of np cells and selective cells in novel (pre/post
%learning) and air trials
%also plot figure showing numbers visually
close all
[topDir, figDir] = cs_setPaths();
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region= regions{r};
h = figure,
hold on

%novel_all
load([topDir,'AnalysesAcrossAnimals\npCells_novel_',region]);
np_novel = size(npCells,1);

load([topDir,'AnalysesAcrossAnimals\selectiveCells_novel_',region]);
sel_novel = size(selectivecells,1);

fract_novel = sel_novel/np_novel;


%novel_prelearn
load([topDir,'AnalysesAcrossAnimals\npCells_novel_prelearn_',region]);
    np_pre = size(npCells,1);

load([topDir,'AnalysesAcrossAnimals\selectiveCells_novel_prelearn_',region]);
    sel_pre = size(selectivecells,1);
    
fract_pre = sel_pre/np_pre;


%novel_postlearn
load([topDir,'AnalysesAcrossAnimals\npCells_novel_postlearn_',region]);
    np_post = size(npCells,1);

load([topDir,'AnalysesAcrossAnimals\selectiveCells_novel_postlearn_',region]);
    sel_post = size(selectivecells,1);
    
fract_post = sel_post/np_post;


%air
load([topDir,'AnalysesAcrossAnimals\npCells_air_',region]);
np_air = size(npCells,1);

load([topDir,'AnalysesAcrossAnimals\selectiveCells_air_',region]);
    sel_air = size(selectivecells,1);

fract_air = sel_air/np_air;

ytxt = {[num2str(sel_novel),'/',num2str(np_novel)], [num2str(sel_pre),'/',num2str(np_pre)], [num2str(sel_post),'/',num2str(np_post)], [num2str(sel_air),'/',num2str(np_air)]}; 
    

b = bar([1:4],[fract_novel,fract_pre,fract_post,fract_air]);
xticklabels({'','novel all','novel early learning','novel late learning','air'})
xtickangle(45)
ylabel('Fraction of odor-responsive cells'); 
%set(0,'defaultaxesfontsize',12);
ylim([0 1])
xlim([0 4.5])
set(gcf,'Position',[744 630 560 420]);

hold on

hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles

  hT=[hT text(b(1).XData+b(1).XOffset,b(1).YData,ytxt, ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
                      

figfile = [figDir, 'Spiking\fractSel_NovelAir_',region];
        print('-dpdf',figfile);
        print('-djpeg',figfile);
        
end

