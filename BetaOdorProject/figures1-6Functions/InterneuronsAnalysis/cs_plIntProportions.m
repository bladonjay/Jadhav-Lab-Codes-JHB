%cs_plCellProportions

%plot bar graphs with percentages of NP cells that are phase locked to beta
%vs RR
close all

[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
figDir = [figDir,'Interneurons\'];
freqbands = {'beta','resp'};
eegregions = {'CA1','PFC'};
regions = {'CA1','PFC'};
figure(1)
hold on
for r = 1:length(regions)
    region = regions{r};
    
        freqband = 'beta';

        load([dataDir, 'npInt_',region,'.mat']);
        allplcells = [];
        for er = 1:length(eegregions)
            eegregion = eegregions{er};
            load([dataDir, 'PhaseLocking\plINs_',freqband,'_',region,'-',eegregion,'.mat']);
            allplcells = [allplcells; plcells];
        end
        plcells_beta = unique(allplcells,'rows');
        x_beta = size(plcells_beta,1);
        
        
        freqband = 'resp';

        allplcells = [];
        for er = 1:length(eegregions)
            eegregion = eegregions{er};
            load([dataDir, 'PhaseLocking\plINs_',freqband,'_',region,'-',eegregion,'.mat']);
            allplcells = [allplcells; plcells];
        end
        plcells_resp = unique(allplcells,'rows');
        x_resp = size(plcells_resp,1);
        
        n = size(npInt,1);
        
        [p,z] = ztestprop(x_beta,x_resp,n,n);
        p_beta = x_beta/n;
        p_resp = x_resp/n;
        figure(1)
        bar([r,r+0.25],[p_beta,p_resp])
        nl = newline;
        text(r,0.1,['z = ',num2str(z),nl,'p = ',num2str(p)]);
        ylim([0 1])

        
    %saveas(gcf,figfile,'fig');
end
set(gca,'XTick',[1 1.2 2 2.2],'XTickLabel',{'Beta','RR'});
xlabel(regions);





%%
figfile = [figDir, 'pl_proportions'];
        print('-djpeg', figfile);
    print('-dpdf', figfile);