%cs_plVenn
clear
cellregions = {'CA1','PFC'};
betaregions = {'CA1','PFC','OB'};
[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\PhaseLocking\'];
for r = 1:length(cellregions)
    region = cellregions{r};
    
    load([dataDir, 'plCells_beta_',region,'-CA1']);
    CA1 = plcells;
    
    load([dataDir, 'plCells_beta_',region,'-PFC']);
    PFC = plcells;
    
    load([dataDir, 'plCells_beta_',region,'-OB']);
    OB = plcells;
    
    CA1_PFC_OB = intersect(intersect(CA1,PFC,'rows'),OB,'rows');
    
    %remove cells that intersect all three
%     CA1 = CA1(~ismember(CA1, CA1_PFC_OB,'rows'),:);
%     PFC = PFC(~ismember(PFC, CA1_PFC_OB,'rows'),:);
%     OB = OB(~ismember(OB, CA1_PFC_OB,'rows'),:);
    
    CA1_PFC = intersect(CA1,PFC,'rows');
    CA1_OB = intersect(CA1,OB,'rows');
    PFC_OB = intersect(PFC,OB,'rows');
    
%     allintersect = [CA1_PFC_OB;CA1_PFC;CA1_OB;PFC_OB];
%     
%     CA1 = CA1(~ismember(CA1, allintersect,'rows'),:);
%     PFC = PFC(~ismember(PFC, allintersect,'rows'),:);
%     OB = OB(~ismember(OB, allintersect,'rows'),:);
    
    venn ([size(CA1,1), size(PFC,1), size(OB,1)], [size(CA1_PFC,1),size(CA1_OB,1), size(PFC_OB,1),size(CA1_PFC_OB,1)],'ErrMinMode', 'ChowRodgers')
    legend({'CA1','PFC','OB'})
    title([region, 'cells'])
    
end