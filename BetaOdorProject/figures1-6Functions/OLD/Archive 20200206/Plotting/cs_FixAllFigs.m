%trigtypes = {'allTriggers'};
trigtype = '';
 
for t = 1:length(trigtypes)
    trigtype = trigtypes{t};
cs_fixFigure(['Spectrogram_CA1_',trigtype,'low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','CA1')
cs_fixFigure(['Spectrogram_PFC_',trigtype,'low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','PFC')
cs_fixFigure(['Spectrogram_OB_',trigtype,'low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','OB')
end
pause
close all


for t = 1:length(trigtypes)
    trigtype = trigtypes{t};
    cs_fixFigure(['Cohgram_CA1-OB_',trigtype,'Tapers2-3.fig'], 'cohgram', 'win', [1 2])
    cs_fixFigure(['Cohgram_PFC-OB_',trigtype,'Tapers2-3.fig'], 'cohgram', 'win', [1 2])
    cs_fixFigure(['Cohgram_CA1-PFC_',trigtype,'Tapers2-3.fig'], 'cohgram', 'win', [1 2])
end
pause
close all

cs_plotPowerBins(['betaPowerOB_binned_-1000-1000ms.mat'])
cs_plotPowerBins(['betaPowerCA1_binned_-1000-1000ms.mat'])
cs_plotPowerBins(['betaPowerPFC_binned_-1000-1000ms.mat'])
close all

cs_plotPowerBins(['betaPowerOB_binned_-600-600ms.mat'])
cs_plotPowerBins(['betaPowerCA1_binned_-600-600ms.mat'])
cs_plotPowerBins(['betaPowerPFC_binned_-600-600ms.mat'])
close all

cs_plotPhaseLocking_specifiedCells('betaMod_CA1.mat','allTriggers',[1, 32, 6; 3, 21, 10; 3, 20, 6; 4, 29, 7], [0.8594 0.0781 0.2344])
cs_plotPhaseLocking_specifiedCells('betaMod_OB.mat','allTriggers',[3, 8, 2; 4, 19, 14; 4, 20, 8; 4, 28, 7], [0.1328 0.5430 0.1328])
cs_plotPhaseLocking_specifiedCells('betaMod_PFC.mat','allTriggers',[1, 27, 2; 2, 7, 1; 3, 18, 4; 3, 10, 1], [0.1172 0.5625 1])

cs_plotPhaseLocking_specifiedCells('betaMod_CA1.mat','allTriggers',[3, 21, 10], [0.1172 0.5625 1])
close all



%cs_popPhaseLockAnalysis('E:\AnalysesAcrossAnimals\', 'E:\', {'CA1','PFC'}, {'CA1','PFC','OB'}, 50, 1)
cs_popPhaseLockAnalysis('E:\AnalysesAcrossAnimals\', 'E:\', {'CA1','PFC'}, {'CA1'}, 50, 1)
close all


cs_plotPhaseLockingFractions('betaPhaseLockingAll.mat')

close all
cs_plotRasterPSTH_specifiedCells('TriggeredSpiking_CA1_-500-1000ms_100msBins.mat', [3,18,7;4,3,2;4,17,7;4,28,4;4,28,8],0)
close all
cs_plotRasterPSTH_specifiedCells('TriggeredSpiking_PFC_-500-1000ms.mat', [3,14,1;3,21,6;3,10,3;4,21,1;4,26,5],1)
close all

cs_plotSelectivity({'CA1','PFC'}, 'E:\Figures\')
close all

cs_plotPlSelectivityOverlap('E:\AnalysesAcrossAnimals\','E:\Figures\', {'CA1','PFC'})
close all

cs_plotEnsembleTraj('E:\AnalysesAcrossAnimals\', 'CS34', 'CA1', 4, {'leftTriggers','rightTriggers'}, 0)
cs_plotEnsembleTraj('E:\AnalysesAcrossAnimals\', 'CS34', 'PFC', 4, {'leftTriggers','rightTriggers'}, 0)
close all

cs_plotEnsembleDistance('ensembleDistance_CA1_-0.2-0.6ms.mat')
cs_plotEnsembleDistance('ensembleDistance_PFC_-0.2-0.6ms.mat')
cs_plotEnsembleDistance('ensembleDistance2_CA1_-200-600ms_50msBins.mat');
cs_plotEnsembleDistance('ensembleDistance3_CA1_-200-600ms_50msBins.mat');

cs_stateSpacePerf_plotAllAnimals({'CS31','CS33','CS34','CS35'}, {[1,2,3],[1,2,3,4],[1,2,3,4],[1,2,3]})
%cs_stateSpacePerf_plotAllAnimals({'SG8'}, {[1,2,3,4,5,6]})



 cs_cellSelectivityBins('TriggeredSpiking_CA1_-100-500ms_100msBins.mat');
 cs_cellSelectivityBins('TriggeredSpiking_PFC_-100-500ms_100msBins.mat');
 
 cs_plotEnsembleTraj(dataDir, animal, region, day, trigtypes, plotall)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_plotSpecgram('eventTrigSpecgramData_CA1_-1500-1500ms_tapers3_5.mat')
cs_plotSpecgram('eventTrigSpecgramData_CA1_-1500-1500ms_tapers2_3.mat')
cs_plotSpecgram('eventTrigSpecgramData_CA1_-1500-1500ms_tapers1_1.mat')

cs_plotSpecgram('eventTrigSpecgramData_PFC_-1500-1500ms_tapers3_5.mat')
cs_plotSpecgram('eventTrigSpecgramData_PFC_-1500-1500ms_tapers2_3.mat')
cs_plotSpecgram('eventTrigSpecgramData_PFC_-1500-1500ms_tapers1_1.mat')

cs_plotSpecgram('eventTrigSpecgramData_OB_-1500-1500ms_tapers3_5_30ms_bins.mat')

cs_plotSpecgram('eventTrigSpecgramData_OB_-1500-1500ms_tapers1_1_20ms_bins.mat')
cs_plotSpecgram('eventTrigSpecgramData_PFC_-1500-1500ms_tapers1_1_20ms_bins.mat')
cs_plotSpecgram('eventTrigSpecgramData_CA1_-1500-1500ms_tapers1_1_20ms_bins.mat')

cs_plotSpecgram('eventTrigSpecgramData_OB_-1000-1000ms_tapers1_1_20ms_bins.mat')
cs_plotSpecgram('eventTrigSpecgramData_PFC_-1000-1000ms_tapers1_1_20ms_bins.mat')
cs_plotSpecgram('eventTrigSpecgramData_CA1_-1000-1000ms_tapers1_1_20ms_bins.mat')

cs_stateSpacePerf('CS31', [4 5])
cs_stateSpacePerf('CS33', [1 2 3 4])
cs_stateSpacePerf('CS34', [1 2 3 4])
cs_stateSpacePerf('JS11', [1 2 3 4])


