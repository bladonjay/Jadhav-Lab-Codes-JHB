cs_fixFigure(['Spectrogram_CA1_low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','CA1')
cs_fixFigure(['Spectrogram_PFC_low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','PFC')
cs_fixFigure(['Spectrogram_OB_low_-1000-2000ms_tapers1-1_20msBins.fig'], 'specgram','region','OB')

cs_fixFigure(['Cohgram_CA1-OB_Tapers2-3.fig'], 'cohgram', 'win', [1 2])
cs_fixFigure(['Cohgram_PFC-OB_Tapers2-3.fig'], 'cohgram', 'win', [1 2])
cs_fixFigure(['Cohgram_CA1-PFC_Tapers2-3.fig'], 'cohgram', 'win', [1 2])


cs_plotPhaseLocking_specifiedCells('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\NicePPTFigures\', 'CA1','CA1',[4 3 1], rgb('MediumPurple'))
cs_plotPhaseLocking_specifiedCells('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\NicePPTFigures\', 'PFC','PFC',[2 23 2], rgb('MediumAquamarine'))

cs_plotPhaseLockingFractions_v2('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\',{'CA1','PFC','OB'})

cs_popPhaseLockAnalysis('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\',{'CA1','PFC','OB'}, 50, 1)

cs_plotRasterPSTH_specifiedCells('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\', 'CA1', [.5 1], 0, 1);
cs_plotRasterPSTH_specifiedCells('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\', 'PFC', [.5 1], 0, 1);

cs_plotSelectivity({'CA1','PFC'}, [.5 1], 'D:\Figures\')

cs_plotPlSelectivityOverlap('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\', 'D:\Figures\')