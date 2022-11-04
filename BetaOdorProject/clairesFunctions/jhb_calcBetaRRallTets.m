% jhb_calcBetaRRallTets


topDir = cs_setPaths();
filtPath='C:\Users\Jadhavlab\Documents\gitRepos\LFP-analysis\JadhavEEGFilter\';
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};


for animNo=1:length(animals)
    animID=animals{animNo};
    
    dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory
    dayepochmatrix = cs_getRunEpochs(dataDir, animal, 'odorplace');
    days = unique(dayepochmatrix(:,1));
    for d=1:length(days)
        
        jhb_LFPdayprocess(dataDir, dataDir, animID, d, 'f', [filtPath 'thetafilter.mat'],'band','theta');
        disp(['doing day ',num2str(d), ' theta']);
        jhb_LFPdayprocess(dataDir, dataDir, animID, d, 'f', [filtPath 'betafilter.mat'],'band','beta');
        disp(['doing day ',num2str(d), ' beta']);
    end
end