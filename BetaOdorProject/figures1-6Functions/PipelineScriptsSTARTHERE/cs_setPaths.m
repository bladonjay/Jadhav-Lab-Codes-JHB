%Sets path to data folder depending on computer. Make sure there is only
%one "OdorPlaceAssociation" folder on computer.
function  [topDir, figDir] = cs_setPaths(driveLetter)
if ~exist('driveLetter','var')
    driveLetter='A';
end
topDir=[];
figDir=[];
for i=1:26 %nletters in alphabet
    if exist(sprintf('%s:\\Data\\OdorPlaceAssociation\\',char(driveLetter)),'dir') == 7
        topDir = sprintf('%s:\\Data\\OdorPlaceAssociation\\',char(driveLetter));
        figDir = sprintf('%s:\\Figures\\',char(driveLetter));
        return
    elseif exist(sprintf('%s:\\Brandeis datasets\\OdorPlaceAssociation\\',char(driveLetter)),'dir') == 7
        topDir = sprintf('%s:\\Brandeis datasets\\OdorPlaceAssociation\\',char(driveLetter));
        figDir = sprintf('%s:\\Brandeis datasets\\OdorPlaceAssociation\\Figures\\',char(driveLetter));
        return
    elseif exist(sprintf('%s:\\BrandeisData\\SymanskiData\\',char(driveLetter)),'dir') == 7
        topDir = sprintf('%s:\\BrandeisData\\SymanskiData\\',char(driveLetter));
        figDir = sprintf('%s:\\BrandeisData\\SymanskiData\\',char(driveLetter));
        return
    else
        driveLetter=driveLetter+1;
    end
end