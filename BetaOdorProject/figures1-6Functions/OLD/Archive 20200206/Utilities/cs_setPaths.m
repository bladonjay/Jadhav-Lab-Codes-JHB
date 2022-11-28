%Sets path to data folder depending on computer. Make sure there is only
%one "OdorPlaceAssociation" folder on computer.
function  [topDir, figDir] = cs_setPaths()

if exist('F:\Data\OdorPlaceAssociation\') == 7
    topDir = 'F:\Data\OdorPlaceAssociation\';
    figDir = 'F:\Figures\';
elseif exist('D:\OdorPlaceAssociation\') == 7
    topDir = 'D:\OdorPlaceAssociation\';
    figDir = 'D:\Figures\';
elseif exist('E:\Brandeis datasets\Claire Data','dir')==7
    topDir='E:\Brandeis datasets\Claire Data\';
    figDir='G:\Brandeis datasets\Claire Data\ClaireGeneratedFigures\';
end