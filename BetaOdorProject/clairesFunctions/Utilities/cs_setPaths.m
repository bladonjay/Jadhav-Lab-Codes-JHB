%Sets path to data folder depending on computer. Make sure there is only
%one "OdorPlaceAssociation" folder on computer.
function  [topDir, figDir] = cs_setPaths()

%topDir = 'D:\OdorPlaceAssociation\';
%    figDir = 'D:\Figures\';

if exist('F:\Data\OdorPlaceAssociation\') == 7
    topDir = 'F:\Data\OdorPlaceAssociation\';
    figDir = 'F:\Figures\';
elseif exist('D:\OdorPlaceAssociation\') == 7
    topDir = 'D:\OdorPlaceAssociation\';
    figDir = 'D:\Figures\';
end