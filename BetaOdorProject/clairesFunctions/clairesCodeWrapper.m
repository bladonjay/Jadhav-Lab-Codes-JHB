%%clairesCodeWrapper
%
%
%
%
%
%
clairesPath=genpath(cd);
addpath(clairesPath);

% start with the pipeline
edit cs_Pipeline
%% first lets generate our taskResponsive cells and Selective cells

% these are the cells that fire during run epochs (at least 1 spike)
edit cs_listRunCells.m


% this produces 'npCells' which are all the nosepoke PYRAMIDAL cells
edit cs_listNPCells_v2
%1168 total pyrams in CA1
%170 active pyrams in CA1
%138 np pyrams in CA1
%664 total pyrams in PFC
%234 active pyrams in PFC
%185 np pyrams in PFC

edit cs_listNPcells


% then run odor selective
edit cs_cellSelectivityTag

% then all cells
edit cs_cellPopulations

%%

% this function creates a cellpops struct and lists the cell counts using
% her sankey filter: runCells, pyrCells, NP cells, Selective Cells
%{
the comments from shantanu are as follows
1. shantanu wants to use claires numbers and not change our cell counts
2. i need to find the number of neurons that fire at least as many spikes
as trials, and mention that this is a second filter for active cells, thus
the proportions of cells that are selective are really the proportion of
cells WHO SPIKE AT ALL that are selective and these numbers are in line
with prevfious reports

3. i need to recalculate the overlap between these cells and the coherence
results, especially for cross-region coherence.  I'll have to use claires
code because now we are using claires dataset instead of mine.