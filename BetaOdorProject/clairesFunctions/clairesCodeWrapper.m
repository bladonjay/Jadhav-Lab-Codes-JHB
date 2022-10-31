%%clairesCodeWrapper
%
%
%
%
%
%


%% first lets generate our taskResponsive cells and Selective cells



cs_listNPCells_v2


% then run odor selective
cs_cellSelectivityTag_v2


%{
the comments from shantanu are as follows
1. shantanu wants to use claires numbers and not change our cel counts
2. i need to find the number of neurons that fire at least as many spikes
as trials, and mention that this is a second filter for active cells, thus
the proportions of cells that are selective are really the proportion of
cells WHO SPIKE AT ALL that are selective and these numbers are in line
with prevfious reports

3. i need to recalculate the overlap between these cells and the coherence
results, especially for cross-region coherence.  I'll have to use claires
code because now we are using claires dataset instead of mine.