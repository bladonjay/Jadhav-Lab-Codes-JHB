function eps = cs_findGoodEpochs(struct,fields,cell)
%fields = array of strings, with names of struct fields to check

filt = [];
%get epochs during which cell was actually present- avoids errors when cell
%was only present on certain epochs 
for f = 1:length(fields)
    filt = [filt,'~isempty($',fields{f},')'];
    if f < length(fields)
        filt = [filt, ' && '];
    end
end
%filt = ['~isempty($',field,') && ~isempty($meanrate)' ];
cells = evaluatefilter(struct,filt);

if isempty(cells)
    eps = [];
else
eps = cells(ismember(cells(:,2:3),cell,'rows'),1);
end