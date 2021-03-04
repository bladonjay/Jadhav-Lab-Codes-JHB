mymap=[3	, nan; 2	,42; 1	,41; 4	,44; 6	,43; 5	,46; 8	,48; 7	,45; 10	,47; 9	,57; 12	,33;
11	,34; 13	,35; 14	,36; 15	,38; 16	,37; 24	,40; 23	,39; 22	,64; 21	,63; 20	,62; 19	,61; 18	,60;
17	,59; 32	,58; 31	,56; 30	,54; 29	,55; 28	,53; 25	,52; 26	,51; 27	,50; nan	,49];


top=mymap(:,1);
bottom=mymap(:,2);

newmap=[circshift(top,1) circshift(bottom,-1)];

finalmap=linearize(newmap');
finalmap(isnan(finalmap))=[];

cellmap={};
for i=1:length(finalmap)
    cellmap{i,1}=['Mux ' num2str(finalmap(i)) ' =  ',...
        num2str(100 + 50*ceil(i/2)) ', ',...
        num2str(140 +60*rem(i,2))];
end
