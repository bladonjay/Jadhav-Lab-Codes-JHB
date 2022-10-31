CA1   = [3 4 17 18 19 20 27 28 29 30 31 32];
 hpcRef=  1;
 PFC   = [7 8 9 12 13 14 15 16 21 22 23 25 26];
 pfcRef=   10;
 OB = 24;
 obRef = 24; 
 %sjcs_baselinespecgram('CS31', 1:3, [2 4 6], [CA1, PFC, OB], 0,'fpass',[0 100]);
%sjcs_baselinespecgram('CS31', 1:3, [2 4 6], [OB], 0,'fpass',[0 100]);
sjcs_baselinespecgram('CS31', 1:3, [2 4 6], [CA1, PFC, OB], 1,'fpass',[0 20]);




CA1   = [1 3 4 17 18 20 27 28 29 30 31];
 hpcRef=  32;
 PFC   = [7 8 9 10 11 13 14 16 22 23 25 26];
 pfcRef=   12;
 OB = 24;
 obRef = 24;
%sjcs_baselinespecgram('CS33', 1:4, 2, [CA1, PFC, OB], 0,'fpass',[0 100]);
sjcs_baselinespecgram('CS33', 1:4, 2, [CA1, PFC, OB], 1,'fpass',[0 40]);


 CA1   = [1 3 4 17 18 19 20 28 29 30 31 32];
 hpcRef=  2;
 PFC   = [7 8 9 10 13 14 15 16 21 22 23 24 25 26];
 pfcRef=   11;
 OB = [12];
 obRef = 12;
%sjcs_baselinespecgram('CS34', 1:4, [2, 4], [CA1, PFC, OB], 0,'fpass',[0 100]);
sjcs_baselinespecgram('CS34', 1:4, [2, 4], [CA1, PFC, OB], 1,'fpass',[0 20]);


CA1   = [1 2 3 4 17 18 20 27 28 29 30 31 32];
 hpcRef=  19;
 PFC   = [7 8 10 11 13 14 15 16 21 23 24 25 26];
 pfcRef=   22;
 OB = [12];
 obRef = 12;
%sjcs_baselinespecgram('CS35', 1:3, [2, 4], [CA1, PFC, OB], 0,'fpass',[0 100]);
sjcs_baselinespecgram('CS35', 1:3, [2, 4], [CA1, PFC, OB], 1,'fpass',[0 20]);