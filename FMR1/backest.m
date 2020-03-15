function [qnew, signewsq, a] = backest(q, qold, sigsq, sigsqold)

%backward filter 
%variables

T = size(q,2);

qnew(T)     = q(T);
signewsq(T) = sigsq(T);
%a=zeros(T-2,1);
for i = T-1 :-1: 2
   a(i)        = sigsq(i)/sigsqold(i+1);
   qnew(i)     = q(i) + a(i)*(qnew(i+1) - qold(i+1));
   signewsq(i) = sigsq(i) + a(i)*a(i)*(signewsq(i+1)-sigsqold(i+1));
end

t=1;


end


