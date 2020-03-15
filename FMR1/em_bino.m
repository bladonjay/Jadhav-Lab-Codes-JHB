function [newsigsq] = em_bino(I, qnew, signewsq, a, muold, startflag);

%computes sigma_eps squared 

N           = size(qnew,2);

qnewt      = qnew(3:N);
qnewtm1    = qnew(2:N-1);
signewsqt  = signewsq(3:N);
a          = a(2:end);

covcalc    = signewsqt.*a;


term1      = sum(qnewt.^2) + sum(signewsqt);
term2      = sum(qnewt.^2) + sum(signewsqt);
term3      = -2*( sum(covcalc) + sum(qnewt.*qnewtm1) );

if startflag == 1
 term6      = 1.5*qnew(2)*qnew(2) + 2.0*signewsq(2) - qnew(end)^2 - signewsq(end);
elseif( startflag == 0)
 term6      = 2*qnew(2)*qnew(2) + 2*signewsq(2) - qnew(end)^2 - signewsq(end);
elseif( startflag == 2)
 term6      = 1*qnew(2)*qnew(2) + 2*signewsq(2) - qnew(end)^2 - signewsq(end);
 N = N-1;
end


newsigsq   = (term1+term2+term3+term6)/N;
end


