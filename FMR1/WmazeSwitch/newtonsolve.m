function [q, timefail] = newtonsolve(muone,  qold, sigoldsq, mm, ll)
    							

%Solve the posterior mode equation using Newtons method
%variable  it(i)     is the estimate of posterior mode


it(1) = qold + sigoldsq*(mm - ll*exp(muone)*exp(qold)/(1 ... 
                                  + exp(muone)*exp(qold)));

for i = 1:40
   g(i)     = qold + sigoldsq*(mm - ll*exp(muone)*exp(it(i))/...
                              (1+exp(muone)*exp(it(i)))) - it(i);
   gprime(i)= -ll*sigoldsq*exp(muone)*exp(it(i))/(1+exp(muone)*exp(it(i)))^2 - 1;
   it(i+1)  = it(i) - g(i)/gprime(i);
   q        = it(i+1);
   if abs(q-it(i))<1e-14
      %fprintf(2,'cvrged (1e-14) in %d iters \n', i);
      timefail = 0; 
      return
   end
end
if(i==40) 
   %fprintf(2,'failed to converge')
   %abs(q-it(i))
   timefail = 1;
        q-it(i)
   return 
end

end
