function  [p, qhat, sigsq, qhatold, sigsqold] ... 
		= recfilter(I, sige, qguess, sigsqguess, muone)

%implements the forward recursive filtering algorithm
%on the spike train data I
%variables:
%        qhatold      one-step prediction
%        sigsqold     one-step prediction variance
%        qhat         posterior mode
%        sigsq        posterior variance


T                  = size(I,2);
mm                 = I(1,:);
ll                 = I(2,:);

%set up some initial values

qhat=zeros(1,T+1);
qhat(1)   = qguess;    
sigsq(1)  = sigsqguess; 


%loop through all time
%calls newtonsolve to find solution to nonlinear posterior prediction estimate

count = 1;

%number_fail saves the time steps if Newton method fails
number_fail = [];

for t=2:T+1
   qhatold(t)  = qhat(t-1);
   sigsqold(t) = sigsq(t-1) + sige^2;

   [qhat(t),flagfail] = newtonsolve(muone,  qhatold(t), sigsqold(t), mm(t-1), ll(t-1));

   if flagfail>0
      number_fail = [number_fail t];
   end
   denom       = -1/sigsqold(t) - ll(t-1)*exp(muone)*exp(qhat(t))/ ...
                                   (1+exp(muone)*exp(qhat(t)))^2 ;
   sigsq(t)    = -1/denom;
end

if isempty(number_fail)<1
   error('Newton convergence failed at times %d \n', number_fail)
end

p = exp(muone)*exp(qhat)./(1+exp(muone)*exp(qhat));



