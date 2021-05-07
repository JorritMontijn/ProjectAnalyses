function [E_maxDprime,S_maxDprime,var_Dprime,E_maxD,S_maxD,dblZetaP]=getD_theory(L_s,L_b,T,tau,m)

%%
if ~exist('m','var') || isempty(m)
L_s=10; %rate during stimulus
L_b=12; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
m = 100;%number of trials
end

%% theoretical d'
%dblConst = 1/12;
N = m*T*L_s+m*(tau-T)*L_b; %sample size (number of spikes); effective sample size may be smaller than real number due to dependence

lambda = 1/N;%2*N;
%dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
%var_D = -1/(lambda^2) + (log(N) + dblEulerMascheroni + N)/(N*lambda);
%var_D = ((N-1)*(N+1)*lambda^2)/(12*N);
var_D = 1/(12*N);
var_Dprime= var_D;%dblConst*var_D;

%% E[max(d')]
alpha = pi/8;
N_e = sqrt(N);

E_maxX = norminv((N_e-alpha)./(N_e-2*alpha+1)); %expected maximum of a single d' sample with X~N(0,1)
E_maxDprime = (sqrt(var_Dprime).*E_maxX);

%% s[max(d')]
N = m*T*L_s+m*(tau-T)*L_b; %sample size (number of spikes); effective sample size may be smaller than real number due to dependence
N_v = N;%nthroot(N,4);

V_maxX = ((N_v) ./ (((N_v+1).^2).*(N_v+2))).*(1./((normpdf(norminv(N_v./(N_v+1)))).^2));
S_maxDprime = sqrt(var_Dprime.*V_maxX);
	

if nargout > 3
%% E[max(d)]

N_b = m*(tau - T)*L_b;
N_s = m*T*L_s;
N = N_b + N_s;

%signal in large N limit
E_EmaxD = abs(N_s/(N_s+N_b )-(L_b*N_s)/(N_b*L_s  + L_b*N_s ))/2;

%use E_maxDprime as gumbel mean & S_maxDprime as gumbel var to compute the
%probability that it will exceed E_EmaxD; then scale by that factor
dblBeta = (sqrt(6).*S_maxDprime)./(pi);

%define Euler-Mascheroni constant
dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)

%derive mode from mean, beta and E-M constant
dblMode = E_maxDprime - dblBeta.*dblEulerMascheroni;

%define Gumbel cdf
fGumbelCDF = @(x) exp(-exp(-((x(:)-dblMode)./dblBeta)));

dblProb = 1-fGumbelCDF(E_EmaxD);

%E_maxD = dblProb*E_maxDprime + E_EmaxD;
%S_maxD = S_maxDprime;

%% new E[max(d)]
%E_maxD = dblProb*E_maxDprime + (1-dblProb)*S_maxDprime + E_EmaxD;
E_maxD = E_EmaxD;
S_maxD = S_maxDprime;

%% zeta-p
%calculate statistical significance using Gumbel distribution
[dblZetaP,dblZETA] = getGumbel(E_maxDprime,S_maxDprime.^2,E_maxD);

end
