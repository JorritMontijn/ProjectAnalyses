%% build erlang from mean-centered delta-dots and compare with non-mean-centered delta-dots erlang

L_s=1; %rate during stimulus
L_b=1; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
m = 100;%number of trials

n_s = m*T*L_s;
n_b = m*(tau-T)*L_b;
n=n_s + n_b;

intSavePoints = round(n);
intNumR = 1000;
intReps=100;

matMaxD = nan(intReps,intNumR);
matVarD_i = nan(intSavePoints,intNumR);
matVarOverD = nan(intReps,intNumR);
matVarDelta_i = nan(intSavePoints,intNumR);
matMeanDelta_i = nan(intSavePoints,intNumR);

matVarOverDelta= nan(intReps,intNumR);
matMeanOverDelta= nan(intReps,intNumR);
for intR=1:intNumR
	intR



matDeltadot_i = nan(intReps,intSavePoints);
matDelta_i = nan(intReps,intSavePoints);
matDeltaU_i = nan(intReps,intSavePoints);

for intRep=1:intReps
%d
[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);


% prepare interpolation points
dblUseMaxDur = tau;
vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur);

% run normal
%pre-allocate
vecStimUseOnTime = vecEventStarts(:,1);
vecThisSpikeTimes = unique(getSpikeT(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur));
vecThisSpikeFracs = linspace(0,1,numel(vecThisSpikeTimes))';
vecThisFrac = interp1(vecThisSpikeTimes,vecThisSpikeFracs,vecSpikeT);

%get linear fractions
vecThisFracLinear = (vecSpikeT./dblUseMaxDur);

%calc difference
vecDelta = vecThisFrac - vecThisFracLinear;
vecDeltaDot = diff(vecDelta);
vecDelta = [0; cumsum(vecDeltaDot)];
varDeltaDotEmp = var(vecDeltaDot);
varDeltaDotThe = 1/n^2;

%% delta_i is not erlang-distributed over repetitions, but S_i is (or should be).
n_e = numel(vecThisSpikeTimes);
S = ((1/(n_e-1)) - vecDeltaDot); %v_i normalized to [0 1]

varDeltaDotEmp2 = nanvar(S); %variance is the same, but mean of S is 1/n

vecDelta_uncentered = [0; cumsum(S)]; %v_i normalized to [0 1] => vecThisFracLinear

%matDeltadot_i(intRep,1:intSavePoints) = S(1:intSavePoints); %all iid exponentials

vecDelta = interp1(linspace(0,1,n_e),vecDelta,linspace(0,1,n));
vecDelta_uncentered = interp1(linspace(0,1,n_e),vecDelta_uncentered,linspace(0,1,n));
matDelta_i(intRep,1:intSavePoints) = vecDelta(1:intSavePoints);
matDeltaU_i(intRep,1:intSavePoints) = vecDelta_uncentered(1:intSavePoints);
end

% predict brownian bridge
t=1:n;
W = (1:n)./(n^2);
W_sub=(t./(n))*(1/n);
%x=((n-t)./sqrt(n))*W_alt

%B Klar (2014/2015), gamma difference distributions:
%?=?1/?1 - ?2/?2;
%?^2=?1/?1^2 + ?2/?2^2;

%moreover,
%if X~Gamma(a,b)
%then cX~Gamma(a,b/c)
%
vecE = nan(1,n);
vecS = nan(1,n);

vecS_diff = nan(1,n);
vecE_diff = nan(1,n);
vecM_T = nan(1,n);
vecV_T = nan(1,n);
for intI=1:n
i = intI;
a1 = i;%shape (i)
b1 = 1/n;%rate (lambda)
deltaU_i_reps = gamrnd(a1,b1,[1 intReps]);

vecE(intI) = mean(deltaU_i_reps);
vecS(intI) = var(deltaU_i_reps);

B1 = (n-i)/(n);
B2 = 1-B1;%4*(0.25-((n/2-i+1)/n)^2);

a2 = n-i;%shape (i)
b2 = 1/n;%1/i;rate (lambda)
brown_rem = gamrnd(a2,b2,[1 intReps]);

delta_i_reps = B1*deltaU_i_reps-B2*brown_rem;
%delta_i_reps = 2.5*(1-dblBrownFactor)*brown_rem;
vecE_diff(intI) = mean(delta_i_reps);
vecS_diff(intI) = var(delta_i_reps);

vecM_T(intI) = B1*a1*b1 - B2*a2*b2;
vecV_T(intI) = (B1*(a1/((1/b1)^2)) + B2*(a2/((1/b2)^2)))/2;
end

%% d
matD_i = matDelta_i - mean(matDelta_i,2);
vecVarD_i = var(matD_i,[],1);
vecVarDelta_i = var(matDelta_i,[],1);
vecMeanDelta_i = mean(matDelta_i,1);
vecVarOverD = var(matD_i,[],2);
vecVarOverDelta = var(matDelta_i,[],2);

matVarD_i(:,intR) = vecVarD_i;
matVarDelta_i(:,intR) = vecVarDelta_i;
matMeanDelta_i(:,intR) = vecMeanDelta_i;
matMaxD(:,intR) = max(abs(matD_i),[],2);

matVarOverD(:,intR) = vecVarOverD;
matVarOverDelta(:,intR) = vecVarOverDelta;
matMeanOverDelta(:,intR) = mean(matDelta_i,2);


end

%% plot
figure
subplot(2,3,1)
%imagesc(matDelta_i-nanmean(matDelta_i))
plot(matDelta_i(1:10,:)')
ylabel('delta_i');
xlabel('i');
fixfig;

subplot(2,3,2)
%imagesc(matDelta_i-nanmean(matDelta_i))
hold on
plot(vecE_diff)
plot(nanmean(matDelta_i))
plot(vecM_T)
hold off
ylabel('E(delta_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run','Theory'},'location','best');
fixfig;

subplot(2,3,3)
hold on
plot(vecS_diff)
plot(nanvar(matDelta_i))
plot(vecV_T)
hold off
ylabel('Var(delta_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run','Theory'},'location','best');
fixfig;

%% theory
i = 1:n;
a1 = i;%shape (i)
b1 = 1/n;%rate (lambda)
lambda = 1/n;

B1 = (n-i)./(n);
B2 = 1-B1;%4*(0.25-((n/2-i+1)/n)^2);

a2 = n-i;%shape (i)
b2 = 1/n;%1/i;rate (lambda)

vecM_T = B1.*a1.*b1 - B2.*a2.*b2;
vecV_T = (B1.*(a1./((1./b1)^2)) + B2.*(a2./((1./b2)^2)))/2;

var_delta_i = (n.*i.*(lambda.^2) - (i.^2).*(lambda.^2))./n;

%Vd=sum((n.*i.*(lambda.^2) - (i.^2).*(lambda.^2))./(2*n.^2))
var_d = ((n-1)*(n+1)*lambda^2)/(12*n);
var_d2 = (n^2 - 1)/(12*n^3);

(sum(var_delta_i)/(2*n))

%%
%{
matY1 = [];
matY2 = [];
vecX = 1/n:1/n:1;
matY = nan(n,n);
for i=1:(n-1)
	%gamterm = gamma(n)/(gamma(i)*gamma(n-i));
	a=i;
	b=n-i;
	m=a/(a+b);
	
	%gamterm = exp(gammaln(n)-gammaln(i)-gammaln(n-i));
	%gamterm = 1/beta(n-i,i)
	%y1 = gamterm.*(vecX.^(i-1)).*((1-vecX).^(n-i-1));
	%y1b = gamterm.*(vecX.^(n-i-1)).*((1-vecX).^(i-1));
	y1 = betapdf(vecX,a,b);
	
matY1(:,i) = y1;
matY2(:,i) = y2;
end
subplot(2,3,1)
imagesc(1:n,vecX,log(1+matY1))
axis xy

subplot(2,3,2)
imagesc(1:n,vecX,log(1+matY2))
axis xy

subplot(2,3,3)
imagesc(1:n,vecX,matY2-matY1)
axis xy
%}
%%
%figure
subplot(2,3,4)
hold on
plot(mean(matMeanDelta_i,2),'b-');
plot(vecM_T,'r--')
plot(mean(matMeanDelta_i,2)-std(matMeanDelta_i,[],2),'b--')
plot(mean(matMeanDelta_i,2)+std(matMeanDelta_i,[],2),'b--')
hold off
title('Mean[delta_i]')
ylabel('Mean(delta_i)');
xlabel('i');
legend({'Empirical','Theory'},'location','best');
fixfig;


subplot(2,3,5)
hold on
plot(mean(matVarDelta_i,2),'b-');
plot(vecV_T,'r--')
plot(mean(matVarDelta_i,2)-std(matVarDelta_i,[],2),'b--')
plot(mean(matVarDelta_i,2)+std(matVarDelta_i,[],2),'b--')
hold off
title('Var[delta_i]')
ylabel('Var(delta_i)');
xlabel('i');
legend({'Empirical','Theory'},'location','best');
fixfig;

subplot(2,3,6)
hold on
plot(mean(matVarD_i,2),'b-');
plot(get(gca,'xlim'),var_d*[1 1],'r--');
plot(mean(matVarD_i,2)-std(matVarD_i,[],2),'b--')
plot(mean(matVarD_i,2)+std(matVarD_i,[],2),'b--')
hold off
vecLimY = [0 max(get(gca,'ylim'))];
ylim(vecLimY);
title('Var[D_i]')
ylabel('Var(d_i)');
xlabel('i');
legend({'Empirical','Theory'},'location','best');
fixfig;
maxfig;

return
%%
figure
subplot(2,3,1)
vecVarOverAllDelta = matVarOverDelta(:);
[phat,pci] = gamfit(vecVarOverAllDelta);
[mu,s,muci,sci]=normfit(log(vecVarOverAllDelta));
k = phat(1)
th = phat(2)
var_d
%mean must be var_d
k_t =pi;
mean_dv = k*th%=var_d
th_t = var_d/k_t;
dblMax = var_d*10;
dblStep = dblMax/100;
%vecVarOverD = vecVarOverD;
vecBins=dblStep:dblStep:dblMax;
vecBinsPlot=vecBins(2:end)-diff(vecBins(1:2))/2;
vecV = histcounts(vecVarOverAllDelta,vecBins)/numel(vecVarOverAllDelta);
hold on
plot(vecBinsPlot,vecV);
plot(vecBinsPlot,dblStep*gampdf(vecBinsPlot,k,th));
plot(vecBinsPlot,dblStep*gampdf(vecBinsPlot,k_t,th_t));
mLog = mean(log(vecVarOverAllDelta/1.4));
sLog = std(log(vecVarOverAllDelta/1.4));
plot(vecBinsPlot,(1.4/10)*normpdf(log(vecBinsPlot),mLog,sLog));
hold off

subplot(2,3,6)
vecVarOverAllD = matVarOverD(:);
[phat,pci] = gamfit(vecVarOverAllD);
[mu,s,muci,sci]=normfit(log(vecVarOverAllD));
k = phat(1)
th = phat(2)
var_d
%mean must be var_d
k_t =pi;
mean_dv = k*th%=var_d
th_t = var_d/k_t;
dblMax = var_d*10;
dblStep = dblMax/100;
%vecVarOverD = vecVarOverD;
vecBins=dblStep:dblStep:dblMax;
vecBinsPlot=vecBins(2:end)-diff(vecBins(1:2))/2;
vecV = histcounts(vecVarOverAllD,vecBins)/numel(vecVarOverAllD);
hold on
plot(vecBinsPlot,vecV);
plot(vecBinsPlot,dblStep*gampdf(vecBinsPlot,k,th));
plot(vecBinsPlot,dblStep*gampdf(vecBinsPlot,k_t,th_t));
mLog = mean(log(vecVarOverAllD/1.4));
sLog = std(log(vecVarOverAllD/1.4));
plot(vecBinsPlot,(1.4/10)*normpdf(log(vecBinsPlot),mLog,sLog));
hold off
return
1;		k = 2.9702;th =1.3967e-04;var_d =4.1666e-04;mean_dv =4.1486e-04
5;		k = 3.1648;th =2.6763e-05;var_d =8.3333e-05;mean_dv =8.4698e-05
10;		k = 3.2094;th =1.2637e-05;var_d =4.1667e-05;mean_dv =4.0557e-05
50;		k = 2.9017;th =2.8811e-06;var_d =8.3333e-06;mean_dv =8.3600e-06
100;	k = 3.1026;th =1.3928e-06;var_d =4.1667e-06;mean_dv =4.2788e-06

dblMV = var(vecVarOverD)/mean(vecVarOverD).^2

1
dblMV =0.3969
10
dblMV =0.3639
100
dblMV =0.4138
    
%%
EmaxX = mean(matMaxD(:));
VarMaxX = var(matMaxD(:));

dblAlpha=pi/8;
%fMaxE = @(N) norminv((1-dblAlpha)./((N*2)-2*dblAlpha+1));
%dblE = -fMaxE(intN);
fMaxE = @(N) norminv((N-dblAlpha)./(N-2*dblAlpha+1));

%calculate approximate variance for sample size N; Baglivo (2005)
fMaxV = @(N) ((N) ./ (((N+1).^2).*(N+2))).*(1./((normpdf(norminv(N./(N+1)))).^2));

	
dblE = fMaxE(sqrt(n));
sd_d = sqrt(var_d);

E_maxX_t = dblE*sd_d

	
dblV = fMaxV(nthroot(n,4));
sd_d = sqrt(var_d);

Var_maxX_t = dblV*var_d
