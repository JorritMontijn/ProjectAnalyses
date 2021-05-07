%%
clear all
L_s=5; %rate during stimulus
L_b=L_s; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
m = 100;%number of trials
intReps = 10;
mMax = 300;%number of trials
vecM=5:5:mMax;
mNum = numel(vecM);

matSum_s = nan(intReps,mNum);
matSum_b = nan(intReps,mNum);
vecE_ns = nan(1,mNum);
vecE_nb = nan(1,mNum);
vecSem_ns = nan(1,mNum);
vecSem_nb = nan(1,mNum);

%%{
hTic=tic;
for i_m=1:mNum
	m=vecM(i_m);
	for intRep=1:intReps
		if toc(hTic) > 5
			m
			hTic=tic;
		end
		[vecSpikeTimes,vecEventStarts] = getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);
		
		vecBins = sort(cat(1,vecEventStarts(1)-tau,vecEventStarts,vecEventStarts+T,vecEventStarts(end)+tau,vecEventStarts(end)+2*tau));
		
		[a,b]=histcounts(vecSpikeTimes,vecBins);
		vecR_s = a(2:2:(end-1));
		vecR_b = a(3:2:end);
		
		sum_s = sum(vecR_s);
		sum_b = sum(vecR_b);
		matSum_s(intRep,i_m) = sum_s;
		matSum_b(intRep,i_m) = sum_b;
	end
	Ens = m*T*L_s;
	Enb = m*(tau-T)*L_b;
	
	sem_s = sqrt(m*T*L_s);
	sem_b = sqrt(m*(tau-T)*L_b);
	
	vecE_ns(i_m) = Ens;
	vecE_nb(i_m) = Enb;
	vecSem_ns(i_m) = sem_s;
	vecSem_nb(i_m) = sem_b;
	
end
%
figure
subplot(2,3,1)
hold on
vecM_s = mean(matSum_s,1)./vecM;
vecE_s = std(matSum_s,[],1)./vecM;

%means
plot(vecM,vecE_ns./vecM,'b');
plot(vecM,vecM_s,'r');

%sems
plot(vecM,vecE_ns./vecM-vecSem_ns./vecM,'b--');
plot(vecM,vecE_ns./vecM+vecSem_ns./vecM,'b--');
plot(vecM,vecM_s-vecE_s,'r--');
plot(vecM,vecM_s+vecE_s,'r--');
hold off
xlabel('# of trials');
ylabel('Mean spiking rate');
title('Mean rate +/- sd, stimulus period')
legend({'Theoretical','Empirical'},'Location','best')
fixfig;

%}
%% delta-dot
%intReps = 100;
%mat_EmaxD_simulated = nan(intReps,mNum);
mat_VarDeltaDotEmp = nan(intReps,mNum);
vec_VarDeltaDotThe = nan(1,mNum);

intMaxI = max(vecM)*tau*L_b;
dblUseStepI = 100;
intMaxIdxI = floor(intMaxI/dblUseStepI);
mat_DeltaI_t = nan(intReps,intMaxIdxI);
mat_DeltaI_e = nan(intReps,intMaxIdxI);
mat_DeltaI_e_erlang = nan(intReps,intMaxIdxI);
%generate data
hTic=tic;
for i_m=1:mNum
	m=vecM(i_m);
	for intRep=1:intReps
		if toc(hTic) > 5
			m
			hTic=tic;
		end
		%get spike times
		[vecSpikeTimes,vecEventStarts] = getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);
		dblUseMaxDur=tau;
		% prepare interpolation points
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
		vecThisDiff = vecThisFrac - vecThisFracLinear;
		
		%%%
		S = diff(vecThisDiff);
		%S = diff(vecThisSpikeTimes)/tau;
		n = numel(vecThisSpikeTimes);
		vecDeltaDot = ((1/(n-1)) - S);
		varDeltaDotEmp = var(vecDeltaDot);
		
		
		mat_VarDeltaDotEmp(intRep,i_m) = varDeltaDotEmp;
		
		%% erlang
		intUseMaxI = min([n intMaxI]);
		vecUseI = dblUseStepI:dblUseStepI:intUseMaxI;
		lambda = n;
		vecDeltaI = vecThisDiff(vecUseI);
		
		vecDeltaI_erlang=nan(1,numel(vecUseI));
		for intIdxI = 1:numel(vecUseI)
			%erlang sample
			k = round(vecUseI(intIdxI)/sqrt(2));
			vecRand = rand(1,k);
			S = (-1/lambda)*sum(log(vecRand));
			vecDeltaI_erlang(intIdxI) = S;
			
		end
		
		vecDeltaI_t = (dblUseStepI:dblUseStepI:intMaxI)/(lambda^2);
		mat_DeltaI_t(intRep,:) = vecDeltaI_t;
		mat_DeltaI_e_erlang(intRep,1:numel(vecUseI)) = vecDeltaI_erlang;
		mat_DeltaI_e(intRep,1:numel(vecUseI)) = vecDeltaI;

	end
	varDeltaDotThe =1/(m*tau*L_b).^2;
	vec_VarDeltaDotThe(i_m) = varDeltaDotThe;
end

%%
h=subplot(2,3,2);
%means
vecMeanVDDE = mean(mat_VarDeltaDotEmp,1);
vecSdVDDE = std(mat_VarDeltaDotEmp,[],1);

hold on
plot(vecM,vecMeanVDDE,'r');
plot(vecM,vec_VarDeltaDotThe,'b');

%sems
plot(vecM,vecMeanVDDE-vecSdVDDE,'r--');
plot(vecM,vecMeanVDDE+vecSdVDDE,'r--');
xlabel('# of trials');
ylabel('E[Var[Delta_dot]]');
fixfig;
hold off
title('E[Var[Delta_dot]] +/- sd','interpreter','none')
legend({'Empirical','Theoretical'},'Location','best')
set(gca,'yscale','log')

subplot(2,3,3)
vecI = dblUseStepI:dblUseStepI:intMaxI;
cla
vecDeltaI_e = var(mat_DeltaI_e,[],1);
vecDeltaI_e_erlang = var(mat_DeltaI_e_erlang,[],1);
vecDeltaI_t = mean(mat_DeltaI_t,1);
%vecDeltaI_t = mean(mat_DeltaI_t,1);

hold on
plot(vecI,vecDeltaI_e,'r');
plot(vecI,vecDeltaI_e_erlang,'r--');
plot(vecI,vecDeltaI_t,'b');
xlabel('i');
ylabel('Var[Delta_i]');
fixfig;
hold off
title('Var[Delta_i]','interpreter','none')
legend({'Empirical','Erlang','Theoretical'},'Location','best')
%set(gca,'yscale','log')

%% var[d]
%intReps = 100;
%mat_EmaxD_simulated = nan(intReps,mNum);
vec_EmaxD = nan(1,mNum);
mat_var_D_sim = nan(intReps,mNum);
mat_MeanCovar = nan(intReps,mNum);
vec_var_D_t = nan(1,mNum);

%generate data
hTic=tic;
%figure
%hold on
vecBins = 0:0.0001:0.01;
matMeanR = nan(intReps,mNum);
for i_m=1:mNum
	m=vecM(i_m);
	vecCounts = zeros([1 numel(vecBins)-1]);
	for intRep=1:intReps
		if toc(hTic) > 5
			m
			hTic=tic;
		end
		[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);
		
		dblUseMaxDur = tau;
		intResampNum = 2;
		[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum);
		%% xcorr
		%vecR = nan(size(vecRealDiff));
		%for i=1:numel(vecRealDiff)
		%	vecR(i)=corr(vecRealDiff,circshift(vecRealDiff,i));
		%end
		%matMeanR(intRep,i_m) = mean(abs(vecR));

		%bin
		%vecCounts = vecCounts + histcounts(abs(matRandDiff),vecBins);
		
		%var[d]
		mat_var_D_sim(intRep,i_m) = var(vecRealDiff);
		%mat_var_D_sim(intRep,i_m) = var(vecRealFrac);
		
		%E[max[d]]
		EmaxDprime_simulated = mean(max(matRandDiff,[],1));
		smaxDprime_simulated = std(max(matRandDiff,[],1));
		
		EmaxD_simulated = max(abs(vecRealDiff));
		%{
vecDelta = vecRealFrac - vecRealFracLinear;
max(vecDelta)

subplot(2,3,4)
plot(vecRealFrac)
hold on
plot(vecRealFracLinear)
hold off

subplot(2,3,5)
histx(vecDelta)
		%}
		
		mat_EmaxD_simulated(intRep,i_m) = EmaxD_simulated;
		
		
	end
	
	%% var[d]
	dblConst = 1/6;
	
	%theoretical, direct, Hyper-Erlang
	N = m*L_b*tau;
	%lambda = 2*N;
	lambda = N;
	
	%vecNum = (1:N)+1;
	%vecDenom = N.*(1:N).*lambda;
	%varD = -1/(lambda^2) + sum(vecNum./vecDenom);
	%vec_var_D_t(i_m) = dblConst*varD;
	
	%theoretical, mean of delta_i, corrected
	%vecDeltaI_t = (1:N)./(N^2);
	%var_delta_dot_i_1 = mean(vecDeltaI_t);
	%var_delta_dot_i_1 = (1 + N)/(2*N^2);
	
	dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
	varD_approx = -1/(lambda^2) + (log(N) + dblEulerMascheroni + N)/(N*lambda);
	vec_var_D_t(i_m) = dblConst*varD_approx;
	
	%sum of centered exponentials
	vecN = 1:N;
	varD_approx2 = (2*log(N) + 2*dblEulerMascheroni + N - 3)/(2*N*lambda);
	vec_var_D_t2(i_m) = dblConst*varD_approx2;
	
	%empirical
	varD_real = mean(mat_var_D_sim(:,i_m));
	sdD_real = sqrt(mean(mat_var_D_sim(:,i_m)));
	vec_var_D_e(i_m) = varD_real;
	
	%% d'(x|n) ~ Prod(i=1...n)[ Binom(x|p=i/n,k=i,n=n)]
	
	%% E[max[d]]
	%contants, N
	N = m*L_b*tau;
	alpha = pi/8;
	
	%constants, sqrt(N)
	N_frac = 0.3*exp(-N/500);
	N_alt = N_frac*N + (1-N_frac)*sqrt(N);
	N_alt = N;%sqrt(N);
	
	%elfving
	%E_maxX = norminv((1-alpha)/(2*N-2*alpha+1)); %expected maximum of a single d' sample with X~N(0,1)
	E_maxX = norminv((N-alpha)/(N-2*alpha+1)); %expected maximum of a single d' sample with X~N(0,1)
	E_maxX_alt = norminv((N_alt-alpha)/(N_alt-2*alpha+1)); %expected maximum of a single d' sample with X~N(0,1)
	
	%theoretical, direct, Elfving Gumbel, sqrt(N)
	vec_EmaxD_t(i_m) = sqrt(vec_var_D_t(i_m))*E_maxX_alt;
	
	%theoretical, mean of delta_i, sqrt(N)
	vec_EmaxD_t2(i_m) = sqrt(vec_var_D_t2(i_m))*E_maxX_alt;
	
	%hybrid; empirical sd + Elfving Gumbel
	vec_EmaxD_h(i_m) = mean(sqrt(mat_var_D_sim(:,i_m)))*E_maxX_alt;
	
	%empirical
	vec_EmaxD_e(i_m) = mean(mat_EmaxD_simulated(:,i_m)) ;
	
	%% var[max(d)]
	%baglivo
	N_alt2 = N;%sqrt(sqrt(N));%nthroot(N,2);%N_alt;
	V_maxX = ((N) ./ (((N+1).^2).*(N+2))).*(1./((normpdf(norminv(N./(N+1)))).^2));
	V_maxX_alt = ((N_alt2) ./ (((N_alt2+1).^2).*(N_alt2+2))).*(1./((normpdf(norminv(N_alt2./(N_alt2+1)))).^2));
	
	%theoretical, direct, Elfving Gumbel, sqrt(N)
	vec_VmaxD_t(i_m) = sqrt(vec_var_D_t(i_m)*V_maxX_alt);
	
	%theoretical, mean of delta_i, sqrt(N)
	vec_VmaxD_t2(i_m) = sqrt(vec_var_D_t2(i_m)*V_maxX_alt);
	
	%hybrid; empirical sd + baglivo Gumbel
	vec_VmaxD_h(i_m) = mean(sqrt(mat_var_D_sim(:,i_m)*V_maxX_alt));
	
	%empirical
	vec_VmaxD_e(i_m) = std(mat_EmaxD_simulated(:,i_m)) ;
	
	%% get all
	[E_maxDprime,S_maxDprime,var_Dprime]=getD_theory(L_s,L_b,T,tau,m);
	
	%vec_var_D_t2(i_m) = var_Dprime;
	%vec_EmaxD_t2(i_m) = E_maxDprime;
	%vec_VmaxD_t2(i_m) = S_maxDprime;
	
end
%ylabel('r(d(i),d(i+x))');
%xlabel('x/N');
%title(sprintf('E[Mean r^2]=%.3f',mean(mean(mat_MeanCovar))))
%%% plot xcorr
%plot(mean(matMeanR).^2)
%1/mean(matMeanR(:).^2)

%%
%figure
subplot(2,3,4)
cla;
hold on
vecE_M = mean(sqrt(mat_var_D_sim),1);
vecE_se = std(sqrt(mat_var_D_sim),[],1);

vecT_st = sqrt(vec_var_D_t);
vecT2 = sqrt(vec_var_D_t2);


%means
plot(vecM,vecE_M,'r');
plot(vecM,vecT_st,'b');
plot(vecM,vecT2,'b--');

%sems
plot(vecM,vecE_M-vecE_se,'r--');
plot(vecM,vecE_M+vecE_se,'r--');
hold off
ylabel('s(D)');
xlabel('# of trials');
title('Standard deviation of d-vector)')
legend({'Empirical','Theory, Hyp-Erl 2N','Theory, sum(v[d-dot])'},'Location','best')
fixfig;

%
subplot(2,3,5)
cla;
hold on
vecE_M = mean(mat_EmaxD_simulated,1);
vecE_Sd = std(mat_EmaxD_simulated,[],1);
vecM_H = vec_EmaxD_h;
	
vecM_t = vec_EmaxD_t;
vecM_t2 = vec_EmaxD_t2;
%means
plot(vecM,vecE_M,'r');
plot(vecM,vecM_H,'m');
plot(vecM,vecM_t,'b');
plot(vecM,vecM_t2,'b--');


%sems
plot(vecM,vecE_M-vecE_Sd,'r--');
plot(vecM,vecE_M+vecE_Sd,'r--');
hold off
ylabel('E[Max(D)');
xlabel('# of trials');
title('E[Max(D)]')
legend({'Empirical','Hybrid','Theory, Hyp-Erl','Theory, mu(v[d-dot])'},'Location','best')
fixfig;

%
subplot(2,3,6)
cla;
hold on
vecE_var = vec_VmaxD_e;
vecM_H = vec_VmaxD_h;
	
vecM_t = vec_VmaxD_t;
vecM_t2 = vec_VmaxD_t2;
%means
plot(vecM,vecE_var,'r');
plot(vecM,vecM_H,'m');
plot(vecM,vecM_t,'b');
plot(vecM,vecM_t2,'b--');


%sems
%plot(vecM,vecE_M-vecE_Sd,'r--');
%plot(vecM,vecE_M+vecE_Sd,'r--');
hold off
ylabel('Var[Max(D)');
xlabel('# of trials');
title('Var[Max(D)]')
legend({'Empirical','Hybrid','Theory, Hyp-Erl','Theory, mu(v[d-dot])'},'Location','best')
fixfig;

