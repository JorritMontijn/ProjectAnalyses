%%
L_s=9; %rate during stimulus
L_b=10; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
%m = 10;%number of trials
intReps = 10;
mMax = 100;%number of trials
vecM=10:10:mMax;
mNum = numel(vecM);

%% sims
vecEmaxDprime_sim = nan(1,mNum);
vecS_maxDprime_sim = nan(1,mNum);
vecvar_Dprime_sim = nan(1,mNum);
vecE_maxD_sim = nan(1,mNum);
vecS_maxD_sim = nan(1,mNum);
vecE_zetaP_sim = nan(1,mNum);

hTic=tic;
for i_m=1:mNum
	m=vecM(i_m);
	if toc(hTic) > 5
		m
		hTic=tic;
	end
	
	%simulate
	[E_maxDprime_sim,S_maxDprime_sim,var_Dprime_sim,E_maxD_sim,S_maxD_sim,E_zetap_sim]=getD_simulation(L_s,L_b,T,tau,m,intReps);
	vecEmaxDprime_sim(i_m)=E_maxDprime_sim;
	vecS_maxDprime_sim(i_m)=S_maxDprime_sim;
	vecvar_Dprime_sim(i_m)=var_Dprime_sim;
	vecE_maxD_sim(i_m)=E_maxD_sim;
	vecS_maxD_sim(i_m)=S_maxD_sim;
	
	
	%zeta
	[dblZetaP_sim,dblZETA] = getGumbel(E_maxDprime_sim,S_maxDprime_sim.^2,E_maxD_sim);
	vecE_zetaP_sim(i_m)=dblZetaP_sim;
end

%% theory
hTic=tic;

%theory
vecEmaxDprime = nan(1,mNum);
vecS_maxDprime = nan(1,mNum);
vecvar_Dprime = nan(1,mNum);
vecE_maxD = nan(1,mNum);
vecS_maxD = nan(1,mNum);
vecE_zetaP = nan(1,mNum);
for i_m=1:mNum
	m=vecM(i_m);
	if toc(hTic) > 5
		m
		hTic=tic;
	end
	
	%theory
	[E_maxDprime,S_maxDprime,var_Dprime,E_maxD,S_maxD,E_zetap]=getD_theory(L_s,L_b,T,tau,m);
	vecEmaxDprime(i_m)=E_maxDprime;
	vecS_maxDprime(i_m)=S_maxDprime;
	vecvar_Dprime(i_m)=var_Dprime;
	vecE_maxD(i_m)=E_maxD;
	vecS_maxD(i_m)=S_maxD;
	vecE_zetaP(i_m)=E_zetap;
	
end
% plot
figure
subplot(2,3,1)
hold on
plot(vecM,vecEmaxDprime_sim,'r');
plot(vecM,vecEmaxDprime,'b');
hold off
title('E[max(d'')]')
ylim([0 max(get(gca,'ylim'))]);

subplot(2,3,2)
hold on
plot(vecM,vecS_maxDprime_sim,'r');
plot(vecM,vecS_maxDprime,'b');
hold off
title('S[max(d'')]')
ylim([0 max(get(gca,'ylim'))]);

subplot(2,3,3)
hold on
plot(vecM,vecvar_Dprime_sim,'r');
plot(vecM,vecvar_Dprime,'b');
hold off
title('Var[d'']')
ylim([0 max(get(gca,'ylim'))]);

subplot(2,3,4)
hold on
plot(vecM,vecE_maxD_sim,'r');
plot(vecM,vecE_maxD,'b');
hold off
title('E[max(d)]')
ylim([0 max(get(gca,'ylim'))]);

subplot(2,3,5)
hold on
plot(vecM,vecS_maxD_sim,'r');
plot(vecM,vecS_maxD,'b');
hold off
title('S[max(d)]')
ylim([0 max(get(gca,'ylim'))]);

subplot(2,3,6)
hold on
plot(vecM,vecE_zetaP_sim,'r');
plot(vecM,vecE_zetaP,'b');
hold off