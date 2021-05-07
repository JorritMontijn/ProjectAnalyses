function [E_maxDprime,S_maxDprime,var_Dprime,E_maxD,S_maxD,E_zetap]=getD_simulation(L_s,L_b,T,tau,m,reps)
	
	%%
	if ~exist('m','var') || isempty(m)
		L_s=10; %rate during stimulus
		L_b=20; %rate during baseline
		T = 1; %stim duration
		tau = 3; %total trial duration (s+b)
		m = 100;%number of trials
	end
	if ~exist('reps','var') || isempty(reps)
		reps=10; %rate during stimulus
	end
	vecE_maxD = nan(1,reps);
	vecE_maxDprime = nan(1,reps);
	vecS_maxDprime = nan(1,reps);
	vecvar_Dprime = nan(1,reps);
	vecE_zetap = nan(1,reps);
	
	vecMuD = nan(1,reps);
	vecVarD = nan(1,reps);
	
	for intRep=1:reps
		
		%d
		[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);
		dblUseMaxDur = tau;
		intResampNum = 100;
		[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum);
		
		vecE_maxD(intRep) = max(abs(vecRealDiff));
		vecE_zetap(intRep) = dblZetaP;
		%d'
		L = (T*L_s+(tau-T)*L_b)/tau;
		[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m,T,tau,L,L);
		dblUseMaxDur = tau;
		intResampNum = 100;
		[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff] = calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum);
		
		vecMuD(intRep) = mean(abs(vecRealDiff))^2;
		vecVarD(intRep) = var(abs(vecRealDiff));
		
		vecE_maxDprime(intRep) = mean(max(matRandDiff,[],1));
		vecS_maxDprime(intRep) = std(max(matRandDiff,[],1));
		vecvar_Dprime(intRep) = mean(var(matRandDiff,[],1),2);
		
	end
	%vecRatio = vecMuD./vecVarD
	
	E_maxDprime = mean(vecE_maxDprime);
	S_maxDprime = mean(vecS_maxDprime);
	var_Dprime = mean(vecvar_Dprime);
	E_maxD = mean(vecE_maxD);
	S_maxD = std(vecE_maxD);
	E_zetap = mean(vecE_zetap);