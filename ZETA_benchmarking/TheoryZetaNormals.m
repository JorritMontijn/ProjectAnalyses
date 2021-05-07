%%
L_s=10; %rate during stimulus
L_b=10; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
%m = 10;%number of trials
intReps = 10;
mMax = 100;%number of trials
vecM=10:10:mMax;
m=mMax;
mNum = numel(vecM);

%% sim
hTic=tic;
for intRep=1:intReps
		if toc(hTic) > 5
			m
			hTic=tic;
		end
		%get spike times
		[vecSpikeTimes,vecEventStarts] = getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);
		dblUseMaxDur=tau;
		[vecSpikeTimes,vecEventStarts] = getGeneratedSpikingGaussian(m,T,tau,L_b,L_s)
		
		mu = 1/L_b;
		sd = 1/L_b;
		logMu = log((mu^2)/sqrt(sd^2+mu^2));
		logVar = log(((sd^2)/(m^2))+1);
		x= lognrnd(mu,sd,[1 round(L_b*m*tau*2)]);
		mean(x)
		std(x)
		
vecISIb = normrnd(1/L_b,1/(L_b^2),[1 round(L_b*m*tau*2)]);
vecISIb(vecISIb<0)=[];
std(diff(vecISIb))

		%make normal
		dblMuISI= 1/((L_s*T + L_b*(tau-T))/tau);
		dblSdISI = dblMuISI;
		
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