%% overview
%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
%gain-scaling axis?

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
runHeaderPopTimeCoding;
boolMakeFigs = true;
%vecTimescales = 0.01:0.01:10;%10;1.5;
vecTimescales = [1e-2 1e-1 1e-0]; %timescales for rapid flucts, noise corrs, tuning
vecJitter = [0 1e3];%[0 (2.^(-9:10))];
intPopSize = inf; %24 (smallest pop of all recs) or inf (uses full pop for each rec)

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecOrigStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecOrigStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	%% go through types
	clear sAggData;
	vecRunTypes = 1:numel(cellTypes);
	for intType=vecRunTypes
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		cellOrigSpikeTimes = sSource.cellSpikeTimes;
		vecPseudoStartT = vecOrigStimOnTime;
		
		%% take only period during stimuli
		intTimescaleNum = numel(vecTimescales);
		intNumN = numel(cellOrigSpikeTimes);
		cellSpikeTimes = cell(1,intNumN);
		vecBinNumPerTimescale = floor(dblTotDur./vecTimescales);
		matMean = nan(intNumN,intTimescaleNum);
		matSd = nan(intNumN,intTimescaleNum);
		for intScale=1:intTimescaleNum
			cellCounts{intScale} = nan(intNumN,vecBinNumPerTimescale(intScale));
		end
		for i=1:intNumN
			[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellOrigSpikeTimes{i},vecOrigStimOnTime,dblStimDur,true);
			
			dblStart = vecPseudoStartT(1);
			dblEnd = vecPseudoStartT(end)+dblStimDur;
			dblTotDur = dblEnd-dblStart;
			vecPseudoSpikeTimes(vecPseudoSpikeTimes<dblStart | vecPseudoSpikeTimes>dblEnd) = [];
			vecPseudoSpikeTimes = vecPseudoSpikeTimes-dblStart;
			dblLastSpike = max(vecPseudoSpikeTimes);
			
			cellSpikeTimes{i} = vecPseudoSpikeTimes;
			
			for intScale=1:intTimescaleNum
				vecBins = 0:vecTimescales(intScale):dblTotDur;
				
				%real
				vecCounts = histcounts(cellSpikeTimes{i},vecBins);
				matMean(i,intScale) = mean(vecCounts);
				matSd(i,intScale) = std(vecCounts);
				cellCounts{intScale}(i,:)=vecCounts;
			end
		end
		
		%% analysis
		for intScale=1:intTimescaleNum
			matCounts = cellCounts{intScale};
			[intNumN,intNumT] = size(matCounts);
			vecRealMean = matMean(:,intScale);
			vecRealSd = matSd(:,intScale);
			
			vecRealCV = vecRealSd./vecRealMean;
			
			%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
			%gain-scaling axis?
			
			
			
			
			%% fit gaussian model and compare real distro to expected distro (q-q plot)
			%what do we want?
			vecReqMu = vecRealMean;
			vecReqSd = vecRealSd;
			vecReqVar = vecReqSd.^2;
			
			%log normal
			vecSd = sqrt(log((vecReqSd.^2)./(vecReqMu.^2)+1));
			vecMu = log(vecReqMu) - 0.5*(vecSd.^2);
			vecVar = vecSd.^2;
			matCovar = diag(vecVar);
			
			%expectation for exp(mvnrand):
			vecExpectedMu = exp(vecMu + 0.5*vecVar);
			matExpectedVar = exp(vecMu + vecMu' + 0.5*(vecVar + vecVar'))*(exp(matCovar)-1);
			matExpectedSd = sqrt(matExpectedVar);
			vecExpectedSd = diag(matExpectedSd);
			
			%example gauss; 2*n free params: vecMu, vecVar
			matCountsGauss = exp(mvnrnd(vecMu,matCovar,intNumT));
			
			%% fit gain-scaling model and compare real
			%example gain; n+3 free params: vecMu, dblCV, dblRange
			vecGainAxis = vecReqMu;
			[vecProjectedLocation,matProjectedPoints,vecProjLocDimNorm] = getProjOnLine(matCounts,vecGainAxis);
			dblGainRange = std(vecProjectedLocation);
			dblGainMean = sqrt(sum(mean(matCounts,1).^2));
			matResiduals = matProjectedPoints-matCounts;
			vecOffAxis = sqrt(sum(matResiduals.^2));
			vecCV = vecOffAxis'./vecProjectedLocation;
			dblCV = 0;%mean(vecCV);
			%dblGainMean1 = mean(vecProjectedLocation);
			%dblGainMean2 = mean(sqrt(sum(matActGauss.^2,2)));
			vecRandPopGainPerTrial = normrnd(dblGainMean,dblGainRange,[1 intNumT]);
			
			%decompose sd into elements for each neuron
			matCountsGain = nan(size(matActGauss));
			for intT = 1:intNumT
				dblThisGain = vecRandPopGainPerTrial(intT);
				vecOnAxis =  vecGainAxis.*dblThisGain/dblGainMean;
				vecOffAxis = dblCV*dblThisGain/dblGainMean;
				
				%get lognormal
				vecSd = sqrt(log((vecOffAxis.^2)./(vecOnAxis.^2)+1));
				vecMu = log(vecOnAxis) - 0.5*(vecOffAxis.^2);
				vecVar = vecSd.^2;
				matCovar = diag(vecVar);
				
				%matCountsGain(intT,:) = exp(mvnrnd(vecMu,matCovar,1));
				matCountsGain(intT,:) = mvnrnd(vecOnAxis,0*matCovar,1);
				
				%project on gain axis and set to
				% 		vecRandNormalPerNeuron = randn(size(vecGainAxis));
				% 		vecRandOnAxis = vecGainAxis.*vecRandPopGainPerTrial(intT);
				% 		dblRandSd = abs(normrnd(0,dblGainRange));
				% 		vecRandOffAxis = vecRandOnAxis.*dblRandSd.*vecRandNormalPerNeuron;%(vecRandNormalPerNeuron./norm(vecRandNormalPerNeuron));
				% 		vecTotAct = vecRandOnAxis+vecRandOffAxis;
				% 		matActGain(intT,:) = vecTotAct;
			end
			% q-q plot
			figure
			subplot(2,3,1)
			vecRealS = sort(matCounts(:));
			vecGaussS = sort(matCountsGauss(:));
			vecGainS = sort(matCountsGain(:));
			vecRealQ = cumsum(vecRealS)./sum(vecRealS);
			vecGaussQ = cumsum(vecGaussS)./sum(vecGaussS);
			vecGainQ = cumsum(vecGainS)./sum(vecGainS);
			
			scatter(vecRealQ,vecGaussQ,'.')
			hold on
			plot([0 1],[0 1],'k--');
			ylim([0 1]);
			
			subplot(2,3,2)
			scatter(vecRealQ,vecGainQ,'.')
			hold on
			plot([0 1],[0 1],'k--');
			ylim([0 1]);
			
		end
		%%
		%check
		
		%example gain; n+3 free params: vecMu, dblCV, dblRange
		[vecProjectedLocationG,matProjectedPointsG,vecProjLocDimNormG] = getProjOnLine(matCountsGain',vecGainAxis);
		dblGainRangeG = std(vecProjectedLocationG)
		dblGainMeanG = sqrt(sum(mean(matCountsGain,1).^2))
		%dblGainMean1 = mean(vecProjectedLocation);
		%dblGainMean2 = mean(sqrt(sum(matActGauss.^2,2)));
		
		matActGauss = matActGauss;
		%matActGain = exp(matActGain);
		
		
		
		%plot
		figure
		subplot(2,3,1)
		scatter(matActGauss(:,1),matActGauss(:,2),'.')
		hold on
		scatter(vecExpectedMu(1),vecExpectedMu(2),'kx');
		hold off
		xlim([0 5]);
		ylim([0 5]);
		set(gca,'xscale','log','yscale','log');
		
		subplot(2,3,2)
		scatter(matCountsGain(:,1),matCountsGain(:,2),'.')
		hold on
		plot([0 norm(vecReqMu)*vecGainAxis(1)],[0 norm(vecReqMu)*vecGainAxis(2)],'k');
		hold off
		xlim([0 5]);
		ylim([0 5]);
		set(gca,'xscale','log','yscale','log');
		return
	end
	return
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q3Data%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
	close all;
end
toc
