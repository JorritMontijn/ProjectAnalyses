%% overview
%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
%gain-scaling axis?

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real'};%{'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
runHeaderPopTimeCoding;
boolMakeFigs = true;
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
		dblTotDur=dblStimDur*numel(vecOrigStimOnTime);
		dblStart = 0;
		dblEnd = dblTotDur;
		vecBinNumPerTimescale = floor(dblTotDur./vecTimescales);
		matMean = nan(intNumN,intTimescaleNum);
		matSd = nan(intNumN,intTimescaleNum);
		for intScale=1:intTimescaleNum
			cellCounts{intScale} = nan(intNumN,vecBinNumPerTimescale(intScale));
		end
		for i=1:intNumN
			[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellOrigSpikeTimes{i},vecOrigStimOnTime,dblStimDur,true);
			
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
		close all;
		clear sData;
		for intScale=1:intTimescaleNum
			%%
			dblTimescale = vecTimescales(intScale);
			matCounts = cellCounts{intScale};
			[intNumN,intNumT] = size(matCounts);
			vecRealMean = matMean(:,intScale);
			vecRealSd = matSd(:,intScale);
			
			vecRealCV = vecRealSd./vecRealMean;
			
			%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
			%gain-scaling axis?
			
			
			%% define source parameters
			vecReqMu = vecRealMean;
			vecReqSd = vecRealSd;
			vecReqVar = vecReqSd.^2;
			matReqCov = diag(vecReqVar);
			
			%% fit gaussian model and compare real distro to expected distro (q-q plot)
			%multi-variate gauss; 2*n free params: vecMean, vecSd
			
			%example gauss; 2*n free params: vecMu, vecVar
			matCountsGauss = logmvnrnd(vecReqMu,matReqCov,intNumT)';
			
			%% fit simple gain-scaling model
			%gain with no off-axis noise; n+2 free params: vecGainAxis, dblGainMean, dblGainSd
			vecGainAxis = vecReqMu/norm(vecReqMu);
			[vecProjectedLocation,matProjectedPoints,vecProjLocDimNorm] = getProjOnLine(matCounts,vecGainAxis);
			dblGainRange = std(vecProjectedLocation);
			dblGainMean = mean(vecProjectedLocation);
			
			%generate random gains
			vecPopGainPerTrial = logmvnrnd(dblGainMean,dblGainRange^2,intNumT);
			
			%prediction is on-axis gain
			matCountsGain = nan(size(matCountsGauss));
			for intT = 1:intNumT
				dblThisGain = vecPopGainPerTrial(intT);
				vecOnAxisAct = vecGainAxis * dblThisGain;
				matCountsGain(:,intT) = vecOnAxisAct;
			end
			
			%% fit gain-scaling model
			%gain with off-axis noise; n+2 free params: vecGainAxis, dblGainMean, dblGainSd
			[vecRealS,vecReorder] = sort(sum(matCounts,1));
			matResiduals = matProjectedPoints-matCounts;
			vecOffAxis = sqrt(sum(matResiduals.^2));
			vecCV = vecOffAxis'./(vecProjectedLocation+eps);
			dblCV = median(vecCV);
			
			matScaledResid = (matResiduals./(vecProjectedLocation+eps)');
			dblResidWidth = std(matScaledResid(:));
			
			%prediction is on-axis gain plus off-axis noise
			matCountsGainFull = nan(size(matCountsGauss));
			for intT = 1:intNumT
				dblThisGain = vecPopGainPerTrial(intT);
				vecOnAxisAct = vecGainAxis * dblThisGain;
				
				%create "residuals"
				vecOffAxisScale = dblCV * vecOnAxisAct;
				vecRandOffset = normrnd(0,dblCV,size(vecOnAxisAct));
				vecOffAxisAct = vecOffAxisScale.*vecRandOffset;
				
				matCountsGainFull(:,intT) = max(0,vecOnAxisAct+vecOffAxisAct);
			end
			
			%% plot
			%get statistics
			vecRealS = sort(sum(matCounts,1));
			vecGaussS = sort(sum(matCountsGauss,1));
			vecGainS = sort(sum(matCountsGain,1));
			vecGainFullS = sort(sum(matCountsGainFull,1));
			dblMaxLim = max(cat(1,vecRealS(:),vecGaussS(:),vecGainS(:),vecGainFullS(:)));
			
			dblR2_Gauss = getR2(vecRealS,vecGaussS);
			dblR2_Gain = getR2(vecRealS,vecGainS);
			dblR2_GainFull = getR2(vecRealS,vecGainFullS);
			
			vecRealS_Single = sort(matCounts(:));
			vecGaussS_Single = sort(matCountsGauss(:));
			vecGainS_Single = sort(matCountsGain(:));
			vecGainFullS_Single = sort(matCountsGainFull(:));
			dblMaxLim_Single = max(cat(1,vecRealS_Single(:),vecGaussS_Single(:),vecGainS_Single(:),vecGainFullS_Single(:)));
			
			dblR2_Gauss_Single = getR2(vecRealS_Single,vecGaussS_Single);
			dblR2_Gain_Single = getR2(vecRealS_Single,vecGainS_Single);
			dblR2_GainFull_Single = getR2(vecRealS_Single,vecGainFullS_Single);
			
			if boolMakeFigs
				% q-q plot
				figure;maxfig;
				subplot(2,3,1)
				scatter(vecRealS,vecGaussS,'.')
				hold on
				plot([0 dblMaxLim],[0 dblMaxLim],'k--');
				xlabel('Real pop spike count')
				ylabel('Mvn model pop spike count')
				title(sprintf('Sorted distro, t=%.3fs; Gauss R^2=%.3f',dblTimescale,dblR2_Gauss));
				
				subplot(2,3,2)
				scatter(vecRealS,vecGainS,'.')
				hold on
				plot([0 dblMaxLim],[0 dblMaxLim],'k--');
				xlabel('Real pop spike count')
				ylabel('Gain-only model pop spike count')
				title(sprintf('%s; Gain R^2=%.3f',strRec,dblR2_Gain),'interpreter','none');
				
				subplot(2,3,3)
				scatter(vecRealS,vecGainFullS,'.')
				hold on
				plot([0 dblMaxLim],[0 dblMaxLim],'k--');
				xlabel('Real pop spike count')
				ylabel('Full gain model pop spike count')
				title(sprintf('Sorted distro; Full gain R^2=%.3f',dblR2_GainFull));
				
				% single neurons
				subplot(2,3,4)
				if numel(vecRealS_Single)>1000
					vecPlot=round(linspace(1,numel(vecRealS_Single),1000));
				else
					vecPlot=1:numel(vecRealS_Single);
				end
				
				scatter(vecRealS_Single(vecPlot),vecGaussS_Single(vecPlot),'.')
				hold on
				plot([0 dblMaxLim_Single],[0 dblMaxLim_Single],'k--');
				xlabel('Real single neuron spike count')
				ylabel('Mvn model single neuron spike count')
				title(sprintf('Sorted distro; Gauss R^2=%.3f',dblR2_Gauss_Single));
				
				subplot(2,3,5)
				scatter(vecRealS_Single(vecPlot),vecGainS_Single(vecPlot),'.')
				hold on
				plot([0 dblMaxLim_Single],[0 dblMaxLim_Single],'k--');
				xlabel('Real single neuron spike count')
				ylabel('Gain only model single neuron spike count')
				title(sprintf('Sorted distro; Gain R^2=%.3f',dblR2_Gain_Single));
				
				subplot(2,3,6)
				scatter(vecRealS_Single(vecPlot),vecGainFullS_Single(vecPlot),'.')
				hold on
				plot([0 dblMaxLim_Single],[0 dblMaxLim_Single],'k--');
				xlabel('Real single neuron spike count')
				ylabel('Full gain model single neuron spike count')
				title(sprintf('Sorted distro; Full gain R^2=%.3f',dblR2_GainFull_Single));
				fixfig;
				
				%% save figure
				export_fig(fullpath(strFigurePathSR,sprintf('Q3A_GaussGainPrediction_%s_%s_T%d_%s.tif',strThisRec,strType,intScale,strOnset)));
				export_fig(fullpath(strFigurePathSR,sprintf('Q3A_GaussGainPrediction_%s_%s_T%d_%s.pdf',strThisRec,strType,intScale,strOnset)));
			end
			%save data
			sData(intTimescaleNum).strRec = strThisRec;
			sData(intTimescaleNum).dblTimescale = dblTimescale;
			sData(intTimescaleNum).strType = strType;
			sData(intTimescaleNum).strOnset = strOnset;
			sData(intTimescaleNum).dblR2_Gauss = dblR2_Gauss;
			sData(intTimescaleNum).dblR2_Gain = dblR2_Gain;
			sData(intTimescaleNum).dblR2_GainFull = dblR2_GainFull;
			sData(intTimescaleNum).dblR2_Gauss_Single = dblR2_Gauss_Single;
			sData(intTimescaleNum).dblR2_Gain_Single = dblR2_Gain_Single;
			sData(intTimescaleNum).dblR2_GainFull_Single = dblR2_GainFull_Single;
		end
		
		%% agg data
		eval(['sAggData.sData' strType ' = sData;']);
	end
	
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q3Data%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
end
toc
