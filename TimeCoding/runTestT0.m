
%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real','Poiss','ShuffTid','ShuffTxClass','Uniform'};%,'Poiss','ShuffTid','Shuff','PoissGain','Uniform','ShuffTxClass'};
runHeaderPopTimeCoding;

%% plot spikes over time per trial for all stimulus classes
for intRec=1:intRecNum %19 || weird: 11
	%prep
	runRecPrepNpx;
	strThisRec = strRec;
	strDataPathT0 = strTargetDataPath;
	
	%get data matrix
	[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
		NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	
	%% prep figs
	figure
	h6 = subplot(2,3,6);hold on;xlabel('Real pop act per trial');ylabel('Pop act per trial');
	cellLegendScat = {};
	%vecH = [];
	for intStim=1:intOriNum
		%vecH(intStim) = subplot(3,4,intStim);hold on;
	end
	matCol = lines(numel(cellTypes));
	dblBinDur = 0.025;
	vecBinEdges = 0:dblBinDur:dblStimDur;
	vecBinC = vecBinEdges(2:end)-diff(vecBinEdges(1:2))/2;
	intBinNum = numel(vecBinC);
	matAggAct = nan(intRepNum,intBinNum,intStim,numel(cellTypes));
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType}
		
		%get cell props
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		strRealType = sSource.strType
		
		% get matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(sSource.cellSpikeTimes,vecStimOnTime,vecStimOffTime,vecOri180);
		if strcmp(strType,'Real')
			vecRealAct = mean(matData);
		else
			scatter(h6,vecRealAct,mean(matData),'.');
			cellLegendScat(end+1) = {strType};
		end
		subplot(2,3,intType)
		imagesc(matData);
		colorbar
		title(strType);
		continue;
		%% plot
		for intStim=1:intOriNum
			intStim
			vecStimTrials = find(vecOriIdx==intStim);
			for intRep=1:intRepNum
				vecStimSpikes = cell2vec(cellSpikeTimesPerCellPerTrial(:,vecStimTrials(intRep)));
				vecSpikeRates = histcounts(vecStimSpikes,vecBinEdges)./dblBinDur;
				matAggAct(intRep,:,intStim,intType) = vecSpikeRates;
			end
			vecMean = mean(matAggAct(:,:,intStim,intType),1);
			vecSd = std(matAggAct(:,:,intStim,intType),[],1);
			
			errorbar(vecH(intStim),vecBinC,vecMean,vecSd./sqrt(intRepNum),'color',matCol(intType,:));
		end
	end
	break
	for i=1:numel(vecH)
		set(vecH(i),'ylim',[0 1000]);
		xlabel(vecH(i),'Time after onset (s)');
		ylabel(vecH(i),'Pop rate (Hz)');
		title(vecH(i),sprintf('Stim %d',round(vecUnique(i))));
	end
	legend(vecH(1),cellTypes,'location','best');fixfig;
	return
end
legend(h6,cellLegendScat);