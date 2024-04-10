%% aim
%{
make plots of IFR distros for real, ISI shuffle, shufftid and uniform
%}

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform'};
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

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
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%% move onset
	%remove first x ms
	vecStimOnTime = vecStimOnTime + dblRemOnset;
	
	%% load prepped data
	strTarget = fullpath(strDataPathT0,[sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real') '.mat']);
	if ~exist(strTarget,'file')
		fprintf('Prepped T0 file did not exist for %s; skipping...\n',strThisRec);
		continue;
	end
	
	%check fr
	sSource = load(strTarget);
	cellSpikeTimes = sSource.cellSpikeTimes;
	intNumN = numel(cellSpikeTimes);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
	for intN=1:intNumN
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
		for intTrial=1:numel(vecStimOnTime)
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	
	%% get ori vars
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	%types: Real, UniformTrial, ShuffTid, PoissGain
	
	%vecPopActEdges = logspace(-5,-1,100);
	vecPopActEdges = linspace(1e-5,5e-2,100);
	vecISI_Centers = vecPopActEdges(2:end)-diff(vecPopActEdges(1:2))/2;
	matISIsByType = nan(numel(cellTypes),numel(vecISI_Centers));
	vecMedianISI = nan(numel(cellTypes),1);
	vecMeanISI = nan(numel(cellTypes),1);
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		if numel(cellSpikeTimes) ~= intNumN
			error('neuron # mismatch!')
		end
		if isempty(vecTime),continue;end
		
		%% take only period during stimuli
		for i=1:numel(cellSpikeTimes)
			[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellSpikeTimes{i},vecStimOnTime,dblStimDur,true);
			cellSpikeTimes{i} = vecPseudoSpikeTimes;
		end
		
		%% calculate ISIs
		vecAllSpikeT = sort(cell2vec(cellSpikeTimes));
		vecTimeNoise = vecAllSpikeT;% + (1e-5)*rand(size(vecAllSpikeT));
		vecISIs = diff(vecTimeNoise);
		vecCounts=histcounts(vecISIs,vecPopActEdges);
		matISIsByType(intType,:) = vecCounts;
		vecMedianISI(intType) = median(vecAllSpikeT);
		vecMeanISI(intType) = mean(vecAllSpikeT);
	end
	% plot
	matNormISIsByType = (1+matISIsByType) ./ (sum(matISIsByType,2)+1);
	figure
	plot(vecISI_Centers*1000,1+matISIsByType)
	set(gca,'yscale','log')
	xlabel('Inter-spike interval (ms)');
	
	%calc log-integral
	cellLegend = cellTypes;
	vecLogIntegral = sum(log10(1+matISIsByType).*diff(vecPopActEdges),2);
	for i=1:numel(cellLegend)
		cellLegend{i} = [cellLegend{i} sprintf('=%.3f',vecLogIntegral(i))];
	end
	legend(cellLegend)
	strRecTit =  strrep(strrep(strrep(strrep(strThisRec,'_t0',''),'_g0',''),'R01',''),'Rec','');
	title([strRecTit '; pop rate var (log10-integral of ISIs)'],'interpreter','none');
	ylabel('Number of spikes');
	fixfig;
	
	%%
	export_fig(fullpath(strFigurePathSR,sprintf('D1_PopActVar_%s_%s.tif',strThisRec,strOnset)));
	export_fig(fullpath(strFigurePathSR,sprintf('D1_PopActVar_%s_%s.pdf',strThisRec,strOnset)));
	
	%% save data
	save(fullpath(strTargetDataPath,sprintf('D1Data_%s_%s.mat',strThisRec,strOnset)),...
		'cellTypes',...
		'vecPopActEdges',...
		'matISIsByType',...
		'strRec',...
		'dblRemOnset',...
		'vecLogIntegral');close all;
end
toc
