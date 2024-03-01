%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: is interneuron/pyramidal activation ratio different between low/high epochs?

q2: are tuned neurons more specifically activated during low/high epochs?

q3: what is different during initial peak?

%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};

boolSaveFigs = true;
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];
dblStartT = 0;

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
	
	%prep grating data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	if numel(sUseNeuron) == 0, continue;end
	
	%% is cell an interneuron (fast/narrow spiking) or pyramid (regular/broad spiking)?
	%load waveform
	if isfield(sUseNeuron,'Waveform') && ~isempty(sUseNeuron(1).Waveform)
		dblSampRateIM = str2double(sThisSource.sMetaAP.imSampRate);
	else
		[sThisRec,sUseNeuron] = loadWaveforms(sThisRec,sUseNeuron);
		dblSampRateIM = sThisRec.sample_rate;
	end
	
	%calculate waveform properties
	dblRecDur = max(cellfun(@max,{sUseNeuron.SpikeTimes})) - min(cellfun(@min,{sUseNeuron.SpikeTimes}));
	vecSpikeRate = cellfun(@numel,{sUseNeuron.SpikeTimes})/dblRecDur;
	matAreaWaveforms = cell2mat({sUseNeuron.Waveform}'); %[cell x sample]
	intNeurons=size(matAreaWaveforms,1);
	vecSpikeDur = nan(1,intNeurons);
	vecSpikePTR = nan(1,intNeurons);
	for intNeuron=1:intNeurons
		%find trough
		[dblTroughVal,intTrough]=min(matAreaWaveforms(intNeuron,:));
		[dblPeakVal,intTroughToPeak]=max(matAreaWaveforms(intNeuron,intTrough:end));
		intPeak = intTrough + intTroughToPeak - 1;
		
		dblTroughTime = intTrough/dblSampRateIM;
		dblTroughToPeakTime = intTroughToPeak/dblSampRateIM;
		dblPeakTime = intPeak/dblSampRateIM;
		vecSpikeDur(intNeuron) = dblTroughToPeakTime;
		vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
	end
	dblPTT = 0.5;
	dblSWT = 0.5/1000;
	vecNarrow = vecSpikeDur < dblSWT & vecSpikePTR > dblPTT; %Ctx BL6
	vecBroad = vecSpikeDur > dblSWT & vecSpikePTR < dblPTT; %Ctx BL6
	vecOther = ~vecNarrow & ~vecBroad;
	vecCol = vecNarrow + vecBroad*2 + vecOther*3;
	%scatter(vecSpikeDur,vecSpikePTR,[],vecCol)
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		cellSpikeTimes = {sArea1Neurons.SpikeTimes};
		[matData,indTuned,indResp,cellSpikeTimes,sOut] = NpxPrepData(cellSpikeTimes,vecStimOnTime,vecStimOffTime,vecOrientation);
		vecOri180 = mod(vecOrientation,180)*2;
		intTunedN = sum(indTuned);
		intNumN = sum(indResp);
		
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = -dblStartT;%0.3;
		dblPostTime = 0;%0.3;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		dblBinWidth = 0.05;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		%matBNT_shifted = nan(intBinNum,intNumN,intTrialNum);
		
		vecRepCounter = zeros(1,intOriNum);
		%get spikes per trial per neuron
		cellSpikeTimesStitched = cell(intNumN,intTrialNum);
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
		cellSpikeTimesPerCellPerTrial_S = cell(intNumN,intTrialNum); %single-neuron ISIs, shuffled per trial
		for intN=1:intNumN
			% build pseudo data, stitching stimulus periods
			[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			cellSpikeTimesStitched{intN} = vecPseudoSpikeTimes;
			
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
			for intTrial=1:intTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecISI = diff(vecSpikeT);
				if isempty(vecSpikeT)
					vecGenSpikesS = [];
				else
					vecISIS = diff(vecSpikeT);
					vecGenSpikesS = cumsum([vecSpikeT(1);vecISI(randperm(numel(vecISI)))]);
				end
				
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
				cellSpikeTimesPerCellPerTrial_S{intN,intTrial} = vecGenSpikesS;
			end
		end
		vecStimOnStitched = vecPseudoEventT;
		
		%% decode using first spike delay
		%constants
		[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		intStimNr = numel(vecUnique);
		dblLambda = 1;%1
		intTypeCV = 2;
		dblUseStartT = 0;
		dblUseMaxDur = dblMaxDur-dblUseStartT;
		intUseMax = inf;
		
		%simple spike-delay "time code" = identical to hybrid time code form with spike # = 1
		matSpikeDelay = cellfun(@(x) min(x(x>dblUseStartT)),cellSpikeTimesPerCellPerTrial,'UniformOutput',false);
		matSpikeDelay(cellfun(@isempty,matSpikeDelay)) = {dblMaxDur};%what to do here?
		matSpikeDelay = cell2mat(matSpikeDelay);
		[dblPerformanceDelayCV,vecDecodedIndexDelayCV,matPosteriorProbabilityDelay,dblMeanErrorDegsDelay,matConfusionDelay,matWeightsDelay] = ...
			doCrossValidatedDecodingLR(matSpikeDelay,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
	
		%simple "rate code" = identical to transformed hybrid time code with spike # = inf
		matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
		matMeanRate = matSpikeCounts./dblUseMaxDur;
		[dblPerformanceMeanRateCV,vecDecodedIndexRateCV,matPosteriorProbabilityRate,dblMeanErrorDegsRate,matConfusionRate,matWeightsRate] = ...
			doCrossValidatedDecodingLR(matMeanRate,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
	
		%simple "rate-as-time code" = identical to raw hybrid time code with spike # = inf
		matRateAsTimeCode = (1./(matMeanRate+1/dblUseMaxDur));
		[dblPerformanceMeanRateAsTimeCV,vecDecodedIndexRateCV,matPosteriorProbabilityRate,dblMeanErrorDegsRate,matConfusionRate,matWeightsRate] = ...
			doCrossValidatedDecodingLR(matRateAsTimeCode,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
	
		%calculate hybrid form for various spike # inclusions
		%variables
		vecAllSpikeCounts = matSpikeCounts(:);
		[varDataOut,vecSpikesPerTrialPerCell,vecCounts,cellSelect,vecRepetition] = val2idx(vecAllSpikeCounts);
		intTotSpikeNum = sum(vecAllSpikeCounts);
		intUseMaxIdx = find(vecSpikesPerTrialPerCell <= intUseMax,1,'last');
		vecRunThresholds = vecSpikesPerTrialPerCell(2:intUseMaxIdx);
		intMaxIdx = numel(vecRunThresholds);
		vecSpikeFractionUsedPerThreshold = zeros(size(vecRunThresholds));
		vecPerfTimeCode = zeros(size(vecRunThresholds));
		vecPerfTimeAsRateCode = zeros(size(vecRunThresholds));
		vecPerfFusionCode = zeros(size(vecRunThresholds));
		for intIdxT = 1:numel(vecRunThresholds)
			%calc spikes used
			intSpikeNum = vecRunThresholds(intIdxT);
			fprintf('Decoding using max %d spikes (%d/%d) [%s]\n',intSpikeNum,intIdxT,numel(vecRunThresholds),getTime);
			vecSpikeTrainsWithLeftovers = vecCounts.*(vecSpikesPerTrialPerCell > intSpikeNum);
			vecSpikesLeftOver = (vecSpikeTrainsWithLeftovers>0) .* (vecSpikesPerTrialPerCell - intSpikeNum);
			intTotLeftOver = sum(vecSpikeTrainsWithLeftovers.*vecSpikesLeftOver);
			vecSpikeFractionUsedPerThreshold(intIdxT) = (intTotSpikeNum-intTotLeftOver)/intTotSpikeNum;
			
			%get time code
			matTimeCode = cellfun(@(vecSpikeT) getTimeCode(vecSpikeT,dblUseStartT,intSpikeNum,dblUseMaxDur),cellSpikeTimesPerCellPerTrial);
			matTimeAsRateCode = (1./matTimeCode) - 1/dblUseMaxDur;
			%hybrid formulation; progresses from time to rate code
			[dblPerformanceTimeCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(matTimeCode,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			[dblPerformanceTimeAsRateCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(matTimeAsRateCode,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			
			%% decode using both
			[dblPerformanceFusionCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(cat(1,matTimeAsRateCode,matTimeCode),vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			
			%save
			vecPerfFusionCode(intIdxT) = dblPerformanceFusionCV;
			vecPerfTimeCode(intIdxT) = dblPerformanceTimeCV;
			vecPerfTimeAsRateCode(intIdxT) = dblPerformanceTimeAsRateCV;
		end
		%save
		vecRunThresholds = cat(1,0,vecRunThresholds);
		vecSpikeFractionUsedPerThreshold = cat(1,0,vecSpikeFractionUsedPerThreshold);
		vecPerfTimeCode = cat(1,1/intStimNr,vecPerfTimeCode);
		vecPerfTimeAsRateCode = cat(1,1/intStimNr,vecPerfTimeAsRateCode);
		vecPerfFusionCode = cat(1,1/intStimNr,vecPerfFusionCode);
		
		%% plot
		vecColTime = lines(1);
		vecColRate = [0.8 0 0];
		figure;maxfig;
		subplot(2,3,1)
		intMaxSpikes = max(vecRunThresholds);
		plot(vecRunThresholds,vecSpikeFractionUsedPerThreshold)
		ylabel('Fraction of spikes used');
		xlabel('Max. # of spikes used per cell per stim');
		fixfig;
		
		subplot(2,3,2)
		hold on
		plot([0 intMaxSpikes],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 intMaxSpikes],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot(vecRunThresholds,vecPerfTimeCode,'color',vecColTime);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('Max. # of spikes used per cell per stim');
		legend({'Chance','Rate-as-time code','Time code'},'Location','best');
		fixfig;
		
		subplot(2,3,3)
		hold on
		plot([0 intMaxSpikes],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 intMaxSpikes],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecRunThresholds,vecPerfTimeAsRateCode,'color',vecColRate);
		hold off
		ylabel('Decoding accuracy');
		title('Time-as-rate code (mean ISI)');
		xlabel('Max. # of spikes used per cell per stim');
		legend({'Chance','Rate code','Time-as-rate code'},'Location','best');
		fixfig;
		
		subplot(2,3,4)
		intStimNr = numel(vecUnique);
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold)],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot([0 max(vecSpikeFractionUsedPerThreshold)],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold,vecPerfTimeAsRateCode,'color',vecColRate);
		plot(vecSpikeFractionUsedPerThreshold,vecPerfTimeCode,'color',vecColTime);
		plot(vecSpikeFractionUsedPerThreshold,vecPerfFusionCode,'color',(vecColRate+vecColTime)/2);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Full time code','Full rate code','Rate code','Time code','Fusion code'},'Location','best');
		fixfig;
		
		
		subplot(2,3,5)
		intStimNr = numel(vecUnique);
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold)],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold,vecPerfTimeCode,'color',vecColTime);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Rate-as-time code','Time code'},'Location','best');
		fixfig;
		
		subplot(2,3,6)
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold)],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold,vecPerfTimeAsRateCode,'color',vecColRate);
		hold off
		ylabel('Decoding accuracy');
		title('Time-as-rate code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Rate code','Time-as-rate code'},'Location','best');
		fixfig;
		
		if boolSaveFigs
			%% save fig
			export_fig(fullpath(strFigurePath,sprintf('2A1_TimeCodingT%s_%s.tif',num2str(dblStartT),strRec)));
			export_fig(fullpath(strFigurePath,sprintf('2A1_TimeCodingT%s_%s.pdf',num2str(dblStartT),strRec)));
		end
		
		%% run same as above, but with time progression; use only spikes in first x ms
		%constants
		[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		intStimNr = numel(vecUnique);
		dblStep = 10/1000;%10ms
		vecEndTimes = (dblStep:dblStep:dblUseMaxDur)+dblUseStartT;
		
		%calculate hybrid form for various spike # inclusions
		%variables
		intTotSpikeNum = sum(matSpikeCounts(:));
		vecSpikeFractionUsedPerThreshold2 = zeros(size(vecEndTimes));
		vecPerfTimeCode2 = zeros(size(vecEndTimes));
		vecPerfTimeAsRateCode2 = zeros(size(vecEndTimes));
		vecPerfFusionCode2 = zeros(size(vecEndTimes));
		for intIdxT = 1:numel(vecEndTimes)
			%calc spikes used
			dblEndTime = vecEndTimes(intIdxT);
			fprintf('Decoding using %.3fs (%d/%d) [%s]\n',dblEndTime,intIdxT,numel(vecEndTimes),getTime);
			
			%select spikes
			cellCroppedSpikeT = cellfun(@(x) x(x<dblEndTime),cellSpikeTimesPerCellPerTrial,'UniformOutput',false);
			vecSpikeFractionUsedPerThreshold2(intIdxT) = sum(flat(cellfun(@numel,cellCroppedSpikeT)))/intTotSpikeNum;
			
			%get time code
			matRateCode = cellfun(@numel,cellCroppedSpikeT)./dblEndTime;
			matTimeCode = cellfun(@(vecSpikeT) getTimeCode(vecSpikeT,dblUseStartT,inf,dblEndTime),cellCroppedSpikeT);
			matTimeAsRateCode = (1./matTimeCode) - 1/dblEndTime;
			%hybrid formulation; progresses from time to rate code
			[dblPerformanceTimeCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(matTimeCode,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			[dblPerformanceTimeAsRateCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(matTimeAsRateCode,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			[dblPerformanceFusionCV,vecDecodedIndexTimeCV,matPosteriorProbabilityTime,dblMeanErrorDegsTime,matConfusionTime,matWeightsTime] = ...
				doCrossValidatedDecodingLR(cat(1,matTimeAsRateCode,matTimeCode),vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			
			%save
			vecPerfTimeCode2(intIdxT) = dblPerformanceTimeCV;
			vecPerfTimeAsRateCode2(intIdxT) = dblPerformanceTimeAsRateCV;
			vecPerfFusionCode2(intIdxT) = dblPerformanceFusionCV;
		end
		
		%% plot
		figure;maxfig;
		subplot(2,3,1)
		dblMaxT= max(vecEndTimes);
		plot(vecEndTimes,vecSpikeFractionUsedPerThreshold2)
		ylabel('Fraction of spikes used');
		xlabel('End of epoch (s)');
		fixfig;
		
		subplot(2,3,2)
		hold on
		plot([0 dblMaxT],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 dblMaxT],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot(vecEndTimes,vecPerfTimeCode2,'color',vecColTime);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('End of epoch (s)');
		legend({'Chance','Rate-as-time code','Time code'},'Location','best');
		fixfig;
		
		subplot(2,3,3)
		hold on
		plot([0 dblMaxT],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 dblMaxT],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecEndTimes,vecPerfTimeAsRateCode2,'color',vecColRate);
		hold off
		ylabel('Decoding accuracy');
		title('Time-as-rate code (mean ISI)');
		xlabel('End of epoch (s)');
		legend({'Chance','Rate code','Time-as-rate code'},'Location','best');
		fixfig;
		
		subplot(2,3,4)
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold2,vecPerfTimeAsRateCode2,'color',vecColRate);
		plot(vecSpikeFractionUsedPerThreshold2,vecPerfTimeCode2,'color',vecColTime);
		plot(vecSpikeFractionUsedPerThreshold2,vecPerfFusionCode2,'color',(vecColRate+vecColTime)/2);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Full time code','Full rate code','Rate code','Time code','Fusion code'},'Location','best');
		fixfig;
		
		subplot(2,3,5)
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],dblPerformanceMeanRateAsTimeCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold2,vecPerfTimeCode2,'color',vecColTime);
		hold off
		ylabel('Decoding accuracy');
		title('Time code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Rate-as-time code','Time code'},'Location','best');
		fixfig;
		
		subplot(2,3,6)
		hold on
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot([0 max(vecSpikeFractionUsedPerThreshold2)],dblPerformanceMeanRateCV*[1 1],'--','color',[0 0 0]);
		plot(vecSpikeFractionUsedPerThreshold2,vecPerfTimeAsRateCode2,'color',vecColRate);
		hold off
		ylabel('Decoding accuracy');
		title('Time-as-rate code (mean ISI)');
		xlabel('Fraction of spikes used');
		legend({'Chance','Rate code','Time-as-rate code'},'Location','best');
		fixfig;
		
		if boolSaveFigs
			%% save fig
			export_fig(fullpath(strFigurePath,sprintf('2A2_TimeCoding2T%s_%s.tif',num2str(dblStartT),strRec)));
			export_fig(fullpath(strFigurePath,sprintf('2A2_TimeCoding2T%s_%s.pdf',num2str(dblStartT),strRec)));
		end
		
	end
end
toc
