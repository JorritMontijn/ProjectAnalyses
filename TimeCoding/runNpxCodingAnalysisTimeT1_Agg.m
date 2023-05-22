%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
boolHome = false;
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
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
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx('','natural',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep nm data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	[vecFrameIdx,vecUniqueFrames,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimIdx);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecStimIdx));
	intRepNum = intTrialNum/intStimNum;
	
	%original movie time
	vecOrigStimOnTime = sThisRec.cellBlock{1}.vecStimOnTime;
	vecOrigStimOffTime = sThisRec.cellBlock{1}.vecStimOffTime;
	dblStimDur = min(diff(vecOrigStimOnTime));
	intOrigTrialNum = numel(vecOrigStimOnTime);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%% load prepro T0 data
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		sLoad = load(fullpath(strTargetDataPath,sprintf('T0Data_%s',strRec)));
		
		dblCutOff = 0.5;%sLoad.dblCutOff;
		vecTime_Real = sLoad.vecTime_Real;
		vecIFR_Real = sLoad.vecIFR_Real;
		vecTime_Shuff = sLoad.vecTime_Shuff;
		vecIFR_Shuff = sLoad.vecIFR_Shuff;
		vecPeakHeight_Real = sLoad.vecPeakHeight_Real;
		vecPeakHeight_Shuff = sLoad.vecPeakHeight_Shuff;
		vecPeakLocs_Real = sLoad.vecPeakLocs_Real;
		vecPeakLocs_Shuff = sLoad.vecPeakLocs_Shuff;
		vecNormTime_Real = sLoad.vecNormTime_Real;
		vecNormTime_Shuff = sLoad.vecNormTime_Shuff;
		vecAllSpikeNeuron_Real = sLoad.vecAllSpikeNeuron_Real;
		vecAllSpikeTime_Real = sLoad.vecAllSpikeTime_Real;
		
		%% rerun peaks
		% filter
		%real
		intLag = round((numel(vecIFR_Real)/numel(vecOrigStimOnTime))/2);
		dblThreshZ = 1;
		dblInfluence = 0.5;
		[signals,avgFilter_Real,stdFilter_Real] = detectpeaks(vecIFR_Real,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Real = (vecIFR_Real - avgFilter_Real)./stdFilter_Real;
		vecNormIFR_Real = (vecIFR_Real - avgFilter_Real)./avgFilter_Real;
		indRem_Real = isinf(vecNormIFR_Real);
		vecNormIFR_Real(indRem_Real) = [];
		vecNormTime_Real  = vecTime_Real(~indRem_Real);
		
		%shuff
		[signals,avgFilter_Shuff,stdFilter_Shuff] = detectpeaks(vecIFR_Shuff,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./stdFilter_Shuff;
		vecNormIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./avgFilter_Shuff;
		indRem_Shuff = isinf(vecNormIFR_Shuff);
		vecNormIFR_Shuff(indRem_Shuff) = [];
		vecNormTime_Shuff = vecTime_Shuff(~indRem_Shuff);
		
		% shuff peaks
		indStimSpikes_Shuff = vecNormTime_Shuff>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime_Shuff<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Shuff = vecNormIFR_Shuff(indStimSpikes_Shuff);
		vecStimTime_Shuff = vecNormTime_Shuff(indStimSpikes_Shuff);
			
		[vecPeakHeight_Shuff,vecPeakLocs_Shuff,w_Shuff,p_Shuff] = findpeaks(vecStimIFR_Shuff);
			
		% real peaks
		indStimSpikes_Real = vecNormTime_Real>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime_Real<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Real = vecNormIFR_Real(indStimSpikes_Real);
		vecStimTime_Real = vecNormTime_Real(indStimSpikes_Real);
			
		[vecPeakHeight_Real,vecPeakLocs_Real,w_Real,p_Real] = findpeaks(vecStimIFR_Real);
		
		%% plot
		%events
		dblCutOff = 0.8;
		dblStartEpoch = vecOrigStimOnTime(1)-dblStimDur;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur;
		
		%retain only stim epoch
		indStimEpoch = vecAllSpikeTime_Real > dblStartEpoch & vecAllSpikeTime_Real < (dblStartEpoch + dblEpochDur);
		vecPopSpikeTime_Real = vecAllSpikeTime_Real(indStimEpoch);
		vecPopNeuronId_Real = vecAllSpikeNeuron_Real(indStimEpoch);
		
		%get event times
		vecPopEventTimes_Real = vecStimTime_Real(vecPeakLocs_Real(vecPeakHeight_Real>dblCutOff));
		intPopEventNum = numel(vecPopEventTimes_Real);
		vecPopEventLocs_Real = nan(1,intPopEventNum);
		vecPopEventLocsIFR_Real = nan(1,intPopEventNum);
		
		for intEvent=1:intPopEventNum
			[dummy,vecPopEventLocs_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecPopSpikeTime_Real));
			[dummy,vecPopEventLocsIFR_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecStimTime_Real));
		end
		
		% plot
		figure
		subplot(2,3,1)
		plot(vecStimTime_Real,vecStimIFR_Real);
		hold on
		scatter(vecPopEventTimes_Real,vecStimIFR_Real(vecPopEventLocsIFR_Real));
		scatter(vecOrigStimOnTime,ones(size(vecOrigStimOnTime)),'kx');
		
		vecTrialBins = 0:0.25:dblStimDur;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes,vecTrialBins,vecOrigStimOnTime,-1);
		
		subplot(2,3,2)
		imagesc(mean(matPET,3));
		
		subplot(2,3,4)
		imagesc(matPET(:,:,1));
		subplot(2,3,5)
		imagesc(matPET(:,:,2));
		subplot(2,3,6)
		imagesc(matPET(:,:,3));
		
		%% merge peaks
		%how?...
		indData = vecStimTime_Real>2050.5 & vecStimTime_Real < 2051;
		vecT = vecStimTime_Real(indData);
		vecV = vecStimIFR_Real(indData);
		
		indPeaks = vecPopEventTimes_Real>2050.5 & vecPopEventTimes_Real < 2051;
		vecP_T = vecPopEventTimes_Real(indPeaks);
		vecP_L = vecPopEventLocsIFR_Real(indPeaks)-find(indData,1)+1;
		vecP_V = vecV(vecP_L);
		
		%merge until prominence no longer drops from merges
		indKeepPeaks = false(size(vecP_T));
		indProcessed = false(size(vecP_T));
		matPeakDomain = nan(numel(vecP_T),3); %peak sample idx, start sample idx, stop sample idx
		matPeakDomain(:,1) = vecP_L; %peak sample idx
		matPeakDomain(:,2) = vecP_L; %start sample idx
		matPeakDomain(:,3) = vecP_L; %stop sample idx
		while ~all(indProcessed)
			%% merge peaks
			%find highest unprocessed peak
			dblHeight = max(vecP_V(~indProcessed));
			intPeak = find(vecP_V==dblHeight);
			intPeakLoc = vecP_L(intPeak);
			dblPeakT = vecP_T(intPeak);
			
			%find nearest peak outside domain
			vecLeftDist = intPeakLoc - matPeakDomain(:,1);
			vecLeftDist(indProcessed | vecLeftDist<=0)=inf;
			[intLeftPeakDist,intLeftPeak]=min(vecLeftDist);
			
			vecRightDist = matPeakDomain(:,1) - intPeakLoc;
			vecRightDist(indProcessed | vecRightDist<=0)=inf;
			[intRightPeakDist,intRightPeak]=min(vecRightDist);
			
			%calculate prominence of original peak
			%left domain
			intLeftLoc = matPeakDomain(intPeak,2);
			intLeftTrough = findtrough(vecV,intLeftLoc,-1);
			matPeakDomain(intPeak,2) = intLeftTrough; %i dont think this will work
			
			%right domain
			intRightLoc = matPeakDomain(intPeak,2);
			intRightTrough = findtrough(vecV,intRightLoc,+1);
			matPeakDomain(intPeak,3) = intRightTrough; %i dont think this will work
			
			%prominence
			dblTrough = max([vecV(intLeftTrough) vecV(intRightTrough)]);
			dblOrigProm = dblHeight-dblTrough;
			
			%calculate prominance of merged peak
			%choose left or right
			[dummy,intLeftRight] = min([intLeftPeakDist intRightPeakDist]);
			if intLeftRight == 1 %left
				intClosestPeak = intLeftPeak;
				%keep right trough
				
				%left domain
				intLeftLoc = matPeakDomain(intClosestPeak,2);
				intLeftTrough = findtrough(vecV,intLeftLoc,-1);
				
			else %right
				intClosestPeak = intRightPeak;
				%keep left trough
				
				%right domain
				intRightLoc = matPeakDomain(intClosestPeak,2);
				intRightTrough = findtrough(vecV,intRightLoc,+1);
				
			end
			%prominence
			dblTrough = max([vecV(intLeftTrough) vecV(intRightTrough)]);
			dblNewProm = dblHeight-dblTrough;
			
			
			%decide to merge or not
			if dblOrigProm >= dblNewProm
				%complete
				indProcessed(intPeak) = true;
				indKeepPeaks(intPeak) = true;
			else
				%merge and continue
				indProcessed(intClosestPeak) = true;
				if intLeftRight == 1 %left
					matPeakDomain(intPeak,2) = vecP_L(intClosestPeak);
				else %right
					matPeakDomain(intPeak,3) = vecP_L(intClosestPeak);
				end
			end
			%%
			%error change to first all left merges until complete, then all right merges until complete
		end
		
		%create new list of merged peaks
		vecPeakLocs = vecP_L(indKeepPeaks);
		
		% plot
		figure
		plot(vecT,vecV);
		hold on
		scatter(vecP_T,vecP_V,'ro')
		scatter(vecP_T(indKeepPeaks),vecP_V(indKeepPeaks),'xb')
		return
	end
end
toc
