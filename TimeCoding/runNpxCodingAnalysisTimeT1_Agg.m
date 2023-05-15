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
		
		%%
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		intTotS = sum(cellfun(@numel,cellSpikeTimes));
		vecAllSpikeTime = nan(1,intTotS);
		vecAllSpikeNeuron = zeros(1,intTotS,'int16');
		intNumN = size(matMeanRate,1);
		intS = 1;
		for intN=1:intNumN
			%add spikes
			intThisS = numel(cellSpikeTimes{intN});
			vecAllSpikeTime(intS:(intS+intThisS-1)) = cellSpikeTimes{intN};
			vecAllSpikeNeuron(intS:(intS+intThisS-1)) = intN;
		end
		%sort
		[vecAllSpikeTime,vecReorder] = sort(vecAllSpikeTime);
		vecAllSpikeNeuron = vecAllSpikeNeuron(vecReorder);
		vecEvents = vecOrigStimOnTime(1)-10;
		dblMaxDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur+10;
		[vecTime,vecIFR] = getIFR(vecAllSpikeTime,vecEvents,dblMaxDur,[],[],[],0); %takes about 1 minute
		vecTime = vecTime + vecEvents(1);
		
		%% go through trials
		figure;maxfig;
		subplot(2,3,1)
		hold on
		dblBinDur = 0.05;
		vecBinEdges = 0:dblBinDur:dblStimDur;
		vecBinC = vecBinEdges(2:end)-dblBinDur/2;
		intBinNum = numel(vecBinC);
		matAvgIFR = nan(intStimNum,intBinNum);
		[vecMean,vecSEM,vecWindowBinCenters,matAvgPET] = doPEP(vecTime,vecBinEdges,vecOrigStimOnTime,-1);
		for intRep=1:intStimNum
			indSpikesInTrial = vecTime>(vecOrigStimOnTime(intRep)-10) & vecTime<(vecOrigStimOnTime(intRep)+dblStimDur+10);
			vecTrialTime = vecTime(indSpikesInTrial)-vecOrigStimOnTime(intRep);
			vecTrialIFR = vecIFR(indSpikesInTrial);
			try
				matAvgIFR(intRep,:) = interp1(vecTrialTime,vecTrialIFR,vecBinC);
			end
			if intRep==50
			plot(vecTrialTime,vecTrialIFR,'color',[0.5 0.5 0.5]);
			end
		end
		xlim([0 dblStimDur]);
		title(strRec,'interpreter','none');
		xlabel('Time (s)');
		ylabel('Raw IFR of trial #50');
		
		subplot(2,3,2)
		errorbar(vecBinC,vecMean,vecSEM);
		xlabel('Time (s)');
		ylabel('Avg binned spiking rate');
		
		subplot(2,3,3)
		errorbar(vecBinC,nanmean(matAvgIFR,1),nanstd(matAvgIFR,[],1)./sqrt(intStimNum));
		xlabel('Time (s)');
		ylabel('Avg binned IFR');
		
		% make distro plots
		subplot(2,3,4)
		histx(matAvgIFR(:));
		xlabel('Binned IFR');
		ylabel('# of bins');
		
		subplot(2,3,5)
		histx(matAvgPET(:));
		xlabel('Binned spike rates');
		ylabel('# of bins');
		
		subplot(2,3,6)
		histx(vecIFR(vecTime>(vecOrigStimOnTime(1)) & vecTime<(vecOrigStimOnTime(end)+dblStimDur)));
		xlabel('Raw IFR');
		ylabel('# of spikes');
		
		%% detect peaks
		indStimSpikes = vecTime>(vecOrigStimOnTime(1)-dblStimDur) & vecTime<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR = vecIFR(indStimSpikes);
		vecStimTime = vecTime(indStimSpikes);
		
		vecT = vecTime(1:2000);
		vecI = vecIFR(1:2000);
		[pks,locs,w,p] = findpeaks(vecI);
		figure
		plot(vecI);
		hold on
		scatter(locs,pks,'rx')
		tsiqr = iqr(vecI)
		
		intLag = round((numel(vecStimIFR)/numel(vecOrigStimOnTime))/2);
		dblThreshZ = 1.5;
		dblInfluence = 0;
		[signals,avgFilter,stdFilter] = detectpeaks(vecI,intLag,dblThreshZ,dblInfluence);
		vecPeakLocs = find(signals==1);
		scatter(vecPeakLocs,vecI(vecPeakLocs),'gx')
		
	end
end
toc
