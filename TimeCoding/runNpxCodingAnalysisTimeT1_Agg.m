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
		vecAllSpikeTime_Shuff = nan(1,intTotS);
		vecAllSpikeNeuron_Shuff = zeros(1,intTotS,'int16');
		intNumN = size(matMeanRate,1);
		intS = 1;
		for intN=1:intNumN
			%add spikes
			intThisS = numel(cellSpikeTimes{intN});
			vecSpikeT = cellSpikeTimes{intN};
			vecAllSpikeTime(intS:(intS+intThisS-1)) = vecSpikeT;
			vecAllSpikeNeuron(intS:(intS+intThisS-1)) = intN;
			
			%
			vecISI = diff(sort(vecSpikeT));
			vecGenSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecSpikeT(1);
			vecAllSpikeTime_Shuff(intS:(intS+intThisS-1)) = vecGenSpikes;
			vecAllSpikeNeuron_Shuff(intS:(intS+intThisS-1)) = intN;
		end
		%events
		vecEvents = vecOrigStimOnTime(1)-10;
		dblMaxDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur+10;
		
		%sort
		[vecAllSpikeTime,vecReorder] = sort(vecAllSpikeTime);
		vecAllSpikeNeuron = vecAllSpikeNeuron(vecReorder);
		[vecTime_Real,vecIFR_Real] = getIFR(vecAllSpikeTime,vecEvents,dblMaxDur,[],[],[],0); %takes about 1 minute
		vecTime_Real = vecTime_Real + vecEvents(1);
		
		%shuffled
		[vecAllSpikeTime_Shuff,vecReorder] = sort(vecAllSpikeTime_Shuff);
		vecAllSpikeNeuron_Shuff = vecAllSpikeNeuron_Shuff(vecReorder);
		[vecTime_Shuff,vecIFR_Shuff] = getIFR(vecAllSpikeTime_Shuff,vecEvents,dblMaxDur,[],[],[],0); %takes about 1 minute
		vecTime_Shuff = vecTime_Shuff + vecEvents(1);
		
		%filter
		intLag = round((numel(vecIFR_Real)/numel(vecOrigStimOnTime))/2);
		dblThreshZ = 1.5;
		dblInfluence = 0.05;
		[signals,avgFilter_Real,stdFilter_Real] = detectpeaks(vecIFR_Real,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Real = (vecIFR_Real - avgFilter_Real)./stdFilter_Real;
		vecIFR_Real = (vecIFR_Real - avgFilter_Real)./avgFilter_Real;
		indRem_Real = isinf(vecIFR_Real);
		vecIFR_Real(indRem_Real) = [];
		vecTime_Real(indRem_Real) = [];
		
		%shuff
		[signals,avgFilter_Shuff,stdFilter_Shuff] = detectpeaks(vecIFR_Shuff,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./stdFilter_Shuff;
		vecIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./avgFilter_Shuff;
		indRem_Shuff = isinf(vecIFR_Shuff);
		vecIFR_Shuff(indRem_Shuff) = [];
		vecTime_Shuff(indRem_Shuff) = [];
		
		%% go through trials
		for intType=[]%1:2
			if intType==1
				vecTime = vecTime_Real;
				vecIFR = vecIFR_Real;
			else
				vecTime = vecTime_Shuff;
				vecIFR = vecIFR_Shuff;
			end
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
			vecBinsE = -1.2:0.1:2;
			vecBinsC = vecBinsE(2:end)-median(diff(vecBinsE))/2;
			vecC = histcounts(vecIFR(vecTime>(vecOrigStimOnTime(1)) & vecTime<(vecOrigStimOnTime(end)+dblStimDur)),vecBinsE);
			plot(vecBinsC,vecC);
			xlabel('Raw IFR');
			ylabel('# of spikes');
			if intType==1
				vecHist_Real = vecC;
			else
				vecHist_Shuff = vecC;
			end
			
		end
		%figure
		%subplot(2,3,1)
		%plot(vecBinsC,vecHist_Real)
		%hold on
		%plot(vecBinsC,vecHist_Shuff)
		%subplot(2,3,4)
		%plot(vecBinsC,(vecHist_Real-vecHist_Shuff)./(vecHist_Real+vecHist_Shuff))
		
		%% shuff peaks
		indStimSpikes_Shuff = vecTime_Shuff>(vecOrigStimOnTime(1)-dblStimDur) & vecTime_Shuff<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Shuff = vecIFR_Shuff(indStimSpikes_Shuff);
		vecStimTime_Shuff = vecTime_Shuff(indStimSpikes_Shuff);
			
		[pks_Shuff,locs_Shuff,w_Shuff,p_Shuff] = findpeaks(vecStimIFR_Shuff);
			
		%% real peaks
		indStimSpikes_Real = vecTime_Real>(vecOrigStimOnTime(1)-dblStimDur) & vecTime_Real<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Real = vecIFR_Real(indStimSpikes_Real);
		vecStimTime_Real = vecTime_Real(indStimSpikes_Real);
			
		[pks_Real,locs_Real,w_Real,p_Real] = findpeaks(vecStimIFR_Real);
		
		%% plot
		vecBinWE = 0:20;
		vecWC_Real = histcounts(w_Real,vecBinWE);
		vecWC_Shuff = histcounts(w_Shuff,vecBinWE);
		vecBinWC = vecBinWE(2:end)-0.5;
		
		%width
		figure
		subplot(2,3,1)
		plot(vecBinWC,vecWC_Real)
		hold on
		plot(vecBinWC,vecWC_Shuff)
		
		subplot(2,3,4)
		plot(vecBinWC,(vecWC_Real-vecWC_Shuff)./(vecWC_Real+vecWC_Shuff));
		
		%prom
		vecBinPP = 0:0.1:2;
		vecPC_Real = histcounts(p_Real,vecBinPP);
		vecPC_Shuff = histcounts(p_Shuff,vecBinPP);
		vecBinPC = vecBinPP(2:end)-0.5;
		
		subplot(2,3,2)
		plot(vecBinPC,vecPC_Real)
		hold on
		plot(vecBinPC,vecPC_Shuff)
		
		subplot(2,3,5)
		plot(vecBinPC,(vecPC_Real-vecPC_Shuff)./(vecPC_Real+vecPC_Shuff));
		
		
		%rate
		subplot(2,3,3);cla
		vecBinEI = 0:0.1:2;
		vecBinPI = vecBinEI(2:end) + mean(diff(vecBinEI))/2;
		vecIC_Real = histcounts(vecStimIFR_Real,vecBinEI);
		vecIC_Shuff = histcounts(vecStimIFR_Shuff,vecBinEI);
		
		plot(vecBinPI,vecIC_Real)
		hold on
		plot(vecBinPI,vecIC_Shuff)
		
		subplot(2,3,6);cla
		plot(vecBinPI,(vecIC_Real-vecIC_Shuff)./(vecIC_Real+vecIC_Shuff));
		
		%height
		vecBinEH = -0.5:0.1:2;
		vecBinPH = vecBinEH(2:end) + mean(diff(vecBinEH))/2;
		vecHC_Real = histcounts(pks_Real,vecBinEH);
		vecHC_Real = vecHC_Real./sum(vecHC_Real(:));
		vecHC_Shuff = histcounts(pks_Shuff,vecBinEH);
		vecHC_Shuff = vecHC_Shuff./sum(vecHC_Shuff(:));
		intRemIdx = find(vecBinPH>0 & (vecHC_Real==0 | vecHC_Shuff==0),1);
		if ~isempty(intRemIdx)
			indRem = false(size(vecHC_Real));
			indRem(intRemIdx:end) = true;
			vecHC_Real(indRem) = [];
			vecHC_Shuff(indRem) = [];
			vecBinPH(indRem) = [];
		end
		%consider only positive peaks
		pks_RealPos = pks_Real(pks_Real>0);
		pks_ShuffPos = pks_Shuff(pks_Shuff>0);
		vecRatio = (vecHC_Real-vecHC_Shuff)./(vecHC_Real+vecHC_Shuff);
		intCutOff = find(vecBinPH>0 & (vecRatio>0 & [vecRatio(2:end)>0 false]),1);
		if isempty(intCutOff) || mean(pks_RealPos) < mean(pks_ShuffPos)
			strTit = sprintf('No cut-off; pos mu real=%.3f,shuff=%.3f',mean(pks_RealPos) ,mean(pks_ShuffPos));
		else
			dblCutOff = vecBinPH(intCutOff);
			strTit = sprintf('Cut-off = %.2f; including %.0f%%',dblCutOff,(sum(pks_Real>dblCutOff)/numel(pks_Real))*100);
		end
		
		
		subplot(2,3,3);cla
		plot(vecBinPH,vecHC_Real)
		hold on
		plot(vecBinPH,vecHC_Shuff)
		title(strRec,'interpreter','none');
		
		subplot(2,3,6);cla
		plot(vecBinPH,vecRatio);
		title(strTit,'interpreter','none');
		
		continue;
		%% scatter properties of peaks
		%figure
		%scatter(pks_Real,p_Real)
		%hold on
		%scatter(pks_Shuff,p_Shuff)
		
		
		%% algo #2
		%vecStimIFR_Shuff = vecIFR_Shuff(indStimSpikes_Shuff);
		%vecStimTime_Shuff = vecTime_Shuff(indStimSpikes_Shuff);
		
		%real
		intLag = round((numel(vecStimIFR_Real)/numel(vecOrigStimOnTime))/2);
			dblThreshZ = 1.5;
			dblInfluence = 0.05;
			[signals,avgFilter_Real,stdFilter_Real] = detectpeaks(vecStimIFR_Real,intLag,dblThreshZ,dblInfluence);
			vecPeakLocs = find(signals==1);
			figure
			subplot(2,3,1)
			plot(vecStimIFR_Real)
			
			subplot(2,3,2)
			scatter(vecPeakLocs,vecStimIFR_Real(vecPeakLocs),'gx')
			hold on
			plot(avgFilter_Real)
			plot(avgFilter_Real+stdFilter_Real)
			
			subplot(2,3,3)
			vecISN_Real = (vecStimIFR_Real - avgFilter_Real)./stdFilter_Real;
			
			plot(vecISN_Real)
			
			%shuff
			intLag = round((numel(vecStimIFR_Shuff)/numel(vecOrigStimOnTime))/2);
			[signals,avgFilter_Shuff,stdFilter_Shuff] = detectpeaks(vecStimIFR_Shuff,intLag,dblThreshZ,dblInfluence);
			vecPeakLocs = find(signals==1);
			
			subplot(2,3,4)
			plot(vecStimIFR_Shuff)
			
			subplot(2,3,5)
			scatter(vecPeakLocs,vecStimIFR_Shuff(vecPeakLocs),'gx')
			hold on
			plot(avgFilter_Shuff)
			plot(avgFilter_Shuff+stdFilter_Shuff)
			
			subplot(2,3,6)
			vecISN_Shuff = (vecStimIFR_Shuff - avgFilter_Shuff)./stdFilter_Shuff;
			
			plot(vecISN_Shuff)
		return
	end
end
toc
