function [sContrast,sMetaData,matPowerSpectrum] = getStimDetectDataAP(ses,sStimAggregate)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get common variables
	intNeurons = numel(ses.neuron);
	
	% get indexing vectors for unique stimulus combinations
	sTypes = getStimulusTypes(ses);
	cellSelect = getSelectionVectors(ses.structStim,sTypes);
	
	%remove empty types
	intTypes = length(cellSelect);
	vecKeep = true(1,intTypes);
	for intType = 1:intTypes
		if sum(cellSelect{intType}) == 0
			vecKeep(intType) = false;
		end
	end
	cellSelect = cellSelect(vecKeep);
	sTypes = getStimulusTypes(ses);
	sTypes.matTypes = sTypes.matTypes(:,vecKeep);
	
	%transform session data to neuropil dF/F
	strFilename = [ses.session 'xyt' sprintf('%02d',ses.recording) '_sesnp'];
	if exist([ses.strImPath strFilename '.mat'],'file')
		fprintf('Found and loaded neuropil annulus ses file [%s]...\n',strFilename);
		sLoad = load([ses.strImPath strFilename]);
		ses_np = sLoad.ses;
	else
		%remove neurons with incorrect types
		indKeepList = true(1,numel(ses.neuron));
		for intNeuron=1:numel(ses.neuron)
			if strcmp(ses.neuron(intNeuron).strPresence,'absent') || strcmp(ses.neuron(intNeuron).strRespType,'silent')
				indKeepList(intNeuron) = false;
			end
		end
		[ses,dummy] = doRecalcdFoF(ses,6,indKeepList);
		
		%calculate neuropil dF/F
		ses_old = ses;
		fprintf('Calculating dF/F for neuropil annuli...\n');
		ses = doRecalcdFoF(ses,4);
		fprintf('\b    Done!\n\n');
		save([ses.strImPath strFilename],'ses');
		ses_np = ses;
		ses = ses_old;
	end
	
	%transform structure-based data to raw dFoF matrix
	intFrames = length(ses_np.neuron(1).dFoF);
	matF_np = zeros(intNeurons,intFrames);
	for intNeuron=1:intNeurons
		matF_np(intNeuron,:) = ses_np.neuron(intNeuron).dFoF;
	end
	
	%perform FFT for all neurons on each trial
	dblSamplingFreq = ses_np.samplingFreq;
	dblT = 1; %seconds
	intPowerFrames = round(dblSamplingFreq*dblT);
	intNumTrials = length(ses_np.structStim.FrameOn);
	intNFFT = 2^nextpow2(intPowerFrames); % Next power of 2 from length of signal
	intNumFrequencies = intNFFT;
	matPowerSpectrum = nan(intNeurons,intNumTrials,intNumFrequencies);
	for intTrial=1:intNumTrials
		
		%get frames for this trial
		intFrameStop = ses_np.structStim.FrameOn(intTrial) - 1;
		intFrameStart = intFrameStop - intPowerFrames + 1;
		
		%get raw data
		matSignal = matF_np(:,intFrameStart:intFrameStop);
		
		%calculate FFT
		matPowerTrial = abs(fft(matSignal,intNFFT,2)/intPowerFrames);
		
		%put in output
		matPowerSpectrum(:,intTrial,:) = matPowerTrial;
	end
	
	
	%general selection vectors
	vecCorrResp = logical(ses.structStim.vecTrialResponse);
	cellFieldsC{1} = 'Contrast';
	sTypesC = getStimulusTypes(ses,cellFieldsC);
	cellSelectC = getSelectionVectors(ses.structStim,sTypesC);
	vecContrasts = nan(1,length(cellSelectC));
	intNumContrasts = length(vecContrasts);
	
	
	%select only full contrast stimuli
	vecTrialsFullC = cellSelectC{end} | cellSelectC{end-1};
	cellSubFields = fieldnames(ses.structStim);
	sesFullC = ses;
	for intField=1:length(cellSubFields)
		strField = cellSubFields{intField};
		sesFullC.structStim.(strField) = sesFullC.structStim.(strField)(vecTrialsFullC);
	end
	
	%get pref stim per neuron
	structStimCorrsFullC = calcStimCorrs(sesFullC);
	matSignalResponse = structStimCorrsFullC.matSignalResponse;
	[vecStimAct,vecPrefStim]=max(matSignalResponse,[],1);
	structStimCorrs = calcStimCorrs(ses);
	
	%create structures
	sContrast = struct;
	sMetaData = struct;
	
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		vecContrasts(intContrastIndex) = dblContrast;
		
		%pre-allocate output matrices
		intPreAllocRespN = intNeurons * sum(vecCorrResp);
		intPreAllocNoRespN = intNeurons * sum(~vecCorrResp);
		matPR = nan(intPreAllocRespN,intNFFT);
		matNPR = nan(intPreAllocRespN,intNFFT);
		matPN = nan(intPreAllocNoRespN,intNFFT);
		matNPN = nan(intPreAllocNoRespN,intNFFT);
		intCounterPR = 1;
		intCounterNPR = 1;
		intCounterPN = 1;
		intCounterNPN = 1;
		
		%plot per type
		for intStimType=1:intTypes
			vecPrefNeurons = vecPrefStim==intStimType;
			intPrefNeurons = sum(vecPrefNeurons);
			vecNonPrefNeurons = ~vecPrefNeurons;
			intNonPrefNeurons = sum(vecNonPrefNeurons);
			if intPrefNeurons>0
				
				%correct trials
				%get selection vectors
				vecRespTrials = structStimCorrs.cellSelect{intStimType} & vecCorrResp & vecSelectC;
				for intTrial=find(vecRespTrials)
					%pref + corr
					matPR(intCounterPR:(intCounterPR+intPrefNeurons-1),:) = matPowerSpectrum(vecPrefNeurons,intTrial,:);
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					matNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1),:) = matPowerSpectrum(vecNonPrefNeurons,intTrial,:);
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
				end
				
				
				%incorrect trials
				%get selection vector
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				for intTrial=1:find(vecNoRespTrials)
					%pref + incorr
					matPN(intCounterPN:(intCounterPN+intPrefNeurons-1),:) = matPowerSpectrum(vecPrefNeurons,intTrial,:);
					intCounterPN = intCounterPN + intPrefNeurons;
					
					%non-pref + incorr
					matNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1),:) = matPowerSpectrum(vecNonPrefNeurons,intTrial,:);
					intCounterNPN = intCounterNPN + intNonPrefNeurons;
				end
			end
		end
		
		%remove trailing nans
		[row,col]=find(isnan(matPR));
		if ~isempty(row),matPR = matPR(1:(row-1),:);end
		[row,col]=find(isnan(matNPR));
		if ~isempty(row),matNPR = matNPR(1:(row-1),:);end
		[row,col]=find(isnan(matPN));
		if ~isempty(row),matPN = matPN(1:(row-1),:);end
		[row,col]=find(isnan(matNPN));
		if ~isempty(row),matNPN = matNPN(1:(row-1),:);end
		
		
		%add data to structure
		sContrast(intContrastIndex).matPR = matPR;
		sContrast(intContrastIndex).matNPR = matNPR;
		sContrast(intContrastIndex).matPN = matPN;
		sContrast(intContrastIndex).matNPN = matNPN;
	end
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
end

