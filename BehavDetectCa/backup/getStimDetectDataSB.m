function [sContrast,sMetaData] = getStimDetectDataSB(ses,sStimAggregate)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get common variables
	intNeurons = numel(ses.neuron);
	dblBaselineSecs = 0.1;
	
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
	
	%transform structure-based data to raw dFoF matrix
	intFrames = length(ses.neuron(1).dFoF);
	matActivity = zeros(intNeurons,intFrames);
	for intNeuron=1:intNeurons
		matActivity(intNeuron,:) = ses.neuron(intNeuron).dFoF;
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
	
	
	%check for stim aggregate; if present, then calculate mean number of
	%frames for detected stimuli
	if exist('sStimAggregate','var') && isstruct(sStimAggregate)
		cellFieldsC_SA{1} = 'Contrast';
		sTypesC_SA = getStimulusTypes(ses,cellFieldsC_SA);
		cellSelectC_SA = getSelectionVectors(sStimAggregate,sTypesC_SA);
		
		vecMeanStimDur = nan(1,intNumContrasts);
		%get responded trials
		vecResponded = logical(sStimAggregate.vecTrialResponse);
		for intContrastIndex=1:intNumContrasts
			%get target trials
			vecSelectTrials = vecResponded & cellSelectC_SA{intContrastIndex};
			
			%calculate mean stim duration
			vecMeanStimDur(intContrastIndex) = round(mean(sStimAggregate.FrameOff(vecSelectTrials)-sStimAggregate.FrameOn(vecSelectTrials)));
		end
	else
		vecMeanStimDur = nan(1,intNumContrasts);
	end
	
	%create structures
	sContrast = struct;
	sMetaData = struct;
	
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		vecContrasts(intContrastIndex) = dblContrast;
		vecCorrRespC = logical(ses.structStim.vecTrialResponse) & vecSelectC;
		
		%pre-allocate output matrices
		intPreAllocRespN = intNeurons * sum(vecCorrRespC);
		intPreAllocNoRespN = intNeurons * sum(~vecCorrRespC);
		vecPR = nan(1,intPreAllocRespN);
		vecNPR = nan(1,intPreAllocRespN);
		vecDiffR = nan(1,intPreAllocRespN);
		vecPN = nan(1,intPreAllocNoRespN);
		vecNPN = nan(1,intPreAllocNoRespN);
		vecDiffN = nan(1,intPreAllocNoRespN);
		intCounterPR = 1;
		intCounterNPR = 1;
		intCounterDiffR = 1;
		intCounterPN = 1;
		intCounterNPN = 1;
		intCounterDiffN = 1;
		
		%plot per type
		for intStimType=1:intTypes
			vecPrefNeurons = vecPrefStim==intStimType;
			intPrefNeurons = sum(vecPrefNeurons);
			vecNonPrefNeurons = ~vecPrefNeurons;
			intNonPrefNeurons = sum(vecNonPrefNeurons);
			if intPrefNeurons>0
				
				%get mean stim dur
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				if isnan(vecMeanStimDur(intContrastIndex))
					vecNoRespOff = ses.structStim.FrameOff(vecNoRespTrials);
				else
					vecNoRespOff = ses.structStim.FrameOn(vecNoRespTrials) + vecMeanStimDur(intContrastIndex);
				end
								
				%correct trials
				%get selection vectors
				vecRespTrials = structStimCorrs.cellSelect{intStimType} & vecCorrResp & vecSelectC;
				vecRespOn = ses.structStim.FrameOn(vecRespTrials);
				vecRespOff = ses.structStim.FrameOff(vecRespTrials);
				vecBaseStart = vecRespOn - round(ses.samplingFreq*dblBaselineSecs) - 1;
				vecBaseStop = vecRespOn - 1;
				
				for intTrial=1:length(vecRespOn)
					%pref + corr
					dblStimAct = std(mean(matActivity(vecPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2));
					vecPrefAct = mean(matActivity(vecPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2);
					vecPrefBaseAct = mean(matActivity(vecPrefNeurons,vecBaseStart(intTrial):vecBaseStop(intTrial)),2);
					dblPrefAct = dblStimAct;% - dblBaseAct;
					
					vecPR(intCounterPR:(intCounterPR+intPrefNeurons-1)) = vecPrefAct(:) - vecPrefBaseAct(:);
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					dblStimAct = std(mean(matActivity(vecNonPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2));
					vecNonPrefAct = mean(matActivity(vecNonPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2);
					vecNonPrefBaseAct = mean(matActivity(vecNonPrefNeurons,vecBaseStart(intTrial):vecBaseStop(intTrial)),2);
					dblNonPrefAct = dblStimAct;% - dblBaseAct;
					
					vecNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1)) = vecNonPrefAct(:) - vecNonPrefBaseAct(:);
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
					
					%difference
					vecDiffR(intCounterDiffR) = dblPrefAct;% - dblNonPrefAct;
					intCounterDiffR = intCounterDiffR + 1;
				end
				
				
				%incorrect trials
				%get selection vector
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				vecNoRespOn = ses.structStim.FrameOn(vecNoRespTrials);
				
				vecBaseStart = vecNoRespOn - round(ses.samplingFreq*dblBaselineSecs) - 1;
				vecBaseStop = vecNoRespOn - 1;
				

				
				for intTrial=1:length(vecNoRespOn)
					%pref + incorr
					%matPN(intCounterPN:(intCounterPN+intPrefNeurons-1),:) = std(mean(matActivity(vecPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2),[],1);
					%dblStimAct = mean(std(matActivity(vecPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),[],2));
					dblStimAct = std(mean(matActivity(vecPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2));
					vecPrefAct = mean(matActivity(vecPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2);
					vecPrefBaseAct = mean(matActivity(vecPrefNeurons,vecBaseStart(intTrial):vecBaseStop(intTrial)),2);
					dblPrefAct = dblStimAct;% - dblBaseAct;
					
					vecPN(intCounterPN:(intCounterPN+intPrefNeurons-1)) = vecPrefAct(:) - vecPrefBaseAct(:);
					intCounterPN = intCounterPN + intPrefNeurons;
					
					%non-pref + incorr
					%matNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1),:) = std(mean(matActivity(vecNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2),[],1);
					%dblStimAct = mean(std(matActivity(vecNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),[],2));
					dblStimAct = std(mean(matActivity(vecNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2));
					vecNonPrefAct = mean(matActivity(vecNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2);
					vecNonPrefBaseAct = mean(matActivity(vecNonPrefNeurons,vecBaseStart(intTrial):vecBaseStop(intTrial)),2);
					dblNonPrefAct = dblStimAct;% - dblBaseAct;
					
					vecNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1)) = vecNonPrefAct(:) - vecNonPrefBaseAct(:);
					intCounterNPN = intCounterNPN + intNonPrefNeurons;
					
					%difference
					vecDiffN(intCounterDiffN) = dblPrefAct;% - dblNonPrefAct;
					intCounterDiffN = intCounterDiffN + 1;
				end
			end
		end
		
		%remove trailing nans
		[indNan]=find(isnan(vecPR),1);
		vecPR = vecPR(1:(indNan-1));
		[indNan]=find(isnan(vecNPR),1);
		vecNPR = vecNPR(1:(indNan-1));
		[indNan]=find(isnan(vecPN),1);
		vecPN = vecPN(1:(indNan-1));
		[indNan]=find(isnan(vecNPN),1);
		vecNPN = vecNPN(1:(indNan-1));
		[indNan]=find(isnan(vecDiffR),1);
		vecDiffR = vecDiffR(1:(indNan-1));
		[indNan]=find(isnan(vecDiffN),1);
		vecDiffN = vecDiffN(1:(indNan-1));
		
		%{
		[row,col]=find(isnan(matPR_mean));
		matPR_mean = matPR_mean(1:(row-1),:);
		[row,col]=find(isnan(matNPR_mean));
		matNPR_mean = matNPR_mean(1:(row-1),:);
		[row,col]=find(isnan(matPN_mean));
		matPN_mean = matPN_mean(1:(row-1),:);
		[row,col]=find(isnan(matNPN_mean));
		matNPN_mean = matNPN_mean(1:(row-1),:);
		%}
		%add data to structure
		sContrast(intContrastIndex).vecPR = vecPR;
		sContrast(intContrastIndex).vecNPR = vecNPR;
		sContrast(intContrastIndex).vecPN = vecPN;
		sContrast(intContrastIndex).vecNPN = vecNPN;
		sContrast(intContrastIndex).vecDiffR = vecDiffR;
		sContrast(intContrastIndex).vecDiffN = vecDiffN;
		%{
		sContrast(intContrastIndex).matPR_mean = matPR_mean;
		sContrast(intContrastIndex).matNPR_mean = matNPR_mean;
		sContrast(intContrastIndex).matPN_mean = matPN_mean;
		sContrast(intContrastIndex).matNPN_mean = matNPN_mean;
		%}
	end
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
end

