function [indTuned,vecNeuronPrefStim,vecOSI,cellSes] = getTunedStimDetectionNeurons(cellSes,matRespMat,boolOnlyPresence)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	if ~exist('boolOnlyPresence','var') || isempty(boolOnlyPresence)
		boolOnlyPresence = false;
	end
	%%
	vecOSI = [];
	if numel(cellSes) > 1
		%% prepare data
		vecOrientations = unique(cellSes{1}.structStim.Orientation);
		intOris = length(vecOrientations);
		intSessions = numel(cellSes);
		intNeurons= numel(cellSes{1}.neuron);
		indPresenceList = true(1,intNeurons);
		indResponseList = false(1,intNeurons);
		for intSes=1:intSessions
			for intNeuron=1:intNeurons
				if strcmp(cellSes{intSes}.neuron(intNeuron).strPresence,'absent')
					indPresenceList(intNeuron) = false;
				end
				if ~strcmp(cellSes{intSes}.neuron(intNeuron).strRespType,'silent') || boolOnlyPresence
					indResponseList(intNeuron) = true;
				end
			end
		end
		intNeurons = sum(indPresenceList&indResponseList);
		matPref = nan(intSessions,intNeurons);
		matOriResp = nan(intSessions*2,intNeurons,intOris);
		intCounter = 0;
		for intSes=1:intSessions
			%remove neurons
			cellSes{intSes}.neuron = cellSes{intSes}.neuron(indPresenceList&indResponseList);
			
			%combine opposite directions
			cellSes{intSes}.structStim.Orientation = mod(cellSes{intSes}.structStim.Orientation,180);
			
			
			%prepare data
			structStim = cellSes{intSes}.structStim;
			structStim.FrameOff = structStim.FrameOn+2;
			[matTrialResponseT,cellSelectContrasts] = getTrialResponseData(cellSes{intSes},structStim);
			intContrasts = length(cellSelectContrasts);
			
			%get orientation-based trial selection vectors
			sTypesOri = getStimulusTypes(structStim,{'Orientation'});
			cellSelectOri = getSelectionVectors(structStim,sTypesOri);
			
			%calculate preferred stimulus orientation for all neurons
			vecOriPrefContrasts = 4:6; % 8% - 100%
			
			for intContrast = vecOriPrefContrasts
				intCounter = intCounter + 1;
				matResp = matTrialResponseT(:,cellSelectContrasts{intContrast});
				structStimC{intContrast} = structStim;
				cellFields = fieldnames(structStim);
				for intField=1:length(cellFields)
					strField = cellFields{intField};
					structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
				end
				cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
				sTuning = calcTuningRespMat(matResp,cellSelect,vecOrientations);
				matPref(intCounter,:) = sTuning.vecPrefIndex;
				for intOri=1:intOris
					matOriResp(((intSes*2-1):intSes*2),:,intOri) = matResp(:,cellSelect{intOri})';
				end
			end
		end
		
		matOriPrefCount = nan(intOris,intNeurons);
		for intOri=1:intOris
			matOriPrefCount(intOri,:) = sum(matPref == intOri,1);
		end
		%check for majority vote
		[vecMax,vecNeuronPrefStim] = max(matOriPrefCount,[],1);
		indTuned = sum(bsxfun(@eq,matOriPrefCount,vecMax),1)==1;
		
		%calc OSI
		vecOrth = mod(vecNeuronPrefStim+round(intOris/2),intOris);
		vecOrth(vecOrth==0)=intOri;
		vecMeanPrefResp = nan(1,intNeurons);
		vecMeanOrthResp = nan(1,intNeurons);
		for intOri=1:intOris
			vecMeanOriResp = mean(matOriResp(:,:,intOri),1);
			vecMeanPrefResp(vecNeuronPrefStim==intOri) = vecMeanOriResp(vecNeuronPrefStim==intOri);
			vecMeanOrthResp(vecOrth==intOri) = vecMeanOriResp(vecOrth==intOri);
		end
		vecMeanPrefResp(vecMeanPrefResp < vecMeanOrthResp | vecMeanPrefResp < 0) = nan;
		vecMeanOrthResp(vecMeanOrthResp<0) = 0;
		vecOSI = (vecMeanPrefResp - vecMeanOrthResp)./(vecMeanPrefResp + vecMeanOrthResp);
		indTuned = ~isnan(vecOSI) & indTuned;
		
	else
		%take opposite directions as the same
		cellSes.structStim.Orientation = mod(cellSes.structStim.Orientation,180);
		vecOrientations = unique(cellSes.structStim.Orientation);
		intNeurons = numel(cellSes.neuron);
		
		%get neuronal responses per trial
		if exist('matRespMat','var') && ~isempty(matRespMat)
			cellFieldsC_SA{1} = 'Contrast';
			sTypesC_SA = getStimulusTypes(cellSes.structStim,cellFieldsC_SA);
			cellSelectContrasts = getSelectionVectors(cellSes.structStim,sTypesC_SA);
			matTrialResponse = matRespMat;
		else
			[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellSes,cellSes.structStim);
		end
		intContrasts = length(cellSelectContrasts);
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(cellSes,{'Orientation'});
		cellSelectOri = getSelectionVectors(cellSes.structStim,sTypesOri);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		%calculate preferred stimulus orientation for all neurons
		vecOriPrefContrasts = 4:6; % 8% - 100%
		matPref = nan(length(vecOriPrefContrasts),intNeurons);
		intCounter = 0;
		for intContrast = vecOriPrefContrasts
			intCounter = intCounter + 1;
			matResp = matTrialNormResponse(:,cellSelectContrasts{intContrast});
			structStimC{intContrast} = cellSes.structStim;
			cellFields = fieldnames(cellSes.structStim);
			for intField=1:length(cellFields)
				strField = cellFields{intField};
				structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
			end
			cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
			sTuning{intContrast} = calcTuningRespMat(matResp,cellSelect,vecOrientations);
			matPref(intCounter,:) = sTuning{intContrast}.vecPrefIndex;
		end
		vecNeuronPrefStim = nan(1,intNeurons);
		for intOri=1:length(vecOrientations);
			vecNeuronPrefStim(sum(matPref == intOri,1) > floor(length(vecOriPrefContrasts)/2)) = intOri;
		end
		%vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = randi(length(vecOrientations),[1 sum(isnan(vecNeuronPrefStim))]);
		indTuned = ~isnan(vecNeuronPrefStim);
		fprintf('%d neurons total; %d tuned [%.1f%%]\n',length(vecNeuronPrefStim),sum(~isnan(vecNeuronPrefStim)),(sum(~isnan(vecNeuronPrefStim))/length(vecNeuronPrefStim))*100);
	end
	%vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = randi(length(vecOrientations),[1 sum(isnan(vecNeuronPrefStim))]);
		
end

