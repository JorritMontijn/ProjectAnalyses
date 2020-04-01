function sTuning = calcTuningRespMat(matResponse,cellSelect,vecOrientations)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	%depends on:
	%ang2rad
	%circ_var
	%getStimulusTypes [depends on getUniqueVals.m]
	%getSelectionVectors
	%getNeuronResponse
	%imnorm
	
	%check input
	vecNeurons = 1:size(matResponse,1);
	if nargin < 3
		vecOrientations = 0:180/length(cellSelect):179;
	end
	if ~iscell(cellSelect)
		vecClasses = cellSelect;
		vecUniqueClasses = unique(vecClasses);
		intClasses = length(vecUniqueClasses);
		cellSelect = cell(1,intClasses);
		%loop through classes
		for intClass=1:intClasses
			dblClass = vecUniqueClasses(intClass);
			cellSelect{intClass} = vecClasses == dblClass;
		end
	end
	
	%check if ori or dir
	boolOri = (max(vecOrientations) - min(vecOrientations)) < 180;
	
	% loop through stimulus types to get nr of reps
	intMaxReps = 0;
	for intS=1:length(cellSelect)
		intMaxReps = max(intMaxReps,sum(cellSelect{intS}));
	end
	
	%pre-allocate response matrix
	intOris = length(vecOrientations);
	intNeuronMax = max(vecNeurons);
	matStimResponse = doMatRespTransform(matResponse,cellSelect);
	
	%pre-allocate output
	vecOSI = nan(1,intNeuronMax);
	vecDSI = nan(1,intNeuronMax);
	vecPrefIndex = nan(1,intNeuronMax);
	vecPrefAngle = nan(1,intNeuronMax);
	vecMeanPrefResp = nan(1,intNeuronMax);
	vecSDPrefResp = nan(1,intNeuronMax);
	vecMeanOrthResp = nan(1,intNeuronMax);
	vecSDOrthResp = nan(1,intNeuronMax);
	
	%transform to ori (if dir)
	if boolOri
		matOriResponse = matStimResponse;
		vecOris = 2*ang2rad(vecOrientations);
	else
		intDirs = intOris;
		intOris = intDirs/2;
		vecFirstHalf = 1:intOris;
		vecSecondHalf = (intOris+1):intDirs;
		matOriResponse = (matStimResponse(vecFirstHalf,:,:) + matStimResponse(vecSecondHalf,:,:))/2;
		
		vecOris = 2*ang2rad(vecOrientations(vecFirstHalf));
	end
	
	%get tuning properties
	for intNeuron = vecNeurons
		%get data and normalize
		matOriResp = matOriResponse(:,:,intNeuron);
		vecMeanResp = nanmean(matOriResp,2);
		vecStdResp = nanstd(matOriResp,[],2);
		
		%get pref stim
		[dblPrefResp,intPrefIndex] = max(vecMeanResp);
		dblPrefAngle = vecOrientations(intPrefIndex);
		vecPrefIndex(intNeuron) = intPrefIndex;
		vecPrefAngle(intNeuron) = dblPrefAngle;
		vecMeanPrefResp(intNeuron) = dblPrefResp;
		vecSDPrefResp(intNeuron) = vecStdResp(intPrefIndex);
		
		%get orth stim
		intOrth = mod(intPrefIndex+round(length(vecMeanResp)/2),length(vecMeanResp));
		if intOrth == 0,intOrth = length(vecMeanResp);end
		vecMeanOrthResp(intNeuron) = vecMeanResp(intOrth);
		vecSDOrthResp(intNeuron) = vecStdResp(intOrth);
		
		%get osi
		%dblCircVar = circ_var(vecOris', vecMeanResp);
		%vecOSI(intNeuron) = 1 - dblCircVar;
		
		dblOrthResp = vecMeanResp(intOrth);
		vecOSI(intNeuron) = (dblPrefResp - dblOrthResp)/(dblPrefResp + dblOrthResp);
		
		%get dsi
		if ~boolOri
			%get data and normalize
			matTempResp = matStimResponse(:,:,intNeuron);
			vecTempResp = nanmean(matTempResp,2);
			
			intOppIndex = mod(intPrefIndex-1+round(intDirs/2),intDirs)+1;
			dblOppResp = vecTempResp(intOppIndex);
			vecDSI(intNeuron) = 1 - abs(dblOppResp/dblPrefResp);
		end
	end
	
	%get Fano factors
	vecRatioMean = vecMeanPrefResp ./ vecMeanOrthResp;
	vecRatioVar = vecSDPrefResp.^2 ./ vecSDOrthResp.^2;
	vecFano = vecRatioVar ./ vecRatioMean;
	
	%assign to output
	sTuning.vecOSI = vecOSI;
	sTuning.vecDSI = vecDSI;
	sTuning.vecPrefIndex = vecPrefIndex;
	sTuning.vecPrefAngle = vecPrefAngle;
	sTuning.vecMeanPrefResp = vecMeanPrefResp;
	sTuning.vecSDPrefResp = vecSDPrefResp;
	sTuning.vecMeanOrthResp = vecMeanOrthResp;
	sTuning.vecSDOrthResp = vecSDOrthResp;
	sTuning.vecRatioMean = vecRatioMean;
	sTuning.vecRatioVar = vecRatioVar;
	sTuning.vecFano = vecFano;
end

