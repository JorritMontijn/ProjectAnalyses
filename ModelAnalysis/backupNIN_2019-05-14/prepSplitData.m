function [vecTrialStimTypeUnsplit,matRespUnsplit,vecTrialStimTypeSplit,matRespSplitNorm,matRespSplit,matSplitTrialIdx] = prepSplitData(matModelRespTS3,vecTrialStimType,intUseIndepTrials)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	
	%% transform model resp
	%get trial data; save original
	vecTrialStimTypeTrials = vecTrialStimType;
	
	%data
	matModelRespTS3 = double(matModelRespTS3); %[neurons * bins * trials]; data
	[intNeurons,intBinsPerTrial,intTrials] = size(matModelRespTS3);
	
	%build splits by stimtype
	vecUniqueStimTypes = unique(vecTrialStimType);
	intStimTypes = numel(vecUniqueStimTypes);
	if exist('intUseIndepTrials','var')
		intRepetitions = intUseIndepTrials;
	else
		intRepetitions = intTrials/intStimTypes;
	end
	
	matTrialIdxTS4 = nan(intNeurons,intBinsPerTrial,intRepetitions,intStimTypes);
	matModelRespTS4 = nan(intNeurons,intBinsPerTrial,intRepetitions,intStimTypes);
	matModelRespTS4Mean = nan(intNeurons,intBinsPerTrial,1,intStimTypes);
	matTrialTypeTS4 = nan(intNeurons,intBinsPerTrial,intRepetitions,intStimTypes);
	for intStimTypeIdx=1:intStimTypes
		intStimType = vecUniqueStimTypes(intStimTypeIdx);
		vecTheseTrials = find(vecTrialStimTypeTrials==intStimType);
		vecUseTrialSubIdx = randperm(numel(vecTheseTrials),intRepetitions);
		vecUseTheseTrials = vecTheseTrials(vecUseTrialSubIdx);
		%remove additional independent trials
		matThisData = matModelRespTS3(:,:,vecUseTheseTrials);
		matTrialIdxTS4(:,:,:,intStimTypeIdx) = repmat(reshape(vecUseTheseTrials,[1 1 intRepetitions]),[intNeurons intBinsPerTrial]);
		matModelRespTS4(:,:,:,intStimTypeIdx) = matThisData;
		matModelRespTS4Mean(:,:,1,intStimTypeIdx) = xmean(matThisData,3);
		matTrialTypeTS4(:,:,:,intStimTypeIdx) = intStimType;
	end
	
	%remove means
	matModelRespTS4Norm=nan(size(matModelRespTS4));
	for intStimTypeIdx=1:intStimTypes
		matModelRespTS4Norm(:,:,:,intStimTypeIdx) = bsxfun(@minus,matModelRespTS4(:,:,:,intStimTypeIdx),mean(matModelRespTS4(:,:,:,intStimTypeIdx),3));
	end
	
	%define data to be used
	intBinTrialStimCombs = intRepetitions*intBinsPerTrial*intStimTypes;
	matDataNorm = reshape(matModelRespTS4Norm,[intNeurons intBinTrialStimCombs]);
	matDataFull = reshape(matModelRespTS4,[intNeurons intBinTrialStimCombs]);
	matOrigTrials = reshape(matTrialIdxTS4,[intNeurons intBinTrialStimCombs]);
	vecTrialStimTypeTS = reshape(matTrialTypeTS4,[intNeurons intBinTrialStimCombs]);
	vecTrialStimType = vecTrialStimTypeTS(1,:);
	
	%split
	matRespSplit = matDataFull;
	matRespSplitNorm = matDataNorm;
	vecTrialStimTypeSplit = vecTrialStimType;
	%non-split
	matRespUnsplitStim = squeeze(sum(matModelRespTS4,2));
	matRespUnsplit = [];
	vecTrialStimTypeUnsplit = [];
	for intStim=1:size(matRespUnsplitStim,3)
		matRespUnsplit = cat(2,matRespUnsplit,matRespUnsplitStim(:,:,intStim));
		vecTrialStimTypeUnsplit = cat(2,vecTrialStimTypeUnsplit,intStim*ones(1,size(matRespUnsplitStim(:,:,intStim),2)));
	end
	%split trial idx
	matSplitTrialIdx = bsxfun(@plus,((1:numel(vecTrialStimTypeUnsplit))-1)*intBinsPerTrial,(1:intBinsPerTrial)')';
end

