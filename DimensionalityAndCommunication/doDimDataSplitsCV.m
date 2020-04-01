function [cellMatX1,cellMatX2,cellNeuronsX,cellMatY1,cellMatY2,cellNeuronsY,cellTrials,cellTrials1,cellTrials2] = ...
		doDimDataSplitsCV(matData,vecTrialStimType,sParams)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%% check input
	if ~exist('sParams','var'),sParams=struct;end
	if isfield(sParams,'intSizeX'),intSizeX=sParams.intSizeX;else,intSizeX = 110;end
	if isfield(sParams,'intSizeY'),intSizeY=sParams.intSizeY;else,intSizeY = 30;end
	if isfield(sParams,'intResamplings'),intResamplings=sParams.intResamplings;else,intResamplings = 10;end
	if isfield(sParams,'vecCellArea'),vecCellArea=sParams.vecCellArea;else,vecCellArea = ones(1,size(matData,1));end
	if isfield(sParams,'intWithinArea'),intWithinArea=sParams.intWithinArea;else,intWithinArea = 1;end %0, across; 1, within 1; 2, within 2
	if isfield(sParams,'intMaxReps'),intMaxReps=sParams.intMaxReps;else,intMaxReps = inf;end
	if isfield(sParams,'vecUseStimTypes'),vecUseStimTypes=sParams.vecUseStimTypes;
	else
		[x,y,vecUseStimTypes] = find(sort(unique(vecTrialStimType),'ascend'),2,'first');
	end
	if isfield(sParams,'boolVerbose'),boolVerbose=sParams.boolVerbose;else,boolVerbose = true;end
	
	
	%% general data (not class-specific)
	intClasses = numel(vecUseStimTypes);
	intNeurons = size(matData,1);
	if intWithinArea == 3 || intWithinArea == 0
		vecNeuronsSource = find(vecCellArea==1);
		vecNeuronsTarget = find(vecCellArea==2);
		if numel(vecNeuronsTarget)<intSizeY,error([mfilename ':NoArea2'],'Insufficient neurons in target area');end
	else
		vecNeuronsSource = find(vecCellArea==intWithinArea);
		vecNeuronsTarget = [];
	end
	%% pre-allocate
	cellMatX1 = cell(intResamplings,intClasses);
	cellMatX2 = cell(intResamplings,intClasses);
	cellNeuronsX = cell(intResamplings,intClasses);
	cellMatY1 = cell(intResamplings,intClasses);
	cellMatY2 = cell(intResamplings,intClasses);
	cellNeuronsY = cell(intResamplings,intClasses);
	cellTrials = cell(intResamplings,intClasses);
	cellTrials1 = cell(intResamplings,intClasses);
	cellTrials2 = cell(intResamplings,intClasses);
	intMaxReps = min(min(sum(bsxfun(@eq,vecTrialStimType,vecUseStimTypes'),2)),intMaxReps);
	
	%% define data for random selection
	for intResampling=1:intResamplings
		%% select data
		if numel(vecNeuronsSource) == intSizeX
			vecNeuronsX = vecNeuronsSource;
		else
			vecNeuronsX = vecNeuronsSource(randperm(numel(vecNeuronsSource),intSizeX));
		end
		if isempty(vecNeuronsTarget)
			vecOthers = vecNeuronsSource(~ismember(vecNeuronsSource,vecNeuronsX));
			if length(vecOthers) == intSizeY
				vecNeuronsY = vecOthers;
			else
				vecNeuronsY = vecOthers(randperm(length(vecOthers),intSizeY));
			end
		else
			if numel(vecNeuronsTarget) == intSizeY
				vecNeuronsY = vecNeuronsTarget;
			else
				vecNeuronsY = vecNeuronsTarget(randperm(numel(vecNeuronsTarget),intSizeY));
			end
		end
		for intClass=1:intClasses
			%select data
			vecTrials = find(vecTrialStimType==vecUseStimTypes(intClass));
			intSamples = numel(vecTrials);
			intUseSamples = floor(intSamples/2);
			vecRandAssign = randperm(intSamples);
			vecTrials1 = vecRandAssign(1:intUseSamples);
			vecTrials2 = vecRandAssign((intUseSamples+1):(intUseSamples*2));
			if intUseSamples > intMaxReps
				vecRandSelect = randperm(intUseSamples,intMaxReps);
				vecTrials1 = vecTrials1(vecRandSelect);
				vecTrials2 = vecTrials2(vecRandSelect);
			end
			intSamples = numel(vecTrials1);
			
			%get data
			if ndims(matData) == 3
				intBins = size(matData,2);
				matDataZ = reshape(matData(:,:,vecTrials),[intNeurons intSamples*intBins]);
			else
				matDataZ = matData(:,vecTrials);
			end
			
			%get data %% => NOT Z-SCORED ANYMORE! [24/10/2017]
			matX1 = matDataZ(vecNeuronsX,vecTrials1)'; %source population 1 (predictor)
			matX2 = matDataZ(vecNeuronsX,vecTrials2)'; %source population 2 (predictor)
			matY1 = matDataZ(vecNeuronsY,vecTrials1)'; %target population 1 to be predicted
			matY2 = matDataZ(vecNeuronsY,vecTrials2)'; %target population 2 to be predicted
			
			%save data
			cellMatX1{intResampling,intClass} = matX1;
			cellMatX2{intResampling,intClass} = matX2;
			cellNeuronsX{intResampling,intClass} = vecNeuronsX;
			cellMatY1{intResampling,intClass} = matY1;
			cellMatY2{intResampling,intClass} = matY2;
			cellNeuronsY{intResampling,intClass} = vecNeuronsY;
			cellTrials{intResampling,intClass} = vecTrials;
			cellTrials1{intResampling,intClass} = vecTrials1;
			cellTrials2{intResampling,intClass} = vecTrials2;
		end
	end
