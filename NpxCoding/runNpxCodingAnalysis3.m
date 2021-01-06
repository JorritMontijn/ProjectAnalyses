%% aim
%{
"distinct subnetworks in mouse visual system for visual and non-visual signals"

Predict noise in area 2 from noise in area 1, calculate angle between
prediction and real, and compare with shuffled. If real diff is lower, then
we found inter-areal structure of high-d noise. Repeat after removing
independent noise; does difference hold?   

%}
%% aim 2
%{
"non-visual activity in visual system is globally invariant to stimulus orientation"

1) Predict noise in area y from noise in area x for stimulus z using stimulus
z; compare with using stimulus z+1,+2... 
2) compare this baseline to shuffled where generalizability is destroyed (how?)
3) compare also to predictability where z+1 is rotated to z, using
difference in f-primes between z and z+1

Possibly need to reduce to predictive subspace + linear regression as in Semedo et al., fig S6
%}

%% define qualifying areas
vecN = [];
vecT = [];
	
cellUseAreas = {...
	'Subiculum',...
	'Posterior complex of the thalamus',...
	'Anterior pretectal nucleus',...
	'Superior colliculus',...
	'Lateral posterior nucleus of the thalamus',...
	'Nucleus of the optic tract',...
	'Dorsal part of the lateral geniculate complex',...
	'Primary visual area',...
	'posteromedial visual area',...
	'Anteromedial visual area',...
	'Retrosplenial area',...
	'Anterior area',...
	};


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2);

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Rec;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Rec}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellStim(cellfun(@(x) x.structEP.intStimTypes,sThisRec.cellStim) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellStim);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = structStim.Orientation;
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Rec},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	%% select area 1
	for intArea1=1:intAreas
		strArea1 = cellUseAreas{intArea1};
		indArea1Neurons = contains({sUseNeuron.Area},strArea1,'IgnoreCase',true);
		if sum(indArea1Neurons) < 40 || numel(vecStimOnTime) < 500, continue;end
		vecN(end+1) = sum(indArea1Neurons);
		vecT(end+1) = numel(vecStimOnTime);
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
		matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
		
		%% split data
		%prep
		vecUseStimTypes = unique(vecOrientation);
		sParams = struct;
		sParams.intSizeX = 35;
		sParams.intSizeY = 5;
		sParams.vecUseStimTypes = vecUseStimTypes;
		matData = matSpikeCountsArea1;
		vecTrialStimType = vecOrientation;
		%remove nulls
		matNullNeurons = doDimRemoveNulls(matData,vecTrialStimType);
		indRem = any(matNullNeurons,2);
		matData(indRem,:) = [];
		%split
		[cellMatX,	cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sParams);
		%dblLambda = 1;
		%matPredictionsR2 = doDimPredFull(cellMatX,cellMatY,dblLambda);
		%mean(matPredictionsR2(:))
		%matPredDimDepR2 = doDimPredRR(cellMatX,cellMatY)
		
		%% find optimal regularisation parameter for ridge regression
		%folds
		intFoldK = 10;
		
		%line search
		vecLambda = 2.^[0:0.2:20];%[0:(ceil(log2(dblLambda))+2)];
		vecLambdaPerf = nan(size(vecLambda));
		for intL=1:numel(vecLambda)
			dblL = vecLambda(intL);
			[vecR2_Ridge_CV, sRidge] = doRidge_CV(cellMatX{1,1}, cellMatY{1,1}, intFoldK, dblL);
			vecLambdaPerf(intL) = mean(vecR2_Ridge_CV);
		end
		intW = 10;
		vecF = [nan([1 intW]) vecLambdaPerf nan([1 intW])];
		vecLPF = filter((1/intW)*ones(1,intW),1,vecF);
		vecLPF = vecLPF((intW+1):(end-intW));
		[dummy,intP] = max(vecLPF);
		dblLambda = vecLambda(intP);
		
		%plot
		%scatter(vecLambda,vecLambdaPerf);
		%set(gca,'xscale','log')
		%% decode
		matDecData = cell2mat(cellMatX(1,:)');
		vecDecStimType = cell2vec(cellfun(@(x,y) y*ones(1,x),vec2cell(cellfun(@(x) size(x,1),cellMatX(1,:))),vec2cell(1:size(cellMatX,2)),'UniformOutput',false));
		
		[dblPerformance1,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion1] = ...
			doCrossValidatedDecodingLR(matDecData,vecDecStimType,2,dblLambda);

		[dblPerformance2,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion2] = ...
			doCrossValidatedDecodingLR(matData,vecTrialStimType,2,dblLambda);

		
figure
subplot(1,2,1)
imagesc(matConfusion1)
title(sprintf('%s, hit=%.3f',strArea1,dblPerformance1))

subplot(1,2,2)
imagesc(matConfusion2)
title(sprintf('%s, hit=%.3f',strArea1,dblPerformance2))

continue
		%% run self-prediction
		intFoldK = 10;
		%dblLambda = 1/3000;
		vecRank = 1:sParams.intSizeY;
		tic
		warning('off','MATLAB:sqrtm:SingularMatrix');
		matPredictionsR2 = doDimPredRRR_CV(cellMatX,cellMatY,vecRank,intFoldK,dblLambda);
		warning('on','MATLAB:sqrtm:SingularMatrix');
		toc
		
		figure
		%plot(squeeze(mean(matPredDimDepR2,1))')
		plot(squeeze(mean(matPredictionsR2,1))')
		
		%vecR2 = getRegInSpace(matX,matSubspace,matY,dblLambda)
		return
		%% select area 2
		for intArea2=(intArea1):intAreas
			strArea2 = cellUseAreas{intArea2}
			indArea2Neurons = contains({sUseNeuron.Area},strArea2,'IgnoreCase',true);
			if sum(indArea2Neurons) == 0, continue;end
			
			%% get orientation responses & single-trial population noise
			sArea2Neurons = sUseNeuron(indArea2Neurons);
			
			%% get spike times
			cellSpikeTimes2 = {sArea2Neurons.SpikeTimes};
			if intArea1 == intArea2
				%real
				matSpikeCountsArea2 = matSpikeCountsArea1;
				% shuffle
				matSpikeCountsShuffled2 = matSpikeCountsShuffled1;
			else
				%real
				matSpikeCountsArea2 = getSpikeCounts(cellSpikeTimes2,vecStimOnTime,vecStimOffTime);
				% shuffle
				matSpikeCountsShuffled2 = matSpikeCountsArea2(:,randperm(size(matSpikeCountsArea2,2)));
			end
			
		end
	end
end
