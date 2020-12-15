%% aim
%{
"distinct subnetworks in mouse visual system for visual and non-visual signals"

At single trial basis, split which info information limiting noise
(projected unto f') and non-limiting noise (orthogonal directions). Then
look at correlation of those noise types between areas. E.g., lgn diff
corrs correlate with v1, but SC non-diff corrs correlate with v1
%}
%% define qualifying areas
cellUseAreas = {'Anterior area',...
	'Anterior pretectal nucleus',...
	'Anteromedial visual area',...
	'Dorsal part of the lateral geniculate complex',...
	'Lateral posterior nucleus of the thalamus',...
	'Nucleus of the optic tract',...
	'Posterior complex of the thalamus',...
	'Primary visual area',...
	'Retrosplenial area',...
	'Subiculum',...
	'Superior colliculus',...
	'posteromedial visual area',...
	};

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim));
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim));
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim));
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim));

%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end

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
		if sum(indArea1Neurons) == 0, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
		matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
		
		%% select area 2
		for intArea2=(intArea1):intAreas
			strArea2 = cellUseAreas{intArea2}
			indArea2Neurons = contains({sUseNeuron.Area},strArea2,'IgnoreCase',true);
			if sum(indArea2Neurons) == 0, continue;end
			%% get orientation responses & single-trial population noise
			sArea2Neurons = sUseNeuron(indArea2Neurons);
			
			%% get spike times
			cellSpikeTimes2 = {sArea2Neurons.SpikeTimes};
			matSpikeCountsArea2 = getSpikeCounts(cellSpikeTimes2,vecStimOnTime,vecStimOffTime);
			
			%% analysis: calculate
			vecR = getAnalysisNpxMD(matSpikeCountsArea1,matSpikeCountsArea2,vecOrientation);
			
			%add to matrices
			matR_PP_All(intArea1,intArea2,intRec) = vecR(1);
			matR_OP_All(intArea1,intArea2,intRec) = vecR(2);
			matR_OO_All(intArea1,intArea2,intRec) = vecR(3);
			matR_PO_All(intArea1,intArea2,intRec) = vecR(4);
			
			matR_PP_All(intArea2,intArea1,intRec) = vecR(1);
			matR_OP_All(intArea2,intArea1,intRec) = vecR(2);
			matR_OO_All(intArea2,intArea1,intRec) = vecR(3);
			matR_PO_All(intArea2,intArea1,intRec) = vecR(4);
		end
	end
end
toc
%% average & remove empty rows
matR_PP = nanmean(matR_PP_All,3);
matR_OP = nanmean(matR_OP_All,3);
matR_PO = nanmean(matR_PO_All,3);
matR_OO = nanmean(matR_OO_All,3);

vecRemR = all(isnan(matR_PP),1);
vecRemC = all(isnan(matR_PP),2);

matR_PP(vecRemR,:) = [];
matR_PP(:,vecRemC) = [];
matR_OP(vecRemR,:) = [];
matR_OP(:,vecRemC) = [];
matR_OO(vecRemR,:) = [];
matR_OO(:,vecRemC) = [];
matR_PO(vecRemR,:) = [];
matR_PO(:,vecRemC) = [];

vecLimC = max(flat(abs([matR_OP matR_PO])))*[-1 1];
subplot(2,3,1)
imagesc(matR_PP,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_PP,vecLimC,redblue);
title(sprintf('Corr Parallel-parallel var'))

subplot(2,3,2)
imagesc(matR_OP,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_OP,vecLimC,redblue);
title(sprintf('Corr Orthogonal-parallel var'))

subplot(2,3,4)
imagesc(matR_OO,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_OO,vecLimC,redblue);
title(sprintf('Corr Orthogonal-Orthogonal var'))

subplot(2,3,5)
imagesc(matR_PO,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_PO,vecLimC,redblue);
title(sprintf('Corr Parallel-Orthogonal var'))

