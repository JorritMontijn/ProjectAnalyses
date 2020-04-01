%load data
vecPEs = zeros(1,8);
vecRPEs = zeros(1,8);
vecNPEs = zeros(1,8);
vecTotDur = zeros(1,8);
vecTotDurPE = zeros(1,8);
vecTotDurRPE = zeros(1,8);
vecTotDurNPE = zeros(1,8);
vecFracOnePE = nan(1,8);
vecFracOnePE_Lick = nan(1,8);
vecNumNeurons = nan(1,8);
vecNumNeurons_Frac = nan(1,8);
vecDuration = nan(1,8);
vecFracNeuronsPE = nan(1,8);

for intMouse=1:8
	close all
	drawnow;
	clearvars -except vecFrac intMouse vecPEs vecRPEs vecNPEs vecTotDur vecTotDurPE vecTotDurRPE vecTotDurNPE vecFracOnePE vecFracOnePE_Lick vecNumNeurons vecNumNeurons_Frac vecDuration vecFracNeuronsPE
	
	%use neuropil subtraction?
	boolUseNeuropilSubtraction = false;
	boolOnlyTuned = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == 4
			strSes = '20140507';
		elseif intMouse == 5
			strSes = '20140530';
		elseif intMouse == 6
			strSes = '20140604';
		elseif intMouse == 7
			strSes = '20140711';
		elseif intMouse == 8
			strSes = '20140715';
		end
		if boolUseNeuropilSubtraction
			strSes = ['NPS' strSes];
		end
		if boolOnlyTuned
			strSes = [strSes 'OT'];
		end
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['D:\Data\Results\spikeAnalysis\dataPreProAssemblies' strSes '.mat']);
	end
	
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	%#ok<*AGROW>
	strSes = ['AA' strSes];
	
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.1;
		indTooFast = (cellMultiSes{intPopulation}.structStim.vecTrialRespSecs-cellMultiSes{intPopulation}.structStim.SecsOn)<dblRemSecs & cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		%structStim.vecTrialResponse(indTooFast) = 0;
		cellMultiSes{intPopulation}.structStim = remel(cellMultiSes{intPopulation}.structStim,~indTooFast);
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%get stimulus variables
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Contrast'});
		cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		intContrasts = length(cellSelectC);
		
		sTypesO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation'});
		cellSelectO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesO);
		intOrientations = length(cellSelectO);
		
		sTypesCO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation','Contrast'});
		cellSelectCO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesCO);
		
		% get stimulus data
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		intOris = length(vecOrientations);
		[vecStimOris,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Contrast,unique(cellMultiSes{intPopulation}.structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn<3;
		
		
		%get behavioral variables
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indHits = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs -  cellMultiSes{intPopulation}.structStim.SecsOn;
		vecRTs(indMiss) = nan;
		indFast = false(size(indMiss));
		indSlow = false(size(indMiss));
		for intContrast=unique(vecStimContrasts)
			indTheseContrastTrials = vecStimContrasts==intContrast;
			indFast = (vecRTs <= nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indFast;
			indSlow = (vecRTs > nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indSlow;
		end
		
		%get assembly variables
		matReorderedAssemblyCorrs=sAssemblies{intPopulation}.matAssemblyCorrelations(sAssemblies{intPopulation}.vecReorder,sAssemblies{intPopulation}.vecReorder);
		matM = getBlockMeans(matReorderedAssemblyCorrs,sAssemblies{intPopulation}.vecAssemblyIdentity);
		
		%assign random durations to miss trials
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		intTrialDur = round(2 * cellMultiSes{intPopulation}.samplingFreq);
		%intMissDur = round(3 * cellMultiSes{intPopulation}.samplingFreq);
		%sStimTemp.FrameOff(indMiss) = sStimTemp.FrameOn(indMiss) + intMissDur;
		
		%get response matrix
		[matTrialResponse,cellSelectContrasts,vecTrialDur] = getTrialResponseData(cellMultiSes{intPopulation},sStimTemp);
		
		%detect spikes
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		matSpikeCountsHD = sAssemblies{intPopulation}.matSpikeCountsHD;
		matSpikeCounts = getNeuronResponseRespMat(matSpikeCountsHD,sStimTemp);
		matSpikeCounts = round(matSpikeCounts*dblSamplingFreq);
		
		%% create miss/slow/fast selection vectors
		intFrames = size(matSpikeCountsHD,2);
		%miss
		vecMissOn = sStimTemp.FrameOn(indMiss);
		vecMissOff = vecMissOn + intTrialDur;
		%vecMissOff = sStimTemp.FrameOff(indMiss);
		indMissHD = false(1,intFrames);
		for intTrial=1:length(vecMissOn)
			indMissHD(vecMissOn(intTrial):vecMissOff(intTrial)) = true;
		end
		vecMissHD = find(indMissHD);
		%hits
		vecHitsOn = sStimTemp.FrameOn(indHits);
		vecHitsOff = vecHitsOn + intTrialDur;
		%vecHitsOff = sStimTemp.FrameOff(indHits);
		indHitsHD = false(1,intFrames);
		for intTrial=1:length(vecHitsOn)
			indHitsHD(vecHitsOn(intTrial):vecHitsOff(intTrial)) = true;
		end
		vecHitsHD = find(indHitsHD);
		%base
		vecBaseOff = sStimTemp.FrameOn;
		vecBaseOn = sStimTemp.FrameOn - intTrialDur;
		indBaseHD = false(1,intFrames);
		for intTrial=1:length(vecBaseOn)
			indBaseHD(vecBaseOn(intTrial):vecBaseOff(intTrial)) = true;
		end
		vecBaseHD = find(indBaseHD);
		
		%% get members of assemblies
		intNeurons = size(sAssemblies{intPopulation}.matAssemblies,1);
		matAssemblyMembers = sAssemblies{intPopulation}.matAssemblies > 0;
		intAssemblies = max(sAssemblies{intPopulation}.vecAssemblies);
		matAssemblyCoreMembers = false(intAssemblies,intNeurons);
		for intAssembly=1:intAssemblies
			%get random memberships
			intRandIters = 1000;
			indAssembly = sAssemblies{intPopulation}.vecAssemblies == intAssembly;
			matRandMembership = nan(intNeurons,intRandIters);
			for intIter=1:intRandIters
				indRandomAssembly = indAssembly(randperm(length(indAssembly)));
				matRandMembership(:,intIter) = sum(matAssemblyMembers(:,indRandomAssembly),2)/sum(indAssembly);
			end
			dblAlpha = 0.01;
			dblCorrectedSDs = qnorm(1-(dblAlpha/2)/intNeurons);
			vecUpperBound = mean(matRandMembership,2)+std(matRandMembership,[],2)*dblCorrectedSDs;
			vecData = sum(matAssemblyMembers(:,indAssembly),2)/sum(indAssembly);
			indHigher = vecData>vecUpperBound;
			
			%save output
			matAssemblyCoreMembers(intAssembly,:) = indHigher;
		end
		%matAssemblyCoreMembers = true(size(matAssemblyCoreMembers));
		%get assembly types
		vecCoreMemberSize = sum(matAssemblyCoreMembers,2);
		vecAssemblyAutoCorr = diag(matM);
		indRecurringClusters = zscore(vecAssemblyAutoCorr)>-2;
		[dummy,intNonRecurringCluster]=min(zscore(vecAssemblyAutoCorr));
		indRecurringClusters(intNonRecurringCluster) = false;
		indMultiNeuronClusters = vecCoreMemberSize>1;
		indRealAssemblies = indRecurringClusters & indMultiNeuronClusters;
		%indRealAssemblies = true(size(indRealAssemblies));
		
		%% get assembly data
		vecAssemblyStarts = sAssemblies{intPopulation}.vecAssemblyStarts;
		vecAssemblyStops = sAssemblies{intPopulation}.vecAssemblyStops;
		vecAssemblies = sAssemblies{intPopulation}.vecAssemblies';
		vecRealAssemblies = ismember(vecAssemblies,find(indRealAssemblies));
		vecNonRecurringAssemblies = vecAssemblies == intNonRecurringCluster;
		vecAssemblyAEs = zeros(size(vecAssemblyStarts));
		dblSecsCCG = (max(vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq);
		vecPlotCCG = (-ceil(dblSamplingFreq*dblSecsCCG):ceil(dblSamplingFreq*dblSecsCCG));
		intMiddle = ceil(dblSamplingFreq)+1;
		
		vecPEs(intMouse) = vecPEs(intMouse) + numel(vecAssemblies);
		vecRPEs(intMouse) = vecRPEs(intMouse) + sum(vecRealAssemblies);
		vecNPEs(intMouse) = vecNPEs(intMouse) + sum(vecNonRecurringAssemblies);
		vecTotDur(intMouse) = vecTotDur(intMouse) + numel(cellMultiSes{intPopulation}.neuron(1).dFoF)/dblSamplingFreq;
		vecTotDurPE(intMouse) = vecTotDurPE(intMouse) + sum((vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq);
		vecTotDurRPE(intMouse) = vecTotDurRPE(intMouse) + sum((vecAssemblyStops(vecRealAssemblies)-vecAssemblyStarts(vecRealAssemblies))/dblSamplingFreq);
		vecTotDurNPE(intMouse) = vecTotDurNPE(intMouse) + sum((vecAssemblyStops(vecNonRecurringAssemblies)-vecAssemblyStarts(vecNonRecurringAssemblies))/dblSamplingFreq);

		
		%likelihood of at least one PE following stim onset and preceding lick
		indStimTrials = sStimTemp.Contrast > 0;
		dblFracAtLeastOne = sum(sum(bsxfun(@gt,vecAssemblyStarts,sStimTemp.FrameOn(indStimTrials)') & bsxfun(@lt,vecAssemblyStarts,(sStimTemp.FrameOn(indStimTrials) + intTrialDur)'),2) > 0) / sum(indStimTrials);
		vecFracOnePE(intMouse) = nanmean([vecFracOnePE(intMouse) dblFracAtLeastOne]);

		dblFracAtLeastOneLick = sum(sum(bsxfun(@gt,vecAssemblyStarts,(sStimTemp.FrameOff(indHits) - intTrialDur)') & bsxfun(@lt,vecAssemblyStarts,sStimTemp.FrameOff(indHits)'),2) > 0) / sum(indHits);
		vecFracOnePE_Lick(intMouse) = nanmean([vecFracOnePE_Lick(intMouse) dblFracAtLeastOneLick]);
		
		
		%overall average number of neurons and duration per PE
		vecNrNeurons = sum(sAssemblies{intPopulation}.matAssemblies,1);
		dblNumNeurons = mean(vecNrNeurons);
		vecNumNeurons(intMouse) = nanmean([vecNumNeurons(intMouse) dblNumNeurons]);
		dblFracNeurons = mean(vecNrNeurons/intNeurons);
		vecNumNeurons_Frac(intMouse) = nanmean([vecNumNeurons_Frac(intMouse) dblFracNeurons]);
		
		dblDur = mean((vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq);
		vecDuration(intMouse) = nanmean([vecDuration(intMouse) dblDur]);
		
		%fraction of PEs with fewer than [X] neurons
		vecFracNeuronsPE(intMouse) = nanmean([vecFracNeuronsPE(intMouse) sum(vecNrNeurons>3)/length(vecNrNeurons)]);
		
		%% get cluster location, OSI, and member tuning data
		matLocationsXY = nan(intNeurons,2);
		for intNeuron = 1:intNeurons
			matLocationsXY(intNeuron,1) = cellMultiSes{intPopulation}.neuron(intNeuron).x;
			matLocationsXY(intNeuron,2) = cellMultiSes{intPopulation}.neuron(intNeuron).y;
		end
		dblPix2Micron = cellMultiSes{intPopulation}.xml.dblActualVoxelSizeX/1000;
		matLocationsXY = matLocationsXY*dblPix2Micron; %transform pixels to microns
		matDiffX = abs(bsxfun(@minus,matLocationsXY(:,1),matLocationsXY(:,1)'));
		matDiffY = abs(bsxfun(@minus,matLocationsXY(:,2),matLocationsXY(:,2)'));
		matDist = sqrt(matDiffX.^2+matDiffY.^2);
		
		%get assembly OSIs
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		matRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp);%[assemblies x trials]
		indStimC = ~cellSelectC{1};%cellSelectC{6} | cellSelectC{5};
		matRespAssStim = matRespAssemblies(:,indStimC);
		cellSelectStimO = cell(size(cellSelectO));
		for intOri=1:length(cellSelectStimO)
			cellSelectStimO{intOri} = cellSelectO{intOri}(indStimC);
		end
		sTuning = calcTuningRespMat(matRespAssStim,cellSelectStimO,vecOrientations);
		vecAssOSI = sTuning.vecOSI;
		
		%get OSIs of neurons
		matSpikeCountsStim = matSpikeCounts(:,indStimC);
		sTuning = calcTuningRespMat(matSpikeCountsStim,cellSelectStimO,vecOrientations);
		vecNeuronPA = 2*ang2rad(sTuning.vecPrefAngle);
		mat_dPA = abs(circ_dist2(vecNeuronPA,vecNeuronPA));
		
		%define recurring clusters
		for intMinCoreMembers = 1:10
			vecRecurringClusters = find(indRecurringClusters' & vecCoreMemberSize'>(intMinCoreMembers-1));
			intRecClusters = numel(vecRecurringClusters);
			
			%compare mean distance of core assembly members to random group of same size
			intIters = 1000;
			vecMeanDist = nan(1,intRecClusters);
			matShuffledDist = nan(intIters,intRecClusters);
			
			vecAssemblyMean_dPA = nan(1,intRecClusters);
			matAssemblyMean_dPA_Shuffled = nan(intIters,intRecClusters);
			
			vecAssembly_OSI = nan(1,intRecClusters);
			matAssembly_OSI_Shuffled = nan(intIters,intRecClusters);
			
			vecPotentialMembers = find(any(matAssemblyCoreMembers,1));
			intPotentialMembers = length(vecPotentialMembers);
			intClusterCounter = 0;
			for intAssembly=vecRecurringClusters
				%get members
				intClusterCounter = intClusterCounter + 1;
				vecMembers = matAssemblyCoreMembers(intAssembly,:);
				intMembers = sum(vecMembers);
				matSelect = tril(true(intMembers),-1);
				
				%get real dist, dPA, and OSI
				matRealDist = matDist(vecMembers,vecMembers);
				vecMeanDist(intClusterCounter) = mean(matRealDist(matSelect));
				
				matReal_dPA = mat_dPA(vecMembers,vecMembers);
				vecAssemblyMean_dPA(intClusterCounter) = mean(matReal_dPA(matSelect));
				
				vecAssResp = sum(matSpikeCountsStim(vecMembers,:),1);
				sTuning = calcTuningRespMat(vecAssResp,cellSelectStimO,vecOrientations);
				vecAssembly_OSI(intClusterCounter) = sTuning.vecOSI;
				
				%get shuffled dPAs
				for intIter=1:intIters
					vecShuffledMembers = vecPotentialMembers(randperm(intPotentialMembers,intMembers));
					%dist
					matShuffledDistTemp = matDist(vecShuffledMembers,vecShuffledMembers);
					matShuffledDist(intIter,intClusterCounter) = mean(matShuffledDistTemp(matSelect));
					%dPA
					matShuffled_dPA = mat_dPA(vecShuffledMembers,vecShuffledMembers);
					matAssemblyMean_dPA_Shuffled(intIter,intClusterCounter) = mean(matShuffled_dPA(matSelect));
					%OSI
					vecAssRespShuf = sum(matSpikeCountsStim(vecShuffledMembers,:),1);
					sTuning = calcTuningRespMat(vecAssRespShuf,cellSelectStimO,vecOrientations);
					matAssembly_OSI_Shuffled(intIter,intClusterCounter) = sTuning.vecOSI;
				end
			end
			
			%transform to z-scores
			cellSaveClusterZ{intPopulation,intMinCoreMembers,1} = (vecMeanDist-mean(matShuffledDist,1))./std(matShuffledDist,[],1);
			cellSaveClusterZ{intPopulation,intMinCoreMembers,2} = (vecAssemblyMean_dPA-mean(matAssemblyMean_dPA_Shuffled,1))./std(matAssemblyMean_dPA_Shuffled,[],1);
			cellSaveClusterZ{intPopulation,intMinCoreMembers,3} = (vecAssembly_OSI-mean(matAssembly_OSI_Shuffled,1))./std(matAssembly_OSI_Shuffled,[],1);
		end
	end
	
	%% save data structures
	if sParams.boolSaveData && exist('cellSaveClusterZ','var')
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnalSix' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end

