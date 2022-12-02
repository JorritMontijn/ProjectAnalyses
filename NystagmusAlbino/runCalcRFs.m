
%% load data and define groups
%strDataPath
%cellUseForEyeTracking
%strTargetPath
%cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%cellUseAreas{2} = {'nucleus of the optic tract'};
%cellUseAreas{3} = {'superior colliculus'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
%cellSubjectGroups = {'BL6','DBA'};
runHeaderNOT;

%% plot
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
indRemoveNoEye = ~contains(strrep(cellExperiment,'-',''),cellUseForEyeTracking);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
for intSubType=1:2
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType) & ~indRemoveNoEye;
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	%vecRunRecs = intBestRec;
	
	
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		vecBlocksNM = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
		%get timing for DG
		intBlock = vecBlocksDG(1);
		sBlock = sRec.cellBlock{intBlock};
		if isfield(sBlock,'vecPupilStimOn')
			vecPupilLatency = sBlock.vecPupilStimOn-sBlock.vecStimOnTime;
		else
			vecPupilLatency = 0;
		end
		vecBlocksRF = find(contains(cellStimType,'receptivefield','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksRF)
			intBlock = vecBlocksRF(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			strName=['RF_' sRec.sJson.subject '_' sRec.sJson.date 'B' num2str(intBlock)];
			
			vecPupilStimOn = sBlock.vecStimOnTime+median(vecPupilLatency);
			vecPupilStimOff = sBlock.vecStimOffTime+median(vecPupilLatency);
			
			%% get pupil data
			dblSampNi = str2double(sRec.sSources.sMeta.niSampRate);
			dblFirstSamp = str2double(sRec.sSources.sMeta.firstSample);
			vecPupilTime = sRec.sPupil.vecTime;
			vecPupilLocX = sRec.sPupil.vecCenterX;
			vecPupilLocY = sRec.sPupil.vecCenterY;
			vecPupilSize = sRec.sPupil.vecRadius;
			vecPupilSync = sRec.sPupil.vecSyncLum;
			if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
				vecBlinks = sRec.sPupil.vecBlinks;
			else
				%filter absvidlum
				dblLowPass = 0.01/(1/median(diff(vecPupilTime)));
				[fb,fa] = butter(2,dblLowPass,'high');
				vecAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum));
				vecBlinks = vecAbsVidLum > 5;
			end
			dblFs = 1/median(diff(vecPupilTime));
			
			%% get rf resps
			%{
			indNOT = contains({sRec.sCluster.Area},'Nucleus of the optic tract','IgnoreCase',true);
			sNeuronsNOT = sRec.sCluster(indNOT);
			
			intNeuron = 1;
			vecSpikes = sNeuronsNOT(intNeuron).SpikeTimes;
			%}
			%Rfs
			sStimObject = sBlock.sStimObject;
			intStims = numel(sStimObject);
			
			%test
			vecBinEdges = [-inf vecPupilStimOn vecPupilStimOn(end)+median(diff(vecPupilStimOn)) inf];
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecPupilTime,double(vecBlinks),vecBinEdges);
			vecBlinkFractionPerTrial = vecMeans(2:(end-1));
			vecBlinkFractionPerTrial(isnan(vecBlinkFractionPerTrial))=0;
			
			%% prep variables
			intNeurons = numel(sRec.sCluster);
			cellCellOnSpikeT = cellfill(cell(size(sStimObject(1).LinLoc)),[1 intNeurons]); %{N}{X,Y}{Rep} = [T_s1 T_s2 ... T_sN]
			cellCellOffSpikeT = cellfill(cell(size(sStimObject(1).LinLoc)),[1 intNeurons]); %{N}{X,Y}{Rep} = [T_s1 T_s2 ... T_sN]
			vecStimOn = sBlock.vecStimOnTime;
			vecStimOff = sBlock.vecStimOffTime;
			matLinLoc = sStimObject(1).LinLoc;
			vecPatches = unique(matLinLoc(:));
			dblDur = median(vecStimOff - vecStimOn);
			matZetaOn = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matZetaOff = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matMeanCountsOn = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matSdCountsOn = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matMeanCountsOff = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matSdCountsOff = nan([size(sStimObject(1).LinLoc) intNeurons]);
			%% go through cells
			for intNeuron=1:intNeurons
				fprintf('Running %s (rec %d/%d), neuron %d/%d\n',strName,intRecIdx,numel(vecRunRecs),intNeuron,intNeurons);
				vecSpikeT = sRec.sCluster(intNeuron).SpikeTimes;
				matTempZetaOn = nan(size(sStimObject(1).LinLoc));
				matTempZetaOff = nan(size(sStimObject(1).LinLoc));
				%cull spikes
				dblMinT = min(vecStimOn) - 10*dblDur;
				dblMaxT = max(vecStimOn) + 10*dblDur;
				vecSpikeT((vecSpikeT < dblMinT) | (vecSpikeT > dblMaxT)) = [];
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeT,vecStimOn,dblDur);
				%go through patches and concatenate repetitions
				for intPatchIdx=1:numel(vecPatches)
					intPatch = vecPatches(intPatchIdx);
					indAssign = matLinLoc==intPatch;
					
					%on
					indOn = cellfun(@ismember,cellfill(intPatch,size({sStimObject.LinLocOn})),{sStimObject.LinLocOn});
					vecPatchIsOn = find(vecBlinkFractionPerTrial(:) < 0.1 & indOn(:));
					for intOnRepIdx=1:numel(vecPatchIsOn)
						intTrial=vecPatchIsOn(intOnRepIdx);
						indTrialSpikes = vecTrialPerSpike==intTrial;
						cellCellOnSpikeT{intNeuron}{indAssign}{intOnRepIdx} = vecTimePerSpike(indTrialSpikes);
					end
					
					%off
					indOff = cellfun(@ismember,cellfill(intPatch,size({sStimObject.LinLocOff})),{sStimObject.LinLocOff});
					vecPatchIsOff = find(vecBlinkFractionPerTrial(:) < 0.1 & indOff(:));
					for intOffRepIdx=1:numel(vecPatchIsOff)
						intTrial=vecPatchIsOff(intOffRepIdx);
						indTrialSpikes = vecTrialPerSpike==intTrial;
						cellCellOffSpikeT{intNeuron}{indAssign}{intOffRepIdx} = vecTimePerSpike(indTrialSpikes);
					end
					
					%% calculate zeta per patch
					%on
					vecDur = vecStimOff - vecStimOn;
					dblBaselineDur = dblDur;
					dblTotDur = dblDur + dblBaselineDur;
					vecUseStimOn = vecStimOn(vecPatchIsOn)-dblBaselineDur/2;
					intResampNum = 100;
					dblZetaP_on = zetatest(vecSpikeT,vecUseStimOn,dblTotDur,intResampNum);
					matTempZetaOn(indAssign) = dblZetaP_on;
					
					%off
					vecUseStimOff = vecStimOn(vecPatchIsOff)-dblBaselineDur/2;
					intResampNum = 100;
					dblZetaP_off = zetatest(vecSpikeT,vecUseStimOff,dblTotDur,intResampNum);
					matTempZetaOff(indAssign) = dblZetaP_off;
					
				end
				matZetaOn(:,:,intNeuron) = matTempZetaOn;
				matZetaOff(:,:,intNeuron) = matTempZetaOff;
				matMeanCountsOn(:,:,intNeuron) = cellfun(@(x) mean(cellfun(@numel,x)),cellCellOnSpikeT{intNeuron});
				matSdCountsOn(:,:,intNeuron) = cellfun(@(x) std(cellfun(@numel,x)),cellCellOnSpikeT{intNeuron});
				matMeanCountsOff(:,:,intNeuron) = cellfun(@(x) mean(cellfun(@numel,x)),cellCellOffSpikeT{intNeuron});
				matSdCountsOff(:,:,intNeuron) = cellfun(@(x) std(cellfun(@numel,x)),cellCellOnSpikeT{intNeuron});
				
			end
			matSdCountsOn(isnan(matSdCountsOn))=inf;
			matSdCountsOff(isnan(matSdCountsOn))=inf;
			
			
			matRepsOn = cellfun(@numel,cellCellOnSpikeT{1});
			matRepsOff = cellfun(@numel,cellCellOffSpikeT{1});
			matOnR = matMeanCountsOn(:,:,intNeuron);
			matOffR = matMeanCountsOff(:,:,intNeuron);
			matR = matOnR;
			
			%% save data
			strSavePath = [strTargetPath filesep 'RF_data'];
			save(fullpath(strSavePath,strName),'sBlock',...
				'intBlock',...
				'matZetaOn',...
				'matZetaOff',...
				'matMeanCountsOn',...
				'matSdCountsOn',...
				'matMeanCountsOff',...
				'matSdCountsOff',...
				'matRepsOn',...
				'matRepsOff');
			
			if 0
				%% plot
				figure,imagesc(matR(:,:,1))
				for intNeuron=1:intNeurons
					subplot(2,4,1)
					imagesc(matMeanCountsOn(:,:,intNeuron))
					title(sprintf('Mean counts on, %d',intNeuron))
					
					subplot(2,4,2)
					imagesc(matMeanCountsOff(:,:,intNeuron))
					title(sprintf('Mean counts off, %d',intNeuron))
					
					subplot(2,4,3)
					imagesc(matMeanCountsOn(:,:,intNeuron) - matMeanCountsOff(:,:,intNeuron))
					title(sprintf('Mean counts on-off, %d',intNeuron))
					
					subplot(2,4,4)
					imagesc(matMeanCountsOn(:,:,intNeuron) + matMeanCountsOff(:,:,intNeuron))
					title(sprintf('Mean counts on+off, %d',intNeuron))
					colorbar;
					
					
					subplot(2,4,5)
					matZetaPlotOn = -norminv( matZetaOn(:,:,intNeuron) /2);
					matZetaPlotOn(matZetaPlotOn<2)=2;
					imagesc(matZetaPlotOn);
					title(sprintf('ZETA on, %d',intNeuron))
					
					subplot(2,4,6)
					matZetaPlotOff = -norminv( matZetaOff(:,:,intNeuron) /2);
					matZetaPlotOff(matZetaPlotOff<2)=2;
					imagesc(matZetaPlotOff)
					title(sprintf('ZETA off, %d',intNeuron))
					
					subplot(2,4,7)
					imagesc(matZetaPlotOn - matZetaPlotOff)
					title(sprintf('ZETA on-off, %d',intNeuron))
					
					subplot(2,4,8)
					imagesc(matZetaPlotOn + matZetaPlotOff)
					title(sprintf('ZETA on+off, %d',intNeuron))
					colorbar;
					
					drawnow;pause
					
					if 0
						%%
						save([strTargetPath filesep 'RF_ZETA']);
						export_fig([strTargetPath filesep 'ExampleRF_' sRec.sJson.subject '_' sRec.sJson.date '_N316.tif']);
						export_fig([strTargetPath filesep 'ExampleRF_' sRec.sJson.subject '_' sRec.sJson.date '_N316.pdf']);
					end
				end
			end
			%% go through stims
			%{
			matRespOn = repmat(nan(size(sStimObject(1).LinLoc)),[1 1 intStims]);
			matRespOff = repmat(nan(size(sStimObject(1).LinLoc)),[1 1 intStims]); %#ok<*REPMAT>
			
			for intStim=1:intStims
				%get spikes this stim
				dblRate = [];
				
				%add to on matrix & off matrix
				matThisOn = nan(size(sStimObject(1).LinLoc));
				matThisOff = nan(size(sStimObject(1).LinLoc));
				matThisOn(sStimObject.LinLocOn) = dblRate;
				matThisOff(sStimObject.LinLocOff) = dblRate;
				matRespOn(:,:,intStim) = matThisOn;
				matRespOff(:,:,intStim) = matThisOff;
			end
			matRespOnMu = mean(matRespOn,3);
			matRespOnSd = std(matRespOn,[],3);
			matRespOffMu = mean(matRespOff,3);
			matRespOffSd = std(matRespOff,[],3);
			%}
			
		end
	end
	
end
