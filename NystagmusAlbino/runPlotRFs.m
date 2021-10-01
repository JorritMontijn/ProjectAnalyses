%% load data
strDataPath = 'F:\Data\Processed\Neuropixels';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		if ~isfield(sLoad.sAP,'sPupil') || isempty(sLoad.sAP.sPupil),continue;end
		sLoad.sAP.name = sFiles(intFile).name;
		if isempty(sExp)
			sExp = sLoad.sAP;
		else
			sExp(end+1) = sLoad.sAP;
		end
	end
end
%MP_20200115 eye tracking remove last stimulus (gunk in eye)
cellUseForEyeTrackingMP = {'20191120','20191121','20191122','20191210','20191211','20191212','20191213','20191216','20191217','20200116','20200116R02'}; %don't forget to set high vid lum as blinks
cellUseForEyeTrackingMA = {'20210212','20210215','20210218','20210220','20210225','20210301'};
cellUseForEyeTracking = cat(2,cellUseForEyeTrackingMA,cellUseForEyeTrackingMP);
strTargetPath = 'F:\Data\Results\AlbinoProject';

%% plot
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
indRemoveNoEye = ~contains(strrep(cellExperiment,'-',''),cellUseForEyeTracking);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
for intSubType=1:2
	if intSubType == 1
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType) & ~indRemoveNoEye
	vecRunRecs = find(indUseRecs);
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false)
		vecBlocksRF = find(contains(cellStimType,'receptivefield','IgnoreCase',true))
		for intBlockIdx=1:numel(vecBlocksRF)
			intBlock = vecBlocksRF(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			
			if isfield(sBlock,'vecPupilStimOn')
				vecPupilStimOn = sBlock.vecPupilStimOn;
				vecPupilStimOff = sBlock.vecPupilStimOff;
			else
				vecPupilStimOn = sBlock.vecStimOnTime;
				vecPupilStimOff = sBlock.vecStimOffTime;
			end
			
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
			
			%% go through patches and concatenate repetitions
			cellCellOnSpikeT = cellfill(cell(size(sStimObject(1).LinLoc)),[1 intNeurons]); %{N}{X,Y}{Rep} = [T_s1 T_s2 ... T_sN]
			cellCellOffSpikeT = cellfill(cell(size(sStimObject(1).LinLoc)),[1 intNeurons]); %{N}{X,Y}{Rep} = [T_s1 T_s2 ... T_sN]
			vecStimOn = sBlock.vecStimOnTime;
			vecStimOff = sBlock.vecStimOffTime;
			intNeurons = numel(sRec.sCluster);
			matLinLoc = sStimObject(1).LinLoc;
			vecPatches = unique(matLinLoc(:));
			dblDur = median(vecStimOff - vecStimOn);
			matMeanCountsOn = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matSdCountsOn = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matMeanCountsOff = nan([size(sStimObject(1).LinLoc) intNeurons]);
			matSdCountsOff = nan([size(sStimObject(1).LinLoc) intNeurons]);
			for intNeuron=1:intNeurons
				intNeuron
				vecSpikeT = sRec.sCluster(intNeuron).SpikeTimes;
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeT,vecStimOn,dblDur);
				for intPatchIdx=1:numel(vecPatches)
					intPatch = vecPatches(intPatchIdx);
					indAssign = matLinLoc==intPatch;
					
					%on
					vecPatchIsOn = find(cellfun(@ismember,cellfill(intPatch,size({sStimObject.LinLocOn})),{sStimObject.LinLocOn}));
					for intOnRepIdx=1:numel(vecPatchIsOn)
						intTrial=vecPatchIsOn(intOnRepIdx);
						indTrialSpikes = vecTrialPerSpike==intTrial;
						cellCellOnSpikeT{intNeuron}{indAssign}{intOnRepIdx} = vecTimePerSpike(indTrialSpikes);
					end
					
					%off
					vecPatchIsOff = find(cellfun(@ismember,cellfill(intPatch,size({sStimObject.LinLocOff})),{sStimObject.LinLocOff}));
					for intOffRepIdx=1:numel(vecPatchIsOff)
						intTrial=vecPatchIsOff(intOffRepIdx);
						indTrialSpikes = vecTrialPerSpike==intTrial;
						cellCellOffSpikeT{intNeuron}{indAssign}{intOffRepIdx} = vecTimePerSpike(indTrialSpikes);
					end
					
				end
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
			figure,imagesc(matR(:,:,1))
			%%
			for intNeuron=1:intNeurons
				subplot(2,3,1)
				imagesc(matMeanCountsOn(:,:,intNeuron))
				title(sprintf('%d',intNeuron))
				
				subplot(2,3,2)
				imagesc(matMeanCountsOff(:,:,intNeuron))
				
				subplot(2,3,3)
				imagesc(matMeanCountsOn(:,:,intNeuron) - matMeanCountsOff(:,:,intNeuron))
				
				drawnow;pause
				
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
