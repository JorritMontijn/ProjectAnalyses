%% load data
strDataPath = 'F:\Data\Processed\Neuropixels';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		if ~isfield(sLoad.sAP,'sPupil') || isempty(sLoad.sAP.sPupil),continue;end
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

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% pre-allocate
vecAggPerfWt = nan(0,numel(cellUseAreas),1);
matAggPerfWt = nan(0,numel(cellUseAreas),intStimNr);
vecAggPerfAlb = nan(0,numel(cellUseAreas),1);
matAggPerfAlb = nan(0,numel(cellUseAreas),intStimNr);
			
%% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1:2
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = -15;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = 0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx)
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
			
			% get blinking
			vecBinEdges = [-inf vecPupilStimOn vecPupilStimOn(end)+median(diff(vecPupilStimOn)) inf];
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecPupilTime,double(vecBlinks),vecBinEdges);
			vecBlinkFractionPerTrial = vecMeans(2:(end-1));
			
			%remove trials with blinking
			indRemTrials = vecBlinkFractionPerTrial > 0.1;
			
			%% prep data
			%split by ori
			sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			vecOrientation(indRemTrials) = [];
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24,continue,end
			
			%get data matrix
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			dblDur = median(vecStimOffTime-vecStimOnTime);
			matData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblDur);
			
			%include?
			vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
			indUseCells = arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
			matUseData = matData(indUseCells,:);
			
			%% split cells into areas
			%cortex
			cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
			%NOT
			cellUseAreas{2} = {'nucleus of the optic tract'};
			%hippocampus
			cellUseAreas{3} = {'Hippocampal formation','CA1','CA2','CA3','subiculum','dentate gyrus'};
			
			%build cell vectors
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			end
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			%% go through adjacent stimuli
			%pre-allocate
			vecPerf = nan(numel(cellUseAreas),1);
			matPerf = nan(numel(cellUseAreas),intStimNr);
			
			%params
			intTypeCV = 2; %leave repetition out
			vecOriNoDir = mod(vecOrientation,180);
			vecTrialTypesNoDir = deg2rad(vecOriNoDir)*2;
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTrialTypesNoDir);
			dblLambda = 100;
			
			%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
			[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			intStimNr = numel(vecUnique);
			for intArea = 1:numel(cellUseAreas)
				%select cells
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				% decode all
				matUseData = matData(vecSelectCells,:);
				%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecOriIdx,intTypeCV,vecPriorDistribution);
				%vecPerfTM(intArea) = dblPerformanceCV;
				[dblPerformanceLR,vecDecodedIndexCV_LR,matPosteriorProbability,matWeights,dblMeanErrorDegsLR,matConfusionLR] = doCrossValidatedDecodingLR(matUseData,vecOrientation,intTypeCV,dblLambda);
				vecPerf(intArea) = dblPerformanceLR;
				
				%adjacent stims
				for intStim1=1:intStimNr
					intStim2 = intStim1-1;
					if intStim2==0
						intStim2=intStimNr;
					end
					indStim1 = cellSelect{intStim1};
					indStim2 = cellSelect{intStim2};
					vecUseTrialTypes = vecTrialTypes(indStim1 | indStim2);
					
					%select data
					matUseData = matData(vecSelectCells,indStim1 | indStim2);
					%decode
					%vecUsePriorDistribution = vecCounts([intStim1 intStim2]);
					%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecUseTrialTypes,intTypeCV,vecUsePriorDistribution);
					%matPerfTM(intArea,intStim1) = dblPerformanceCV;
					[dblPerformanceLR,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseData,vecUseTrialTypes,intTypeCV,dblLambda);
					matPerf(intArea,intStim1) = dblPerformanceLR;
					
				end
			end
			if strcmp(strSubjectType,'BL6')
				vecAggPerfWt(end+1,:,:) = vecPerf;
				matAggPerfWt(end+1,:,:) = matPerf;
			else
				vecAggPerfAlb(end+1,:,:) = vecPerf;
				matAggPerfAlb(end+1,:,:) = matPerf;
			end
		end
	end
	
	%% plot
	
end

vecMeanWtPerfCtx = squeeze(nanmean(matAggPerfWt(:,1,:),1))
vecMeanWtPerfNot = squeeze(nanmean(matAggPerfWt(:,2,:),1))
vecMeanWtPerfHip = squeeze(nanmean(matAggPerfWt(:,3,:),1))

subplot(2,3,1)
polar(deg2rad(vecUnique),vecMeanWtPerfCtx-0.5)


subplot(2,3,2)
polar(deg2rad(vecUnique),vecMeanWtPerfNot-0.5)

subplot(2,3,3)
polar(deg2rad(vecUnique),vecMeanWtPerfHip-0.5)

vecMeanAlbPerfCtx = squeeze(nanmean(matAggPerfAlb(:,1,:),1))
vecMeanAlbPerfNot = squeeze(nanmean(matAggPerfAlb(:,2,:),1))
vecMeanAlbPerfHip = squeeze(nanmean(matAggPerfAlb(:,3,:),1))

subplot(2,3,4)
polar(deg2rad(vecUnique),vecMeanAlbPerfCtx-0.5)


subplot(2,3,5)
polar(deg2rad(vecUnique),vecMeanAlbPerfNot-0.5)

subplot(2,3,6)
polar(deg2rad(vecUnique),vecMeanAlbPerfHip-0.5)

%% save
return
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.jpg']));
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.pdf']));