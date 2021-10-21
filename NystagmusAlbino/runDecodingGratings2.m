%%
%{
edit to run like runDecodingGratings, but compare BL6/DBA properly

%}
%% load data
strDataPath = 'E:\DataPreProcessed';
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
strTargetPath = 'D:\Data\Results\AlbinoProject';

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% define area categories
%cortex
cellUseAreas = [];
cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
%cellUseAreas{3} = {'Hippocampal formation','Field CA1','Field CA2','Field CA3','subiculum','dentate gyrus'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
cellAreaGroups = {'Vis. ctx','NOT'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
cellAreaGroupsAbbr = {'Ctx','NOT'};
cellSubjectGroups = {'BL6','DBA'};

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);
	
%% pre-allocate
intMinCells = 10;%NOT+interaction significant at (1-off diag): 7, 8, 9, 10; not sign. at 6; only diag: 6,7,8 (only interaction), 9,10 (both)
intStimNr = 24;
matDiag = diag(diag(true(intStimNr,intStimNr)));
matIsCorrect = circshift(matDiag,-1) | matDiag | circshift(matDiag,1);
dblChanceP = sum(matIsCorrect(:))/numel(matIsCorrect);
vecAggPerfWt = nan(0,numel(cellUseAreas),1);
matAggPerfWt = nan(0,numel(cellUseAreas),intStimNr);
vecAggPerfAlb = nan(0,numel(cellUseAreas),1);
matAggPerfAlb = nan(0,numel(cellUseAreas),intStimNr);
matAggConfusionWt = nan(0,numel(cellUseAreas),intStimNr,intStimNr);
matAggConfusionAlb = nan(0,numel(cellUseAreas),intStimNr,intStimNr);

cellAggPrefWt = cell(0,numel(cellUseAreas));
cellAggPrefAlb = cell(0,numel(cellUseAreas));

cellAggMeanRespWt = cellfill(nan(0,intStimNr),[1 3]);
cellAggMeanRespAlb = cellfill(nan(0,intStimNr),[1 3]);
		
%% run
vecUseCellNum = [];
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
cellAllSubDataAlb = cell(size(cellAreaGroupsAbbr));
cellAllSubDataBL6 = cell(size(cellAreaGroupsAbbr));
cellAllOriVecsAlb = cell(size(cellAreaGroupsAbbr));
cellAllOriVecsBL6 = cell(size(cellAreaGroupsAbbr));
for intSubType=1:2
	intPopCounter = 0;
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
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
			intPopCounter = intPopCounter + 1;
			
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
			%indUseCells = arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%% split cells into areas
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
			cellPref = cell(numel(cellUseAreas),1);
			
			%params
			intTypeCV = 2; %leave repetition out
			vecOriNoDir = mod(vecOrientation,180);
			vecTrialTypesNoDir = deg2rad(vecOriNoDir)*2;
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTrialTypesNoDir);
			dblLambda = 100;
			
			%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
			[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			vecPriorDistribution = vecCounts;
			intStimNr = numel(vecUnique);
			
			for intArea = 1:numel(cellUseAreas)
				%% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				%% calc tuning curves
				matUseData = matData(vecSelectCells,:);
				sOut = getTuningCurves(matUseData,vecOrientation);
				vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
				[h,d,vecP] = fdr_bh(sOut.vecFitP);
				indKeepCells = vecP<0.05;
				if intArea==3
					indKeepCells = true(size(indKeepCells));
				end
				vecPrefOri(~indKeepCells)=[];
				cellPref{intArea} = vecPrefOri;
				if strcmp(strSubjectType,'BL6')
					cellAggMeanRespWt{intArea} = cat(1,cellAggMeanRespWt{intArea},sOut.matMeanResp(indKeepCells,:));
				else
					cellAggMeanRespAlb{intArea} = cat(1,cellAggMeanRespAlb{intArea},sOut.matMeanResp(indKeepCells,:));
				end
				vecUseCellNum(end+1) = sum(indKeepCells);
				
				if strcmp(strSubjectType,'BL6')
					cellAllSubDataBL6{intArea}(end+1) = {matUseData(indKeepCells,:)};
					cellAllOriVecsBL6{intArea}(end+1) = {vecOrientation};
				else
					cellAllSubDataAlb{intArea}(end+1) = {matUseData(indKeepCells,:)};
					cellAllOriVecsAlb{intArea}(end+1) = {vecOrientation};
				end
				
				%% decode all
				if sum(indKeepCells) < intMinCells
					continue;
				end
				vecSubSelect = find(indKeepCells);
				%run multiple permutations and average
				intRunPerms = 50;
				intNumN = numel(vecSubSelect);
				intChoosePerms = min([intRunPerms nchoosek(intNumN,intMinCells)-1]);
				vecPermPerf = nan(1,intChoosePerms);
				boolChosen = false;
				intChooseIdx=1;
				matPerms = nan(intChoosePerms,intMinCells);
				while ~boolChosen
					vecRandPerm = sort(randperm(intNumN,intMinCells));
					if ~any(all(bsxfun(@eq,matPerms,vecRandPerm),2))
						matPerms(intChooseIdx,:) = vecRandPerm;
						intChooseIdx = intChooseIdx + 1;
					end
					if intChooseIdx > intChoosePerms
						boolChosen = true;
					end
				end
				
				for intPerm=1:intChoosePerms
					vecUseRandPerm = matPerms(intPerm,:);
					vecUseSubSelect = vecSubSelect(vecUseRandPerm);
					matUseSubData = matUseData(vecUseSubSelect,:);
					%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecOriIdx,intTypeCV,vecPriorDistribution);
					%vecPerfTM(intArea) = dblPerformanceCV;
					[dblPerformanceLR,vecDecodedIndexCV_LR,matPosteriorProbability,dblMeanErrorDegsLR,matConfusionLR,matWeights] = ...
						doCrossValidatedDecodingLR(matUseSubData,vecOrientation,intTypeCV,vecPriorDistribution,dblLambda);
					vecPermPerf(intPerm) = sum(matConfusionLR(matIsCorrect))/sum(matConfusionLR(:));
				end
				vecPerf(intArea) = mean(vecPermPerf);
				
				if strcmp(strSubjectType,'BL6')
					matAggConfusionWt(intPopCounter,intArea,:,:) = matConfusionLR;
				else
					matAggConfusionAlb(intPopCounter,intArea,:,:) = matConfusionLR;
				end
			end
			
			if strcmp(strSubjectType,'BL6')
				vecAggPerfWt(end+1,:,:) = vecPerf;
				%matAggPerfWt(end+1,:,:) = matPerf;
				%cellAggPrefWt(end+1,:) = cellPref;
			else
				vecAggPerfAlb(end+1,:,:) = vecPerf;
				%matAggPerfAlb(end+1,:,:) = matPerf;
				%cellAggPrefAlb(end+1,:) = cellPref;
			end
			
		end
	end
	
	
	
end

%%
vecPerfCtxBL6 = vecAggPerfWt(:,1);
vecPerfCtxBL6(isnan(vecPerfCtxBL6)) = [];
vecPerfNOTBL6 = vecAggPerfWt(:,2);
vecPerfNOTBL6(isnan(vecPerfNOTBL6)) = [];

vecPerfCtxAlb = vecAggPerfAlb(:,1);
vecPerfCtxAlb(isnan(vecPerfCtxAlb)) = [];
vecPerfNOTAlb = vecAggPerfAlb(:,2);
vecPerfNOTAlb(isnan(vecPerfNOTAlb)) = [];

%% test
[h,p_Ctx] = ttest2(vecPerfCtxBL6,vecPerfCtxAlb);
[h,p_NOT] = ttest2(vecPerfNOTBL6,vecPerfNOTAlb);

%interaction
vecY = cat(1,vecPerfCtxBL6,vecPerfCtxAlb,vecPerfNOTBL6,vecPerfNOTAlb);
vecG1 = cat(1,ones(size(vecPerfCtxBL6)),2*ones(size(vecPerfCtxAlb)),ones(size(vecPerfNOTBL6)),2*ones(size(vecPerfNOTAlb))); %bl6 vs alb
cellG1 = cellSubjectGroups(vecG1)';
vecG2 = cat(1,ones(size(vecPerfCtxBL6)),ones(size(vecPerfCtxAlb)),2*ones(size(vecPerfNOTBL6)),2*ones(size(vecPerfNOTAlb))); %ctx vs not
cellG2 = cellAreaGroupsAbbr(vecG2)';

[p,tbl,stats,terms] = anovan(vecY,{cellG1,cellG2},'model','full','display','off');

%% plot
figure
errorbar([1 2],[mean(vecPerfCtxBL6) mean(vecPerfNOTBL6)],...
	[std(vecPerfCtxBL6)/sqrt(numel(vecPerfCtxBL6)) std(vecPerfNOTBL6)/sqrt(numel(vecPerfNOTBL6))],...
	'x-','Color',vecColBl6,'CapSize',20);
hold on
errorbar([1 2],[mean(vecPerfCtxAlb) mean(vecPerfNOTAlb)],...
	[std(vecPerfCtxAlb)/sqrt(numel(vecPerfCtxAlb)) std(vecPerfNOTAlb)/sqrt(numel(vecPerfNOTAlb))],...
	'x-','Color',vecColAlb,'CapSize',20);
plot([1 2],[1 1]*dblChanceP,'--','Color',[0.5 0.5 0.5]);
hold off
ylabel(sprintf('%d-tuple decoding accuracy (fraction)',intMinCells));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
legend(cellSubjectGroups,'Location','best');
title(sprintf('Ctx,p=%.4f,NOT,p=%.4f, Interaction,p=%.4f',p_Ctx,p_NOT,p(3)));
fixfig;grid off;

%% pseudo population BL6Ctx
%calculate min rep
intMinRep = inf;

for intRec=1:numel(cellAllOriVecsBL6{1})
	[varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(cellAllOriVecsBL6{1}{intRec});
	intMinRep = min([intMinRep; vecCounts]);
end
vecPseudoOri = repmat(vecUnique(:),[intMinRep 1]);
[vecPseudoStimIdx,vecPseudoUnique,vecPseudoCounts,cellPseudoSelect,vecPseudoRepetition] = label2idx(vecPseudoOri);
intStimNum = numel(vecUnique);
matPseudoData = nan(0,numel(vecPseudoOri));
for intRec=1:numel(cellAllOriVecsBL6{1})
	[varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(cellAllOriVecsBL6{1}{intRec});
	matTemp = nan(size(cellAllSubDataBL6{1}{intRec},1),intMinRep*numel(vecUnique));
	for intRep=1:intMinRep
		vecUseTrials = find(vecRepetition==intRep);
		[dummy,vecReorder] = sort(cellAllOriVecsBL6{1}{intRec}(vecUseTrials));
		matTemp(:,(1:intStimNum)+intStimNum*(intRep-1)) = cellAllSubDataBL6{1}{intRec}(:,vecUseTrials(vecReorder));
	end
	matPseudoData((end+1):(end+size(matTemp,1)),:) = matTemp;
end
matPseudoData(range(matPseudoData,2)==0,:) = [];
%randomize repetition per neuron
matRandPseudoData = matPseudoData;
for intNeuron=1:size(matPseudoData,1)
	for intStimType=1:intStimNum
		vecTrialTemp = find(vecPseudoStimIdx==intStimType);
		vecAssignRandIdx = vecTrialTemp(randperm(numel(vecTrialTemp)));
		matRandPseudoData(intNeuron,vecTrialTemp) = matPseudoData(intNeuron,vecAssignRandIdx);
	end
end
%check if means per stim type are unaltered
[matRandRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matRandPseudoData,vecPseudoOri);
[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matPseudoData,vecPseudoOri);
matMeanRandR = mean(matRandRespNSR,3);
matMeanR = mean(matRespNSR,3);
dblR=corr(matMeanRandR(:),matMeanR(:));
if dblR~=1
	error('something is wrong')
end
subplot(2,3,1)
imagesc(corr(matPseudoData'))
subplot(2,3,2)
imagesc(corr(matRandPseudoData'))

%%
cellAllSubDataBL6
cellAllOriVecsBL6

cellAllSubDataAlb
cellAllOriVecsAlb