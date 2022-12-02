%% creates ori decoding figs, including pseudo pops
%[done/iii) decoding left vs right: is DBA worse than BL6?

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

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);



%% pre-allocate
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
			
			% get pupil data
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
			
			% prep data
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
			
			% split cells into areas
			%build cell vectors
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			end
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			for intArea = 1:numel(cellUseAreas)
				% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				%calc tuning curve
				matUseData = matData(vecSelectCells,:);
				sOut = getTuningCurves(matUseData,vecOrientation);
				vecTuningP_A = sOut.vecOriAnova;
				vecTuningP_R2 = sOut.vecFitP;
				[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
				vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
				indTunedCells = vecTuningP_R2_corr(:)'<0.05;
				
				%get RLR
				vecR = mean(sOut.matFittedResp(:,ismember(sOut.vecUniqueDegs,0)),2);
				vecL = mean(sOut.matFittedResp(:,ismember(sOut.vecUniqueDegs,180)),2);
				vecRLR = vecR ./ (vecR+vecL);
				
				%remove range 0
				vecRangeHz = range(matUseData,2)';
				indRem2=vecRangeHz==0 | ~indTunedCells | isnan(vecRLR(:)');
				matUseData(indRem2,:)=[];
				vecRLR(indRem2)=[];
				vecSelectCells = vecSelectCells(~indRem2);
				
				%assign data
				if strcmp(strSubjectType,'BL6')
					cellAllSubDataBL6{intArea}(end+1) = {matUseData};
					cellAllOriVecsBL6{intArea}(end+1) = {vecOrientation};
				else
					cellAllSubDataAlb{intArea}(end+1) = {matUseData};
					cellAllOriVecsAlb{intArea}(end+1) = {vecOrientation};
				end
			end
		end
	end
end

%% pseudo populations
intTypeCV = 2;
dblLambda = 1000;
intUseCells = 20;
intUseReps = 14;
intRunPerms = 2500;
cellPref = cell(2,2);
matCellNum = nan(2,2);
matPermPerf = nan(2,2,intRunPerms);
for intArea=1:2
	strArea = cellAreaGroupsAbbr{intArea}
	for intSubjectType=1:2
		strSubjectType = cellSubjectGroups{intSubjectType}
		if intSubjectType == 1
			cellData = cellAllSubDataBL6{intArea};
			cellLabels = cellAllOriVecsBL6{intArea};
		else
			cellData = cellAllSubDataAlb{intArea};
			cellLabels = cellAllOriVecsAlb{intArea};
		end
		%get pseudo data
		[vecPseudoOri,matWhitePseudoData,matPseudoData] = buildPseudoData(cellLabels,cellData,intUseReps);
		[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecPseudoOri);
		vecPriorDistribution = vecCounts;
		
		%check if means per stim type are unaltered
		[matRandRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matWhitePseudoData,vecPseudoOri);
		[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matPseudoData,vecPseudoOri);
		matMeanRandR = mean(matRandRespNSR,3);
		matMeanR = mean(matRespNSR,3);
		dblR=corr(matMeanRandR(:),matMeanR(:));
		if roundi(dblR,10)~=1
			figure
			subplot(2,3,1)
			imagesc(corr(matPseudoData'))
			subplot(2,3,2)
			imagesc(corr(matWhitePseudoData'))
			error('something is wrong')
		end
		
		%run analysis
		sOut = getTuningCurves(matWhitePseudoData,vecPseudoOri);
		vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
		[h,d,vecP] = fdr_bh(sOut.vecFitP);
		indKeepCells = vecP<0.05;
		vecPrefOri(~indKeepCells)=[];
		cellPref{intArea} = vecPrefOri;
		matCellNum(intArea,intSubjectType) = sum(indKeepCells);
		
		%% decode all
		vecSubSelect = find(indKeepCells);
		%run multiple permutations and average
		intNumN = numel(vecSubSelect);
		intChoosePerms = min([intRunPerms nchoosek(intNumN,min([intUseCells intNumN]))-1]);
		vecPermPerf = nan(1,intChoosePerms);
		boolChosen = false;
		intChooseIdx=1;
		matPerms = nan(intChoosePerms,intUseCells);
		while ~boolChosen
			vecRandPerm = sort(randperm(intNumN,intUseCells));
			if ~any(all(bsxfun(@eq,matPerms,vecRandPerm),2))
				matPerms(intChooseIdx,:) = vecRandPerm;
				intChooseIdx = intChooseIdx + 1;
			end
			if intChooseIdx > intChoosePerms
				boolChosen = true;
			end
		end
		
		vecUseTrials = vecPseudoOri==0 | vecPseudoOri==180;
		vecPseudoOriLR = vecPseudoOri(vecUseTrials);
		[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecPseudoOriLR);
		vecPriorDistroLR = vecCounts;
		
		for intPerm=1:intChoosePerms
			vecUseRandPerm = matPerms(intPerm,:);
			vecUseSubSelect = vecSubSelect(vecUseRandPerm);
			
			matUseSubData = matWhitePseudoData(vecUseSubSelect,vecUseTrials);
			
			
			%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecOriIdx,intTypeCV,vecPriorDistribution);
			%vecPerfTM(intArea) = dblPerformanceCV;
			[dblPerformanceLR,vecDecodedIndexCV_LR,matPosteriorProbability,dblMeanErrorDegsLR,matConfusionLR,matWeights] = ...
				doCrossValidatedDecodingLR(matUseSubData,vecPseudoOriLR,intTypeCV,vecPriorDistroLR,dblLambda);
			vecPermPerf(intPerm) = dblPerformanceLR;
		end
		matPermPerf(intArea,intSubjectType,:) = vecPermPerf;
	end
end

%% plot
vecLowHighCtxBL6 = getCI(squeeze(matPermPerf(1,1,:)),1,0.05,true);
vecLowHighNotBL6 = getCI(squeeze(matPermPerf(2,1,:)),1,0.05,true);
vecLowHighCtxDBA = getCI(squeeze(matPermPerf(1,2,:)),1,0.05,true);
vecLowHighNotDBA = getCI(squeeze(matPermPerf(2,2,:)),1,0.05,true);


figure;maxfig;
intPoints = intUseReps;
vecBinCenters = linspace(0,1,intPoints+1);
vecBinEdges = [vecBinCenters(1)-1/intPoints vecBinCenters] + 0.5/intPoints;

subplot(2,3,1)
d=getCohensD(flat(matPermPerf(1,1,:)),flat(matPermPerf(1,2,:)));
p_Ctx=1-(normcdf(d,0,1)-normcdf(-d,0,1));

vecCountsCtxBL6 = histcounts(matPermPerf(1,1,:),vecBinEdges);
vecCountsCtxDBA = histcounts(matPermPerf(1,2,:),vecBinEdges);
hold on
plot(vecBinCenters,vecCountsCtxBL6,'color',vecColBl6);
plot(vecBinCenters,vecCountsCtxDBA,'color',vecColAlb);
hold off
legend({'BL6','DBA'},'location','best');
xlabel('Fraction correct');
ylabel('Count (bootstrap sample)');
title(sprintf('Pseudo-pop Ctx, L vs R decoding, z-test,p=%.3f',p_Ctx));
fixfig;

subplot(2,3,2)
d=getCohensD(flat(matPermPerf(2,1,:)),flat(matPermPerf(2,2,:)));
p_NOT=normcdf(d,0,1,'upper')/2;

vecCountsNotBL6 = histcounts(matPermPerf(2,1,:),vecBinEdges);
vecCountsNotDBA = histcounts(matPermPerf(2,2,:),vecBinEdges);
hold on
plot(vecBinCenters,vecCountsNotBL6,'color',vecColBl6);
plot(vecBinCenters,vecCountsNotDBA,'color',vecColAlb);
hold off
legend({'BL6','DBA'},'location','best');
xlabel('Fraction of trials correctly decoded');
ylabel('Count (bootstrap sample)');
title(sprintf('Pseudo-pop NOT, L vs R decoding, z-test,p=%.1e',p_NOT));
fixfig;

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('OriDecodingLeftRight%dperms.tif',intChoosePerms)]);
export_fig([strTargetPath filesep sprintf('OriDecodingLeftRight%dperms.pdf',intChoosePerms)]);
	
