%% creates ori decoding figs, including pseudo pops
%[done/1) how do confusion matrices of ori decoding differ between alb/bl6?

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
intMinCells = 3;%NOT+interaction significant at (1-off diag): 7, 8, 9, 10; not sign. at 6; only diag: 6,7,8 (only interaction), 9,10 (both)
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
fixfig;grid off;drawnow;
export_fig([strTargetPath filesep sprintf('OriDecoding%dtuple.tif',intMinCells)]);
export_fig([strTargetPath filesep sprintf('OriDecoding%dtuple.pdf',intMinCells)]);
	
%% pseudo populations
figure;maxfig;
intUseCells = 20;
intUseReps = 14;
cellPref = cell(2,2);
matCellNum = nan(2,2);
matMeanPerf = nan(2,2);
for intArea=1:2
	strArea = cellAreaGroupsAbbr{intArea};
	for intSubjectType=1:2
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
		intRunPerms = 50;
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
		
		matAggConfusion = zeros(numel(vecUnique),numel(vecUnique));
		for intPerm=1:intChoosePerms
			fprintf('Permutation %d/%d [%s]\n',intPerm,intChoosePerms,getTime);
			vecUseRandPerm = matPerms(intPerm,:);
			vecUseSubSelect = vecSubSelect(vecUseRandPerm);
			matUseSubData = matWhitePseudoData(vecUseSubSelect,:);
			%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecOriIdx,intTypeCV,vecPriorDistribution);
			%vecPerfTM(intArea) = dblPerformanceCV;
			[dblPerformanceLR,vecDecodedIndexCV_LR,matPosteriorProbability,dblMeanErrorDegsLR,matConfusionLR,matWeights] = ...
				doCrossValidatedDecodingLR(matUseSubData,vecPseudoOri,intTypeCV,vecPriorDistribution,dblLambda);
			vecPermPerf(intPerm) = sum(matConfusionLR(matIsCorrect))/sum(matConfusionLR(:));
			matAggConfusion = matAggConfusion + matConfusionLR;
		end
		matMeanPerf(intArea,intSubjectType) = mean(vecPermPerf);
		intPlotLoc = (intSubjectType-1)*3 + intArea;
		subplot(2,3,intPlotLoc)
		imagesc(vecUnique,vecUnique,matConfusionLR,[0 max(matConfusionLR(:))]);
		colormap(redwhite);
		xlabel('Real stim. ori. (degs)')
		ylabel('Decoded stim. ori. (degs)');
		set(gca,'xtick',0:45:360);
		set(gca,'ytick',0:45:360);
		title(sprintf('Example pseudo pop %s %s',cellSubjectGroups{intSubjectType},cellAreaGroupsAbbr{intArea}));
		axis xy;
		fixfig;
		grid off
		drawnow;
	end
end

%% plot means
n = numel(vecPseudoOri);
[p_CtxWt_vs_CtxAlb] = bino2test(round(matMeanPerf(1,1)*n),n,round(matMeanPerf(1,2)*n),n);
[p_NotWt_vs_NotAlb] = bino2test(round(matMeanPerf(2,1)*n),n,round(matMeanPerf(2,2)*n),n,true);
[p_CtxWt_vs_NotWt] = bino2test(round(matMeanPerf(1,1)*n),n,round(matMeanPerf(2,1)*n),n,true);
[p_CtxAlb_vs_NotAlb] = bino2test(round(matMeanPerf(1,2)*n),n,round(matMeanPerf(2,2)*n),n,true);

[dblM_CtxWt, vecCI_CtxWt] = binofit(matMeanPerf(1,1)*n,n);
[dblM_CtxAlb, vecCI_CtxAlb] = binofit(matMeanPerf(1,2)*n,n);
[dblM_NotWt, vecCI_NotWt] = binofit(matMeanPerf(2,1)*n,n);
[dblM_NotAlb, vecCI_NotAlb] = binofit(matMeanPerf(2,2)*n,n);

subplot(2,3,3)
errorbar([1 2],[dblM_CtxWt dblM_NotWt],...
	[dblM_CtxWt-vecCI_CtxWt(1) dblM_CtxWt-vecCI_CtxWt(2)],...
	[dblM_NotWt-vecCI_NotWt(1) dblM_NotWt-vecCI_NotWt(2)],...
	'x-','Color',vecColBl6,'CapSize',20);
hold on
errorbar([1 2],[dblM_CtxAlb dblM_NotAlb],...
	[dblM_CtxAlb-vecCI_CtxAlb(1) dblM_CtxAlb-vecCI_CtxAlb(2)],...
	[dblM_NotAlb-vecCI_NotAlb(1) dblM_NotAlb-vecCI_NotAlb(2)],...
	'x-','Color',vecColAlb,'CapSize',20);
plot([1 2],[1 1]*dblChanceP,'--','Color',[0.5 0.5 0.5]);
hold off
ylabel(sprintf('Pseudo pop decoding accuracy (fraction)'));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
legend(cellSubjectGroups,'Location','best');
title(sprintf('C,p=%.3f,N,p=%.1e,Wt C-N,p=%.3f,A C-N,p=%.1e',p_CtxWt_vs_CtxAlb,p_NotWt_vs_NotAlb,p_CtxWt_vs_NotWt,p_CtxAlb_vs_NotAlb));
fixfig;grid off;

subplot(2,3,6)
imagesc(vecUnique,vecUnique,matIsCorrect);
colorbar
xlabel('Real stim. ori. (degs)')
ylabel('Decoded stim. ori. (degs)');
set(gca,'xtick',0:45:360);
set(gca,'ytick',0:45:360);
title(sprintf('"Correct"'));
axis xy;
fixfig;
grid off
drawnow;
export_fig([strTargetPath filesep sprintf('OriDecodingPseudoPop.tif')]);
export_fig([strTargetPath filesep sprintf('OriDecodingPseudoPop.pdf')]);
	