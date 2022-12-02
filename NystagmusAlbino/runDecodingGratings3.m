%% creates ori decoding figs, including pseudo pops
%3) does info in NOT predict info in V1? [is this doable? sufficient twin recording??]

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


vecCorrectRWt = [];
vecPropOnDiagWt = [];
matAggConfWt = zeros(intStimNr,intStimNr);
vecCorrectRAlb = [];
vecPropOnDiagAlb = [];
matAggConfAlb = zeros(intStimNr,intStimNr);


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
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx);
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
			if numel(vecUnique) ~= 24 || numel(vecRepetition) < 100,continue,end
			
			%get data matrix
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			dblDur = median(vecStimOffTime-vecStimOnTime);
			matData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblDur);
			
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			vecUseCells = find(indUseCells);
			
			%% split cells into areas
			%build cell vectors
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster(:),cellUseAreas{intArea},'IgnoreCase',true);
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
			dblLambda = 100;
			
			%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
			[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			vecPriorDistribution = vecCounts;
			intStimNr = numel(vecUnique);
			
			%select cells in Ctx and NOT
			vecSelectCellsCtx = find(indUseCells(:) & cellCellsPerArea{1}(:));
			vecSelectCellsNOT = find(indUseCells(:) & cellCellsPerArea{2}(:));
			dblMinCells = intMinCells;
			if numel(vecSelectCellsCtx) < dblMinCells || numel(vecSelectCellsNOT) < dblMinCells
				fprintf('Skipping %sB%d: Ctx=%d, NOT=%d\n',strName,intBlock,numel(vecSelectCellsCtx),numel(vecSelectCellsNOT));
				continue;
			else
				fprintf('%s; %sB%d; Ctx=%d, NOT=%d\n',strSubjectType,strName,intBlock,numel(vecSelectCellsCtx),numel(vecSelectCellsNOT));
			end
			
			%% decode Ctx
			[dblPerformanceLR_Ctx,vecDecodedIndexCV_LR_Ctx,matPosteriorProbability_Ctx,dblMeanErrorDegsLR_Ctx,matConfusionLR_Ctx,matWeights_Ctx] = ...
				doCrossValidatedDecodingLR(matData(indUseCells(:) & cellCellsPerArea{1},:),vecTrialTypes,intTypeCV,vecPriorDistribution,dblLambda);
			pBinoCtx=myBinomTest(dblPerformanceLR_Ctx*sum(matConfusionLR_Ctx(:)),sum(matConfusionLR_Ctx(:)),1/intStimNr);
			
			%% decode NOT
			[dblPerformanceLR_NOT,vecDecodedIndexCV_LR_NOT,matPosteriorProbability_NOT,dblMeanErrorDegsLR_NOT,matConfusionLR_NOT,matWeights_NOT] = ...
				doCrossValidatedDecodingLR(matData(indUseCells & cellCellsPerArea{2},:),vecTrialTypes,intTypeCV,vecPriorDistribution,dblLambda);
			pBinoNOT=myBinomTest(dblPerformanceLR_NOT*sum(matConfusionLR_NOT(:)),sum(matConfusionLR_NOT(:)),1/intStimNr);
			if (dblPerformanceLR_NOT < 1/intStimNr || pBinoNOT>0.05) && (dblPerformanceLR_Ctx < 1/intStimNr || pBinoCtx>0.05),continue;end
			
			%% get prob correct per trial
			intTrials = numel(vecDecodedIndexCV_LR_Ctx);
			vecProbCorrectCtx = nan(1,intTrials);
			vecProbCorrectNOT = nan(1,intTrials);
			for intTrial=1:intTrials
				%get probability
				vecProbCorrectCtx(intTrial) = matPosteriorProbability_Ctx(vecTrialTypes(intTrial),intTrial);
				vecProbCorrectNOT(intTrial) = matPosteriorProbability_NOT(vecTrialTypes(intTrial),intTrial);
			end
			
			%remove trials
			indRemTrials = vecProbCorrectCtx == 0 | isnan(vecProbCorrectCtx) | vecProbCorrectNOT == 0 | isnan(vecProbCorrectNOT);
			vecProbCorrectCtx(indRemTrials) = [];
			vecProbCorrectNOT(indRemTrials) = [];
			vecDecodedIndexCV_LR_Ctx(indRemTrials) = [];
			vecDecodedIndexCV_LR_NOT(indRemTrials) = [];
			vecTrialTypes(indRemTrials) = [];
			
			%plot
			%close;
			figure
			subplot(2,3,1)
			[r,p,ul,ll]=corrcoef(vecProbCorrectNOT(:),vecProbCorrectCtx(:));
			scatter(vecProbCorrectNOT(:),vecProbCorrectCtx(:),'.')
			title(sprintf('%s - Ctx/NOT coding correlation, r=%.3f, p=%.3f',strName,r(1,2),p(1,2)));
			xlabel('Cortex, Ori decoding P(correct)');
			ylabel('NOT, Ori decoding P(correct)');
			fixfig;
			
			%select only trials where both are wrong
			indBothWrong = vecDecodedIndexCV_LR_Ctx~=vecTrialTypes & vecDecodedIndexCV_LR_NOT~=vecTrialTypes;
			subplot(2,3,2)
			matConfusion = getFillGrid(zeros(intStimNr),vecDecodedIndexCV_LR_Ctx(indBothWrong),vecDecodedIndexCV_LR_NOT(indBothWrong),ones(intTrials,1));
			imagesc(matConfusion);
			dblOnDiag = sum(diag(matConfusion))/sum(matConfusion(:));
			dblChanceDiag = intStimNr/(intStimNr.^2);
			%dblOnDiag = sum(matConfusion(matIsCorrect))/sum(matConfusion(:));
			%dblChanceDiag = dblChanceP;
			[phat,pci]=binofit(dblOnDiag*sum(matConfusion(:)),sum(matConfusion(:)));
			pBino=myBinomTest(dblOnDiag*sum(matConfusion(:)),sum(matConfusion(:)),dblChanceDiag);
			
			axis xy;
			title(sprintf('Orientation error similarity; n,NOT=%d,Ctx=%d',numel(vecSelectCellsNOT),numel(vecSelectCellsCtx)));
			xlabel('Cortex, Ori decoding P(correct)');
			ylabel('NOT, Ori decoding P(correct)');
			fixfig;grid off;
			
			subplot(2,3,3)
			hold on
			errorbar(1,phat,phat-pci(1),phat-pci(2),'bx','capsize',50)
			plot([0.5 1.5],dblChanceDiag*[1 1],'--','color',[0.5 0.5 0.5])
			hold off
			fixfig;grid off
			title(sprintf('P(Same error) +/- 95-CI, -- chance; Bino-test.p=%.3f',pBino));
			maxfig;drawnow;
			
			if strcmp(strSubjectType,'BL6')
				matAggConfWt = matAggConfWt + matConfusion;
				vecCorrectRWt(end+1) = r(1,2);
				vecPropOnDiagWt(end+1) = dblOnDiag;
			else
				matAggConfAlb = matAggConfAlb + matConfusion;
				vecCorrectRAlb(end+1) = r(1,2);
				vecPropOnDiagAlb(end+1) = dblOnDiag;
			end
		end
	end
	
	
	
end

vecCorrectRWt
vecPropOnDiagWt

vecCorrectRAlb
vecPropOnDiagAlb

%% test
[h,p_correctr]=ttest2(vecCorrectRWt,vecCorrectRAlb);
[h,p_correctrAlb]=ttest(vecCorrectRAlb,0);
[h,p_correctrWt]=ttest(vecCorrectRWt,0);

[h,p_sameoriWt]=ttest(vecPropOnDiagWt,dblChanceDiag);
[h,p_sameoriAlb]=ttest(vecPropOnDiagAlb,dblChanceDiag);

%% plot
figure
subplot(2,3,1)
hold on
errorbar(1,mean(vecCorrectRWt),std(vecCorrectRWt)/sqrt(numel(vecCorrectRWt)),'xb')
errorbar(2,mean(vecCorrectRAlb),std(vecCorrectRAlb)/sqrt(numel(vecCorrectRAlb)),'xr')
hold off
xlim([0 3]);
title(sprintf('t-test,p=%.3f, vs 0; Wt,p=%.3f, Alb,p=%.3f',p_correctr,p_correctrWt,p_correctrAlb))
set(gca,'xtick',[1 2],'xticklabel',{'BL6','DBA'});
ylabel('Correlation correct R(P(Ctx),P(NOT))');
fixfig;grid off;

subplot(2,3,2)
hold on
errorbar(1,mean(vecPropOnDiagWt),std(vecPropOnDiagWt)/sqrt(numel(vecPropOnDiagWt)),'xb')
errorbar(2,mean(vecPropOnDiagAlb),std(vecPropOnDiagAlb)/sqrt(numel(vecPropOnDiagAlb)),'xr')
plot([1 2],[1 1]*dblChanceDiag,'--','color',[0.5 0.5 0.5])
set(gca,'xtick',[1 2],'xticklabel',{'BL6','DBA'});
ylabel('Same-ori decoding error (Ctx,NOT)');
xlim([0 3]);
title(sprintf('t-test vs 0; Wt,p=%.3f; Alb,p=%.3f',p_sameoriWt,p_sameoriAlb))
hold off
fixfig;grid off
maxfig;

subplot(2,3,5)
imagesc(vecUnique,vecUnique,matAggConfWt);
axis xy;
title('Error matrix BL6');
ylabel('Ori. decoded from cortex (degs)');
xlabel('Ori. decoded from NOT (degs)');
set(gca,'xtick',0:45:360);
set(gca,'ytick',0:45:360);
fixfig;grid off

subplot(2,3,6)
imagesc(vecUnique,vecUnique,matAggConfAlb);
axis xy;
title('Error matrix DBA');
ylabel('Ori. decoded from cortex (degs)');
xlabel('Ori. decoded from NOT (degs)');
set(gca,'xtick',0:45:360);
set(gca,'ytick',0:45:360);
fixfig;grid off

%
%error add lines that envelop the diagonal?

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('OriDecodingErrors.tif')]);
export_fig([strTargetPath filesep sprintf('OriDecodingErrors.pdf')]);
	
