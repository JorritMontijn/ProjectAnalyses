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
intMinCells = 1;%NOT+interaction significant at (1-off diag): 7, 8, 9, 10; not sign. at 6; only diag: 6,7,8 (only interaction), 9,10 (both)
intStimNr = 24;
matDiag = diag(diag(true(intStimNr,intStimNr)));
matIsCorrect = circshift(matDiag,-1) | matDiag | circshift(matDiag,1);
dblChanceP = sum(matIsCorrect(:))/numel(matIsCorrect);


vecCorrectRWt = [];
vecPropOnDiagWt = [];
matAggConfWt = zeros(intStimNr,intStimNr);
matAggConf2Wt = zeros(intStimNr,intStimNr);
vecCorrectRAlb = [];
vecPropOnDiagAlb = [];
matAggConfAlb = zeros(intStimNr,intStimNr);
matAggConf2Alb = zeros(intStimNr,intStimNr);


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
			fprintf('Running %sB%d (rec %d/%d for %s) [%s]\n',strName,intBlock,intRecIdx,numel(vecRunRecs),strSubjectType,getTime);
			
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
			indRemTrials = false;%vecProbCorrectCtx == 0 | isnan(vecProbCorrectCtx) | vecProbCorrectNOT == 0 | isnan(vecProbCorrectNOT);
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
			title(sprintf('%s - Ctx/NOT coding correlation, r=%.3f, p=%.3f',strName,r(1,2),p(1,2)),'interpreter','none');
			xlabel('Cortex, Ori decoding P(correct)');
			ylabel('NOT, Ori decoding P(correct)');
			fixfig;
			
			%select only trials where both are wrong
			indBothWrong = vecDecodedIndexCV_LR_Ctx~=vecTrialTypes & vecDecodedIndexCV_LR_NOT~=vecTrialTypes;
			subplot(2,3,2)
			matConfusion = getFillGrid(zeros(intStimNr),vecDecodedIndexCV_LR_Ctx(indBothWrong),vecDecodedIndexCV_LR_NOT(indBothWrong),ones(intTrials,1));
			imagesc(matConfusion);
			matOriX = repmat(vecUnique,[1 size(vecUnique,1)]);
			matOriY = matOriX';
			matOriDiff = abs(circ_dist(deg2rad(matOriX),deg2rad(matOriY)));
			matUseDiag = matOriDiff < (0.5/intStimNr)*2*pi;
			dblOnDiag = sum(matConfusion(matUseDiag))/sum(matConfusion(:));
			dblChanceDiag = sum(matUseDiag(:))/(intStimNr.^2);
			%dblOnDiag = sum(matConfusion(matIsCorrect))/sum(matConfusion(:));
			%dblChanceDiag = dblChanceP;
			[phat,pci]=binofit(dblOnDiag*sum(matConfusion(:)),sum(matConfusion(:)));
			pBino=myBinomTest(dblOnDiag*sum(matConfusion(:)),sum(matConfusion(:)),dblChanceDiag);
			
			axis xy;
			title(sprintf('Orientation error similarity; n,NOT=%d,Ctx=%d',numel(vecSelectCellsNOT),numel(vecSelectCellsCtx)));
			xlabel('Cortex, decoded ori on error');
			ylabel('NOT, decoded ori on error');
			fixfig;grid off;
			
			subplot(2,3,3)
			hold on
			errorbar(1,phat,phat-pci(1),phat-pci(2),'bx','capsize',50)
			plot([0.5 1.5],dblChanceDiag*[1 1],'--','color',[0.5 0.5 0.5])
			hold off
			fixfig;grid off
			title(sprintf('P(Same error) +/- 95-CI, -- chance; Bino-test.p=%.3f',pBino));
			
			%ori diff
			vecOriDiff = round(rad2deg(matOriDiff(:)));
			vecErrorCounts = matConfusion(:);
			[vecTrialTypesDiff,vecUniqueDiff,vecCountsDiff] = val2idx(vecOriDiff);
			
			vecCountsPerOri= accumarray(vecTrialTypesDiff,vecErrorCounts);
			vecErrorFracPerOri = (vecCountsPerOri ./vecCountsDiff);
			vecErrorFracPerOri = vecErrorFracPerOri ./ sum(vecErrorFracPerOri(:));
			[R,P] = corrcoef(vecUniqueDiff,vecErrorFracPerOri);
			
			subplot(2,3,4)
			hold on
			scatter(vecUniqueDiff,vecErrorFracPerOri)
			hold off
			ylabel('Norm. fraction of errors');
			xlabel('Orientation difference Ctx/NOT (\Deltadeg)');
			set(gca,'xtick',vecUniqueDiff(1:2:end));
			xlim([0 max(vecUniqueDiff)]);
			title(sprintf('R=%.3f, p=%.3e',R(1,2),P(1,2)));
			fixfig;grid off
			maxfig;drawnow;
			
			%get error in NOT as function of error in Ctx
			vecRealOri = vecUnique(vecTrialTypes);
			vecCtxOri = vecUnique(vecDecodedIndexCV_LR_Ctx);
			vecNotOri = vecUnique(vecDecodedIndexCV_LR_NOT);
			
			vecCtxError = roundi(rad2deg(circ_dist(deg2rad(vecRealOri),deg2rad(vecCtxOri))),6);
			vecNotError = roundi(rad2deg(circ_dist(deg2rad(vecRealOri),deg2rad(vecNotOri))),6);
			vecCtxError(vecCtxError==180)=-180;
			vecNotError(vecNotError==180)=-180;
			
			vecErrorIdx = -180:median(diff(vecUnique)):179;
			[a,vecCtxErrorIdx]=ismember(vecCtxError,vecErrorIdx);
			[a,vecNotErrorIdx]=ismember(vecNotError,vecErrorIdx);
			
			%2nd comparison
			matConfusion2WithCorrect = getFillGrid(zeros(numel(vecErrorIdx)),vecCtxErrorIdx,vecNotErrorIdx,ones(intTrials,1));
			%remove trials where one is correct
			matConfusion2 = matConfusion2WithCorrect;
			matConfusion2(vecErrorIdx==0,:)=0;
			matConfusion2(:,vecErrorIdx==0)=0;
			matOriX2 = repmat(vecErrorIdx',[1 size(vecErrorIdx',1)])';
			matOriY2 = matOriX2';
			matOriDiff2 = abs(circ_dist(deg2rad(matOriX2),deg2rad(matOriY2)));
			matUseDiag2 = matOriDiff2 < (0.5/intStimNr)*2*pi;
			dblOnDiag;
			dblChanceDiag;
			dblOnDiag2 = sum(matConfusion2(matUseDiag2))/sum(matConfusion2(:));
			dblChanceDiag2 = sum(matUseDiag2(:))/numel(matUseDiag2);
			[phat2,pci2]=binofit(dblOnDiag2*sum(matConfusion2(:)),sum(matConfusion2(:)));
			pBino2=myBinomTest(dblOnDiag2*sum(matConfusion2(:)),sum(matConfusion2(:)),dblChanceDiag2);
			[dummy,vecErrorDegNot,vecErrorCountsNot] = val2idx(vecNotError);
			[dummy,vecErrorDegCtx,vecErrorCountsCtx] = val2idx(vecCtxError);
			
			
			subplot(2,3,5)
			imagesc(vecErrorIdx,vecErrorIdx,matConfusion2WithCorrect);
			axis xy
			colorbar
			ylabel('Ctx decoding error (degs)');
			xlabel('Not decoding error (degs)');
			set(gca,'xtick',vecErrorDegNot(1:6:end));
			set(gca,'ytick',vecErrorDegNot(1:6:end));
			fixfig;grid off
			title(sprintf('method 2, P(Same error) vs chance; Bino-test p=%.3f',pBino2));
			
			subplot(2,3,6)
			hold on
			plot(vecErrorDegCtx,vecErrorCountsCtx,'b');
			plot(vecErrorDegNot,vecErrorCountsNot,'r');
			set(gca,'xtick',vecErrorDegNot(1:6:end));
			hold off
			fixfig;grid off;
			legend({'Ctx','NOT'})
			ylabel('Number of trials (count)');
			xlabel('Decoding error (degs)');
			
			%save fig
			drawnow;
			export_fig([strTargetPath filesep 'single_recs' filesep sprintf('OriDecodingErrors%sB%d.tif',strName,intBlock)]);
			saveas(gcf,[strTargetPath filesep 'single_recs' filesep sprintf('OriDecodingErrors%sB%d.pdf',strName,intBlock)]);
			
			%save data
			if strcmp(strSubjectType,'BL6')
				matAggConfWt = matAggConfWt + matConfusion;
				matAggConf2Wt = matAggConfWt + matConfusion2WithCorrect;
				vecCorrectRWt(end+1) = R(1,2);
				vecPropOnDiagWt(end+1) = dblOnDiag;
			else
				matAggConfAlb = matAggConfAlb + matConfusion;
				matAggConf2Alb = matAggConf2Alb + matConfusion2WithCorrect;
				vecCorrectRAlb(end+1) = R(1,2);
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

subplot(2,3,3)
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


subplot(2,3,4)
matAggConf2WtNoDiag = matAggConf2Wt;
matAggConf2WtNoDiag(vecErrorIdx==0,:)=nan;
matAggConf2WtNoDiag(:,vecErrorIdx==0)=nan;		
imagesc(vecErrorIdx,vecErrorIdx,matAggConf2WtNoDiag);
%nancolorbar(matAggConf2WtNoDiag,[min(matAggConf2WtNoDiag(:)) max(matAggConf2WtNoDiag(:))],'parula',[1 1 1]);
axis xy;
title('Error matrix BL6, v2');
ylabel('Ctx decoding error (degs)');
xlabel('Not decoding error (degs)');
set(gca,'xtick',vecErrorDegNot(1:6:end));
set(gca,'ytick',vecErrorDegNot(1:6:end));
fixfig;grid off

subplot(2,3,5)
matAggConf2AlbNoDiag = matAggConf2Alb;
matAggConf2AlbNoDiag(vecErrorIdx==0,:)=nan;
matAggConf2AlbNoDiag(:,vecErrorIdx==0)=nan;
imagesc(vecErrorIdx,vecErrorIdx,matAggConf2AlbNoDiag);
%nancolorbar(matAggConf2AlbNoDiag,[min(matAggConf2AlbNoDiag(:)) max(matAggConf2AlbNoDiag(:))],'parula',[1 1 1]);
axis xy;
title('Error matrix DBA, v2');
ylabel('Ctx decoding error (degs)');
xlabel('Not decoding error (degs)');
set(gca,'xtick',vecErrorDegNot(1:6:end));
set(gca,'ytick',vecErrorDegNot(1:6:end));
fixfig;grid off

%
%error add lines that envelop the diagonal?

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('OriDecodingErrors.tif')]);
saveas(gcf,[strTargetPath filesep sprintf('OriDecodingErrors.pdf')]);
