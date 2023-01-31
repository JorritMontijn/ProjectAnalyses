%{
ori vs dir tuning
%}
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

%% pre-allocate
cellTuningP_AlbCtx = {};
cellOSI_AlbCtx = {};
cellOPI_AlbCtx = {};
cellDSI_AlbCtx = {};
vecOSI_AlbCtx = [];
vecOPI_AlbCtx = [];
vecDSI_AlbCtx = [];
cellPrefOri_AlbCtx = {};

cellTuningP_Bl6Ctx = {};
cellOSI_Bl6Ctx = {};
cellOPI_Bl6Ctx = {};
cellDSI_Bl6Ctx = {};
vecOSI_Bl6Ctx = [];
vecOPI_Bl6Ctx = [];
vecDSI_Bl6Ctx = [];
cellPrefOri_Bl6Ctx = {};

cellTuningP_Bl6NOT = {};
cellOSI_Bl6NOT = {};
cellOPI_Bl6NOT = {};
cellDSI_Bl6NOT = {};
vecOSI_Bl6NOT = [];
vecOPI_Bl6NOT = [];
vecDSI_Bl6NOT = [];
cellPrefOri_Bl6NOT = {};

cellTuningP_AlbNOT = {};
cellOSI_AlbNOT = {};
cellOPI_AlbNOT = {};
cellDSI_AlbNOT = {};
vecOSI_AlbNOT = [];
vecOPI_AlbNOT = [];
vecDSI_AlbNOT = [];
cellPrefOri_AlbNOT = {};


%% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
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
		vecBlocksNM = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
		%get timing for DG
		intBlock = vecBlocksDG(1);
		sBlock = sRec.cellBlock{intBlock};
		if isfield(sBlock,'vecPupilStimOn')
			vecPupilLatency = sBlock.vecPupilStimOn-sBlock.vecStimOnTime;
		else
			vecPupilLatency = 0;
		end
		
		%% concatenate blocks
		vecAggOri = [];
		vecAggStimOnTime = [];
		vecAggStimOffTime = [];
		for intBlockIdx=1:min(numel(vecBlocksDG),2)
			intBlock = vecBlocksDG(intBlockIdx);
			fprintf('Running %sB%d (rec %d/%d for %s) [%s]\n',strName,intBlock,intRecIdx,numel(vecRunRecs),strSubjectType,getTime);
			sBlock = sRec.cellBlock{intBlock};
			intPopCounter = intPopCounter + 1;
			
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
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			
			vecAggOri = cat(1,vecAggOri,vecOrientation);
			vecAggStimOnTime = cat(2,vecAggStimOnTime,vecStimOnTime);
			vecAggStimOffTime = cat(2,vecAggStimOffTime,vecStimOffTime);
		end
		%get data matrix
		cellSpikeT = {sRec.sCluster(:).SpikeTimes};
		dblStimDur = median(vecAggStimOffTime-vecAggStimOnTime);
		[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecAggOri);
		if numel(vecUnique) ~= 24,continue,end
		matData = getSpikeCounts(cellSpikeT,vecAggStimOnTime,dblStimDur);
		
		%include?
		indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
		
		%% split cells into areas
		%build cell vectors
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		
		%% go through areas
		%pre-allocate
		%vecPerf = nan(numel(cellUseAreas),1);
		%matPerf = nan(numel(cellUseAreas),intStimNr);
		%cellPref = cell(numel(cellUseAreas),1);
		
		%params
		intTypeCV = 2; %leave repetition out
		vecOriNoDir = mod(vecAggOri,180);
		vecTrialTypesNoDir = deg2rad(vecOriNoDir)*2;
		[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTrialTypesNoDir);
		dblLambda = 100;
		
		%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
		[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecAggOri);
		vecPriorDistribution = vecCounts;
		intStimNr = numel(vecUnique);
		
		for intArea = 1:2%numel(cellUseAreas)
			%% select cells
			strAreaGroup =  cellAreaGroupsAbbr{intArea};
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			if isempty(vecSelectCells)
				continue;
			end
			
			%% calc tuning curves & zeta
			%calc tuning curve
			matUseData = matData(vecSelectCells,:);
			sOut = getTuningCurves(matUseData,vecAggOri);
			vecTuningP_A = sOut.vecOriAnova;
			vecTuningP_R2 = sOut.vecFitP;
			[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
			%vecTuningP_R2_corr = vecTuningP_R2;
			vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
			indTunedCells = vecTuningP_R2_corr<0.05;
			
			%get OPI
			%matUseR = sOut.matMeanResp(indTunedCells,:);
			matUseR = sOut.matFittedResp;
			vecOri24 = mod(vecOri24,180);
			vecDir24 = sOut.vecUniqueDegs;
			
			%collapse angles
			vecOPI = getOPI(matUseR,vecOri24);
			vecOSI = getOSI(matUseR,vecOri24);
			
			%get DSI
			%intDirNum = numel(unique(vecDir24));
			%[vecPrefR,vecPrefIdx]=max(matUseR,[],2);
			%vecAntiIdx = modx(vecPrefIdx+intDirNum/2,intDirNum);
			%vecAntiR = getMatVals(matUseR,(1:numel(vecAntiIdx))',vecAntiIdx);
			%vecDSI = (vecPrefR-vecAntiR)./(vecAntiR+vecPrefR);
			
			%get DSI horizontal
			vecPrefR = matUseR(:,vecDir24==0);
			vecAntiR = matUseR(:,vecDir24==180);
			vecDSI = (vecPrefR-vecAntiR)./(vecAntiR+vecPrefR);
			
			%{
				figure;
				subplot(2,3,1)
				scatter(vecOPI,vecOSI)
				xlabel('OPI')
				ylabel('OSI')
				
				subplot(2,3,2)
				scatter(vecOPI,vecDSI)
				xlabel('OPI')
				ylabel('DSI')
				
				subplot(2,3,3)
				scatter(vecOSI,vecDSI)
				xlabel('OSI')
				ylabel('DSI')
			%}
			
			if 0
				%% plot
				for intCellIdx=1:numel(vecSelectCells)
					intCell = vecSelectCells(intCellIdx);
					vecSpikeT = cellSpikeT{intCell};
					
					[vecTime,vecRate] = getIFR(vecSpikeT,vecAggStimOnTime,dblTrialDur);
					
					figure
					subplot(2,3,1)
					errorbar(sOut.vecUniqueDegs,sOut.matMeanResp(intCellIdx,:),sOut.matSDResp(intCellIdx,:)./sqrt(vecCounts(:)'));
					xlabel('Stim ori (degs)');
					ylabel('Mean rate (Hz)');
					title(sprintf('tuning-p=%.3f (%.3e)',vecTuningP_R2_corr(intCellIdx),vecTuningP_R2_corr(intCellIdx)));
					
					subplot(2,3,2)
					plot(vecTime,vecRate);
					xlabel('Time (s)');
					ylabel('Instantaneous rate (Hz)');
					title(sprintf('Zeta-p=%.3f (%.3e)',vecZetaP_corr(intCellIdx),vecZetaP_corr(intCellIdx)));
					
					drawnow;
					
				end
			end
			
			%% save data
			if intArea == 1
				if strcmpi(strSubjectType,'DBA')
					cellTuningP_AlbCtx{end+1} = vecTuningP_R2_corr;
					cellOSI_AlbCtx{end+1} = vecOSI;
					cellDSI_AlbCtx{end+1} = vecDSI;
					cellOPI_AlbCtx{end+1} = vecOPI;
					vecOSI_AlbCtx(end+1) = nanmean(vecOSI);
					vecDSI_AlbCtx(end+1) = nanmean(vecDSI);
					vecOPI_AlbCtx(end+1) = nanmean(vecOPI);
					cellPrefOri_AlbCtx{end+1} = vecPrefOri;
				elseif strcmpi(strSubjectType,'Bl6')
					cellTuningP_Bl6Ctx{end+1} = vecTuningP_R2_corr;
					cellOSI_Bl6Ctx{end+1} = vecOSI;
					cellDSI_Bl6Ctx{end+1} = vecDSI;
					cellOPI_Bl6Ctx{end+1} = vecOPI;
					vecOSI_Bl6Ctx(end+1) = nanmean(vecOSI);
					vecDSI_Bl6Ctx(end+1) = nanmean(vecDSI);
					vecOPI_Bl6Ctx(end+1) = nanmean(vecOPI);
					cellPrefOri_Bl6Ctx{end+1} = vecPrefOri;
				end
			elseif intArea == 2
				if strcmpi(strSubjectType,'DBA')
					cellTuningP_AlbNOT{end+1} = vecTuningP_R2_corr;
					cellOSI_AlbNOT{end+1} = vecOSI;
					cellDSI_AlbNOT{end+1} = vecDSI;
					cellOPI_AlbNOT{end+1} = vecOPI;
					vecOSI_AlbNOT(end+1) = nanmean(vecOSI);
					vecDSI_AlbNOT(end+1) = nanmean(vecDSI);
					vecOPI_AlbNOT(end+1) = nanmean(vecOPI);
					cellPrefOri_AlbNOT{end+1} = vecPrefOri;
				elseif strcmpi(strSubjectType,'Bl6')
					cellTuningP_Bl6NOT{end+1} = vecTuningP_R2_corr;
					cellOSI_Bl6NOT{end+1} = vecOSI;
					cellDSI_Bl6NOT{end+1} = vecDSI;
					cellOPI_Bl6NOT{end+1} = vecOPI;
					vecOSI_Bl6NOT(end+1) = nanmean(vecOSI);
					vecDSI_Bl6NOT(end+1) = nanmean(vecDSI);
					vecOPI_Bl6NOT(end+1) = nanmean(vecOPI);
					cellPrefOri_Bl6NOT{end+1} = vecPrefOri;
				end
			end
		end
		
	end
end

%% plot ori vs dir tuning in bl6 vs albino: is direction selectivity in DBA NOT specifically affected?
%{
cellTuningP_AlbCtx = {};
cellOSI_AlbCtx = {};
cellOPI_AlbCtx = {};
cellDSI_AlbCtx = {};
cellPrefOri_AlbCtx = {};

cellTuningP_Bl6Ctx = {};
cellOSI_Bl6Ctx = {};
cellOPI_Bl6Ctx = {};
cellDSI_Bl6Ctx = {};
cellPrefOri_Bl6Ctx = {};
	
cellTuningP_Bl6NOT = {};
cellOSI_Bl6NOT = {};
cellOPI_Bl6NOT = {};
cellDSI_Bl6NOT = {};
cellPrefOri_Bl6NOT = {};

cellTuningP_AlbNOT = {};
cellOSI_AlbNOT = {};
cellOPI_AlbNOT = {};
cellDSI_AlbNOT = {};
cellPrefOri_AlbNOT = {};
%}
figure;maxfig
for intArea=1:2
	if intArea == 1
		%ctx
		vecB_Ori = cell2vec(cellOPI_Bl6Ctx);
		vecB_Dir = cell2vec(cellDSI_Bl6Ctx);
		vecA_Ori = cell2vec(cellOPI_AlbCtx);
		vecA_Dir = cell2vec(cellDSI_AlbCtx);
		
		strArea = cellAreaGroupsAbbr{intArea};
	else
		%not
		vecB_Ori = cell2vec(cellOPI_Bl6NOT);
		vecB_Dir = cell2vec(cellDSI_Bl6NOT);
		vecA_Ori = cell2vec(cellOPI_AlbNOT);
		vecA_Dir = cell2vec(cellDSI_AlbNOT);
		
		strArea = cellAreaGroupsAbbr{intArea};
	end
	
	%rem nans
	indRem = isnan(vecB_Ori) | isnan(vecB_Dir);
	vecB_Ori(indRem)=[];
	vecB_Dir(indRem)=[];
	indRem = isnan(vecA_Ori) | isnan(vecA_Dir);
	vecA_Ori(indRem)=[];
	vecA_Dir(indRem)=[];
	
	subplot(2,3,1+(intArea-1)*3)
	%bl6
	scatter(vecB_Ori,vecB_Dir);
	[h,pNOTBDir]=ttest(vecB_Dir);
	[h,pNOTBOri]=ttest(vecB_Ori);
	xlim([0 1]);
	ylim([-1 1]);
	title(sprintf('WT %s; t-test vs 0, DSI=%.3e,OSI=%.3e',strArea,pNOTBDir,pNOTBOri));
	xlabel('OSI (1-CV)');
	ylabel('L-R DSI');
	fixfig;
	grid off;
	
	subplot(2,3,2+(intArea-1)*3)
	%alb
	scatter(vecA_Ori,vecA_Dir,[],[1 0 0])
	[h,pNOTADir]=ttest(vecA_Dir);
	[h,pNOTAOri]=ttest(vecA_Ori);
	xlim([0 1]);
	ylim([-1 1]);
	title(sprintf('Alb %s; t-test vs 0, DSI=%.3f,OSI=%.3e',strArea,pNOTADir,pNOTAOri));
	xlabel('OSI (1-CV)');
	ylabel('L-R DSI');
	fixfig;
	grid off;
	
	subplot(2,3,3+(intArea-1)*3)
	%NOT bl6/alb
	errorbar(1:2,[mean(vecB_Ori) mean(vecB_Dir)],[std(vecB_Ori)/sqrt(numel(vecB_Ori)) std(vecB_Dir)/sqrt(numel(vecB_Dir))],'color',lines(1));
	hold on
	errorbar(1:2,[mean(vecA_Ori) mean(vecA_Dir)],[std(vecA_Ori)/sqrt(numel(vecA_Ori)) std(vecA_Dir)/sqrt(numel(vecA_Dir))],'color',[1 0 0]);
	fixfig;grid off
	ylabel(sprintf('Tuning strength %s',strArea));
	set(gca,'xtick',[1 2],'xticklabel',{'OSI (1-CV)','L-R DSI'});
	
	%anova
	vecB_Ori = vecB_Ori(:);
	vecA_Ori = vecA_Ori(:);
	vecB_Dir = vecB_Dir(:);
	vecA_Dir = vecA_Dir(:);
	y=cat(1,vecB_Ori,vecA_Ori,vecB_Dir,vecA_Dir);
	cellAlbG = cat(1,cellfill('WT',size(vecB_Ori)),cellfill('Alb',size(vecA_Ori)),cellfill('WT',size(vecB_Dir)),cellfill('Alb',size(vecA_Dir)));
	cellTuningG = cat(1,cellfill('Ori',size(vecB_Ori)),cellfill('Ori',size(vecA_Ori)),cellfill('Dir',size(vecB_Dir)),cellfill('Dir',size(vecA_Dir)));
	g = {cellAlbG,cellTuningG};
	[vecAnovaP,tbl,stats,terms]=anovan(y,g,'display','off','model','interaction');
	[h,dDSI]=ttest2(vecA_Dir,vecB_Dir);
	[h,dOSI]=ttest2(vecA_Ori,vecB_Ori);
	
	plot([1 2],[0 0],'k--')
	legend({'WT','Alb'},'location','best');
	title(sprintf('Inter, p=%.3e, t-test Alb/WT, OSI=%.3f,DSI=%.3e',vecAnovaP(3),dOSI,dDSI));
	fixfig;
	grid off;
end

% save figure
drawnow;
strFigName = ['OriVsDirTuning'];
export_fig(fullpath(strTargetPath,[strFigName '.tif']));
export_fig(fullpath(strTargetPath,[strFigName '.pdf']));
%% mean per rec
%{
figure;maxfig
subplot(2,3,1)
%ctx bl6
%vecCtxB_Ori = vecOSI_Bl6Ctx;
vecCtxB_Ori = vecOPI_Bl6Ctx;
vecCtxB_Dir = vecDSI_Bl6Ctx;
indRem = isnan(vecCtxB_Ori) | isnan(vecCtxB_Dir);
vecCtxB_Ori(indRem)=[];
vecCtxB_Dir(indRem)=[];
scatter(vecCtxB_Ori,vecCtxB_Dir)

subplot(2,3,2)
%ctx alb
%vecCtxA_Ori = vecOSI_AlbCtx;
vecCtxA_Ori = vecOPI_AlbCtx;
vecCtxA_Dir = vecDSI_AlbCtx;
indRem = isnan(vecCtxA_Ori) | isnan(vecCtxA_Dir);
vecCtxA_Ori(indRem)=[];
vecCtxA_Dir(indRem)=[];
scatter(vecCtxA_Ori,vecCtxA_Dir)

subplot(2,3,3)
%ctx bl6/alb
errorbar(1:2,[mean(vecCtxB_Ori) mean(vecCtxA_Ori)],[std(vecCtxB_Ori)/sqrt(numel(vecCtxB_Ori)) std(vecCtxA_Ori)/sqrt(numel(vecCtxA_Ori))]);
hold on
errorbar(1:2,[mean(vecCtxB_Dir) mean(vecCtxA_Dir)],[std(vecCtxB_Dir)/sqrt(numel(vecCtxB_Dir)) std(vecCtxA_Dir)/sqrt(numel(vecCtxA_Dir))]);
hold off
fixfig;grid off
xlim([0.8 2.2]);
ylabel('Tuning strength Ctx')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Alb'});
legend({'OSI','DSI'});

subplot(2,3,4)
%NOT bl6
%vecNOTB_Ori = vecOSI_Bl6NOT;
vecNOTB_Ori = vecOPI_Bl6NOT;
vecNOTB_Dir = vecDSI_Bl6NOT;
indRem = isnan(vecNOTB_Ori) | isnan(vecNOTB_Dir);
vecNOTB_Ori(indRem)=[];
vecNOTB_Dir(indRem)=[];
scatter(vecNOTB_Ori,vecNOTB_Dir)
title(sprintf('WT NOT; t-test vs 0, DSI=%.3f,OSI=%.3f'));
xlabel('OSI (1-CV)');
ylabel('L-R DSI');

subplot(2,3,5)
%NOT alb
%vecNOTA_Ori = vecOSI_AlbNOT;
vecNOTA_Ori = vecOPI_AlbNOT;
vecNOTA_Dir = vecDSI_AlbNOT;
indRem = isnan(vecNOTA_Ori) | isnan(vecNOTA_Dir);
vecNOTA_Ori(indRem)=[];
vecNOTA_Dir(indRem)=[];
scatter(vecNOTA_Ori,vecNOTA_Dir)
title('Alb NOT');
xlabel('OSI (1-CV)');
ylabel('L-R DSI');

subplot(2,3,6)
%NOT bl6/alb
errorbar(1:2,[mean(vecNOTB_Ori) mean(vecNOTA_Ori)],[std(vecNOTB_Ori)/sqrt(numel(vecNOTB_Ori)) std(vecNOTA_Ori)/sqrt(numel(vecNOTA_Ori))]);
hold on
errorbar(1:2,[mean(vecNOTB_Dir) mean(vecNOTA_Dir)],[std(vecNOTB_Dir)/sqrt(numel(vecNOTB_Dir)) std(vecNOTA_Dir)/sqrt(numel(vecNOTA_Dir))]);
hold off
fixfig;grid off
xlim([0.8 2.2]);
ylabel('Tuning strength NOT')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Alb'});
legend({'OSI (1-CV)','L-R DSI'});

%anova
vecNOTB_Ori = vecNOTB_Ori(:);
vecNOTA_Ori = vecNOTA_Ori(:);
vecNOTB_Dir = vecNOTB_Dir(:);
vecNOTA_Dir = vecNOTA_Dir(:);
y=cat(1,vecNOTB_Ori,vecNOTA_Ori,vecNOTB_Dir,vecNOTA_Dir);
cellAlbG = cat(1,cellfill('WT',size(vecNOTB_Ori)),cellfill('Alb',size(vecNOTA_Ori)),cellfill('WT',size(vecNOTB_Dir)),cellfill('Alb',size(vecNOTA_Dir)));
cellTuningG = cat(1,cellfill('Ori',size(vecNOTB_Ori)),cellfill('Ori',size(vecNOTA_Ori)),cellfill('Dir',size(vecNOTB_Dir)),cellfill('Dir',size(vecNOTA_Dir)));
g = {cellAlbG,cellTuningG};
[p,tbl,stats,terms]=anovan(y,g,'model','interaction');
%}