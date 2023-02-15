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
	vecRunRecs = find(indUseRecs & ~(indRemRecs));% | indRemRecs2));
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
		vecAggDir = [];
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
			if isfield(sRec,'sPupil') && ~isempty(sRec.sPupil)
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
			else
				%remove trials with blinking
				indRemTrials = false(size(vecPupilStimOn));
			end
			%% prep data
			%split by ori
			sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			vecOrientation(indRemTrials) = [];
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24 || numel( sBlock.vecStimOnTime) ~= numel(sTrialObjects),continue,end
			
			%get data matrix
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			
			vecAggDir = cat(1,vecAggDir,vecOrientation);
			vecAggStimOnTime = cat(2,vecAggStimOnTime,vecStimOnTime);
			vecAggStimOffTime = cat(2,vecAggStimOffTime,vecStimOffTime);
		end
		%get data matrix
		cellSpikeT = {sRec.sCluster(:).SpikeTimes};
		dblStimDur = median(vecAggStimOffTime-vecAggStimOnTime);
		[vecAggOriIdx,vecAggUnique,vecAggCounts,cellAggSelect,vecAggRepetition] = val2idx(vecAggDir);
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
		vecOriNoDir = mod(vecAggDir,180);
		vecTrialTypesNoDir = deg2rad(vecOriNoDir)*2;
		[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTrialTypesNoDir);
		dblLambda = 100;
		
		%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
		[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecAggDir);
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
			cellSubSpikes = cellSpikeT(vecSelectCells);
			sOut = getTuningCurves(matUseData,vecAggDir);
			vecTuningP_A = sOut.vecOriAnova;
			vecTuningP_R2 = sOut.vecFitP;
			[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
			%vecTuningP_R2_corr = vecTuningP_R2;
			vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
			indTunedCells = vecTuningP_R2_corr<0.05;
			
			%get OPI
			matUseR = sOut.matMeanResp;%(indTunedCells,:);
			%matUseR = sOut.matFittedResp;
			vecDir24 = sOut.vecUniqueDegs;
			vecOri24 = mod(vecDir24,180);
			
			%collapse angles
			vecOPI = getOPI(matUseR,deg2rad(vecOri24));
			vecOSI = getOSI(matUseR,deg2rad(vecOri24));
			
			%get DSI
			%intDirNum = numel(unique(vecDir24));
			%[vecPrefR,vecPrefIdx]=max(matUseR,[],2);
			%vecAntiIdx = modx(vecPrefIdx+intDirNum/2,intDirNum);
			%vecAntiR = getMatVals(matUseR,(1:numel(vecAntiIdx))',vecAntiIdx);
			%vecDSI = (vecPrefR-vecAntiR)./(vecAntiR+vecPrefR);
			
			%get DSI horizontal
			matUseRD = matUseR(indTunedCells,:);
			vecLeftR = mean(matUseRD(:,vecDir24==0 | vecDir24==345 | vecDir24==15),2);
			vecRightR = mean(matUseRD(:,vecDir24==165 | vecDir24==180 | vecDir24==195),2);
			vecDSI = (vecLeftR-vecRightR)./(vecRightR+vecLeftR);
			
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
				%% plot tuning curve
				[dummy,vecIdx]=sort(vecTuningP_R2_corr);
				%get DSI horizontal
				matUseRD = matUseR;
				vecLeftR = mean(matUseRD(:,vecDir24==0 | vecDir24==345 | vecDir24==15),2);
				vecRightR = mean(matUseRD(:,vecDir24==165 | vecDir24==180 | vecDir24==195),2);
				vecDSI = (vecLeftR-vecRightR)./(vecRightR+vecLeftR);
				
				%% plot
				for intCellIdx=1:numel(vecIdx)
					%%
					%intCellIdx=6
					intCell=vecIdx(intCellIdx);
					%intCell=13;
					figure;maxfig;
					vecDirs = vecAggUnique(:)';
					vecMeanR = sOut.matMeanResp(intCell,:);
					vecMeanSem = sOut.matSDResp(intCell,:)./sqrt(vecAggCounts');
					vecFitR = sOut.matFittedResp(intCell,:);
					intPlotPre = 12;
					intPlotPost = -12;
					vecPlotX = [(numel(vecDirs)-intPlotPre+1):numel(vecDirs) 1:(numel(vecDirs)+intPlotPost)];
					vecDirs_Ext = vecDirs(vecPlotX);
					vecDirs_Ext(1:intPlotPre) = vecDirs_Ext(1:intPlotPre)-360;
					vecMeanR_Ext = vecMeanR(vecPlotX);
					vecMeanSem_Ext = vecMeanSem(vecPlotX);
					vecFitR_Ext = vecFitR(vecPlotX);
					subplot(2,3,1)
					polarplot(deg2rad(vecDirs),vecMeanR);
					hold on
					polarplot(deg2rad(vecDirs),vecFitR);
					hold off
					fixfig;
					title(sprintf('OSI=%.3f, DSI=%.3f',vecOPI(intCell),vecDSI(intCell)));
					
					subplot(2,3,2)
					errorbar(vecDirs_Ext,vecMeanR_Ext,vecMeanSem_Ext);
					hold on
					plot(vecDirs_Ext,vecFitR_Ext)
					hold off
					xlabel('Stimulus direction (degs)');
					ylabel('Spiking rate (Hz)');
					title(sprintf('%d',intCell))
					fixfig;grid off;
					%xlim([0 360]);
					ylim([0 max(get(gca,'ylim'))]);
					
					subplot(2,3,3)
					[dummy,intPrefDir] = max(vecMeanR);
					vecPrefT =  vecAggStimOnTime(cellAggSelect{intPrefDir});
					dblOffset = -0.3;
					dblTrialDur = 1.5-dblOffset;
					dblBinSize = 0.05;
					vecWindow = dblOffset:dblBinSize:(dblTrialDur+dblOffset);
					vecLimX = [vecWindow(1) vecWindow(end)];
					vecSpikeT = cellSubSpikes{intCell};
					intSmoothSd = 5;
					%[vecTime,vecRate] = getIFR(vecSpikeT,vecPrefT+dblOffset,dblTrialDur,intSmoothSd);
					[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikeT,vecWindow,vecPrefT,-1);
					%plot(vecTime+dblOffset,vecRate);
					hold on
					errorbar(vecWindowBinCenters,vecMean,vecSEM);
					hold off
					
					xlabel('Time (s)');
					ylabel('Spiking rate (Hz)');
					title('pref dir')
					fixfig;grid off;
					xlim(vecLimX);
					
					subplot(2,3,6)
					intDirNum = numel(cellAggSelect);
					intOrthDir = modx(intPrefDir + intDirNum/4,intDirNum);
					vecOrthT =  vecAggStimOnTime(cellAggSelect{intOrthDir});
					[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikeT,vecWindow,vecOrthT,-1);
					%plot(vecTime,vecRate);
					hold on
					errorbar(vecWindowBinCenters,vecMean,vecSEM);
					hold off
					xlabel('Time (s)');
					ylabel('Instantaneous rate (Hz)');
					title('orth dir')
					fixfig;grid off;
					xlim(vecLimX);
					
					% save figure
					drawnow;
					strFigName = sprintf('ExampleTuningCurve%s_%s%d',strName,cellAreaGroupsAbbr{intArea},intCell);
					export_fig(fullpath([strTargetPath filesep 'single_cells'],[strFigName '.tif']));
					export_fig(fullpath([strTargetPath filesep 'single_cells'],[strFigName '.pdf']));
					%%
					pause
				end
				
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
figure;maxfig
for intMetric=1:2
	if intMetric == 1
		%% OSI
		strMetric = 'OSI (1-CV)';
		vecBinsE = 0:0.1:1;
		%ctx
		vecBl6_Ctx = cell2vec(cellOPI_Bl6Ctx);
		vecAlb_Ctx = cell2vec(cellOPI_AlbCtx);
		
		%not
		vecBl6_NOT = cell2vec(cellOPI_Bl6NOT);
		vecAlb_NOT = cell2vec(cellOPI_AlbNOT);
	else
		%% DSI
		vecBinsE = -1:0.2:1;
		strMetric = 'L-R DSI';
		vecBl6_Ctx = cell2vec(cellDSI_Bl6Ctx);
		vecAlb_Ctx = cell2vec(cellDSI_AlbCtx);
		vecBl6_NOT = cell2vec(cellDSI_Bl6NOT);
		vecAlb_NOT = cell2vec(cellDSI_AlbNOT);
	end
	vecBinsC = vecBinsE(2:end)-median(diff(vecBinsE))/2;
	%rem nans
	vecBl6_Ctx(isnan(vecBl6_Ctx))=[];
	vecAlb_Ctx(isnan(vecAlb_Ctx))=[];
	
	vecBl6_NOT(isnan(vecBl6_NOT))=[];
	vecAlb_NOT(isnan(vecAlb_NOT))=[];
	
	
	subplot(2,3,1+(intMetric-1)*3)
	%ctx
	vecBl6_Ctx_C = histcounts(vecBl6_Ctx,vecBinsE);
	vecAlb_Ctx_C = histcounts(vecAlb_Ctx,vecBinsE);
	hold on
	plot(vecBinsC,vecBl6_Ctx_C./sum(vecBl6_Ctx_C(:)),'color',lines(1));
	plot(vecBinsC,vecAlb_Ctx_C./sum(vecAlb_Ctx_C(:)),'color',[1 0 0]);
	hold off
	xlabel(strMetric);
	ylabel('Fraction of Ctx cells');
	fixfig;
	
	subplot(2,3,2+(intMetric-1)*3)
	%ctx
	vecBl6_NOT_C = histcounts(vecBl6_NOT,vecBinsE);
	vecAlb_NOT_C = histcounts(vecAlb_NOT,vecBinsE);
	hold on
	plot(vecBinsC,vecBl6_NOT_C./sum(vecBl6_NOT_C(:)),'color',lines(1));
	plot(vecBinsC,vecAlb_NOT_C./sum(vecAlb_NOT_C(:)),'color',[1 0 0]);
	hold off
	xlabel(strMetric);
	ylabel('Fraction of NOT cells');
	fixfig;
	
	subplot(2,3,3+(intMetric-1)*3)
	%WT
	errorbar(1:2,[mean(vecBl6_Ctx) mean(vecBl6_NOT)],[std(vecBl6_Ctx)/sqrt(numel(vecBl6_Ctx)) std(vecBl6_NOT)/sqrt(numel(vecBl6_NOT))],'color',lines(1));
	
	hold on
	%Alb
	errorbar(1:2,[mean(vecAlb_Ctx) mean(vecAlb_NOT)],[std(vecAlb_Ctx)/sqrt(numel(vecAlb_Ctx)) std(vecAlb_NOT)/sqrt(numel(vecAlb_NOT))],'color',[1 0 0]);
	%scatter(ones(size(vecBl6_Ctx))*0.8,vecBl6_Ctx,[],lines(1));
	%scatter(ones(size(vecAlb_Ctx))*0.9,vecAlb_Ctx,[],[1 0 0]);
	%scatter(ones(size(vecBl6_NOT))*2.1,vecBl6_NOT,[],lines(1));
	%scatter(ones(size(vecAlb_NOT))*2.2,vecAlb_NOT,[],[1 0 0]);
	
	%fixfig;grid off
	ylabel(strMetric);
	set(gca,'xtick',[1 2],'xticklabel',{'Ctx','NOT'});
	
	%anova
	vecBl6_Ctx = vecBl6_Ctx(:);
	vecAlb_Ctx = vecAlb_Ctx(:);
	vecBl6_NOT = vecBl6_NOT(:);
	vecAlb_NOT = vecAlb_NOT(:);
	y=cat(1,vecBl6_Ctx,vecAlb_Ctx,vecBl6_NOT,vecAlb_NOT);
	cellMouseG = cat(1,cellfill('WT',size(vecBl6_Ctx)),cellfill('Alb',size(vecAlb_Ctx)),cellfill('WT',size(vecBl6_NOT)),cellfill('Alb',size(vecAlb_NOT)));
	cellAreaG = cat(1,cellfill('Ctx',size(vecBl6_Ctx)),cellfill('Ctx',size(vecAlb_Ctx)),cellfill('NOT',size(vecBl6_NOT)),cellfill('NOT',size(vecAlb_NOT)));
	g = {cellMouseG,cellAreaG};
	[vecAnovaP,tbl,stats,terms]=anovan(y,g,'display','off','model','interaction');
	[h,dCtx]=ttest2(vecBl6_Ctx,vecAlb_Ctx);
	dCtxW=ranksum(vecBl6_Ctx,vecAlb_Ctx);
	[h,dNOT]=ttest2(vecBl6_NOT,vecAlb_NOT);
	dNOTW=ranksum(vecBl6_NOT,vecAlb_NOT);
	
	[h,x,p]=fdr_bh([dCtx dNOT]);
	dCtxC=p(1);
	dNOTC=p(2);
	
	plot([1 2],[0 0],'k--')
	legend({'WT','Alb'},'location','best');
	title(sprintf('Inter, p=%.4f, t-test OSI, Ctx=%.4f,NOT=%.4f',vecAnovaP(3),dCtxC,dNOTC));
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