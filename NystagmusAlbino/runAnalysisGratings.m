%% further analyses
%{
[done/1) how do confusion matrices of ori decoding differ between alb/bl6?

[done/2) is neural code of nat movs more variable during eye movements in NOT than Ctx?

[done/3) does info in NOT predict info in V1? [is this doable? sufficient twin recording??]

[done/4) plot eye-tracking of stimulus 

%to do:
[done/A) plot tuning in NOT as function of location in NOT: is a recording closer to the border more
likely to be tuned? 

[done/B) plot tuning as second-closest pair decoding: is there a difference in decoding between horizontal
and vertical gratings? 

[done/C) plot results as separate recordings (and make selection of recordings based on visual
responsiveness)

[done/D) plot example NOT cells that respond visually but are not orientation tuned

[done/E) remake figure for tuning distribution over areas

[done/F) test whether eye-movement direction bias to temporonasal for DBAs is different from 0.5 and
different from BL6
=> is rightward bias of albinos significant? 
=> can this be explained by stimulation of only the left eye, causing temporonasal (rightward) motion?
==> was actually an error in counting nans

[done/G) spike shape NOT vs Ctx

H) make model explaining results

%to do 2:
[done/i) color waveform plot with firing rate
[done/ii) spike-triggered average of eye movement in NOT
[done/iii) decoding left vs right: is DBA worse than BL6?
[done/iv) plot x position as function of time of grating
v) make ccgs in NOT to check for local connections
vi) channels fast PV NOT?

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
%cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
%cellUseAreas{3} = {'Hippocampal','Field CA1','Field CA2','Field CA3','subiculum','Postsubiculum','Prosubiculum','dentate gyrus'};
cellUseAreas{3} = {'superior colliculus'};
cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
cellSubjectGroups = {'BL6','DBA'};

%% pre-allocate
cellTuningP_AlbCtx = {};
cellZetaP_AlbCtx = {};
cellPrefOri_AlbCtx = {};

cellTuningP_Bl6Ctx = {};
cellZetaP_Bl6Ctx = {};
cellPrefOri_Bl6Ctx = {};


cellTuningP_Bl6NOT = {};
cellPrefOri_Bl6NOT = {};
cellZetaP_Bl6NOT = {};

cellTuningP_AlbNOT = {};
cellZetaP_AlbNOT = {};
cellPrefOri_AlbNOT = {};


cellTuningP_Bl6Hip = {};
cellZetaP_Bl6Hip = {};
cellPrefOri_Bl6Hip = {};

cellTuningP_AlbHip = {};
cellZetaP_AlbHip = {};
cellPrefOri_AlbHip = {};


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
			numel({sTrialObjects.Orientation})
			%get data matrix
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			dblStimDur = median(vecStimOffTime-vecStimOnTime);
			matData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblStimDur);
			
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
				
				%% calc tuning curves & zeta
				%calc tuning curve
				matUseData = matData(vecSelectCells,:);
				sOut = getTuningCurves(matUseData,vecOrientation);
				vecTuningP_A = sOut.vecOriAnova;
				vecTuningP_R2 = sOut.vecFitP;
				%[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
				vecTuningP_R2_corr = vecTuningP_R2;
				vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
				indTunedCells = vecTuningP_R2<0.05;
				
				%calc zeta
				vecZetaP = nan(1,numel(vecSelectCells));
				dblTrialDur = 1.5;
				for intCellIdx=1:numel(vecSelectCells)
					intCell = vecSelectCells(intCellIdx);
					vecSpikeT = cellSpikeT{intCell};
					dblZetaP = getZeta(vecSpikeT,vecStimOnTime,dblTrialDur);
					vecZetaP(intCellIdx) = dblZetaP;
				end
				%[h, crit_p, vecZetaP_corr] = fdr_bh(vecZetaP);
				vecZetaP_corr = vecZetaP;
				if 0
					%% plot
					for intCellIdx=1:numel(vecSelectCells)
						intCell = vecSelectCells(intCellIdx);
						vecSpikeT = cellSpikeT{intCell};
						
						[vecTime,vecRate] = getIFR(vecSpikeT,vecStimOnTime,dblTrialDur);
						
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
						cellZetaP_AlbCtx{end+1} = vecZetaP_corr;
						cellPrefOri_AlbCtx{end+1} = vecPrefOri;
					elseif strcmpi(strSubjectType,'Bl6')
						cellTuningP_Bl6Ctx{end+1} = vecTuningP_R2_corr;
						cellZetaP_Bl6Ctx{end+1} = vecZetaP_corr;
						cellPrefOri_Bl6Ctx{end+1} = vecPrefOri;
					end
				elseif intArea == 2
					if strcmpi(strSubjectType,'DBA')
						cellTuningP_AlbNOT{end+1} = vecTuningP_R2_corr;
						cellZetaP_AlbNOT{end+1} = vecZetaP_corr;
						cellPrefOri_AlbNOT{end+1} = vecPrefOri;
					elseif strcmpi(strSubjectType,'Bl6')
						cellTuningP_Bl6NOT{end+1} = vecTuningP_R2_corr;
						cellZetaP_Bl6NOT{end+1} = vecZetaP_corr;
						cellPrefOri_Bl6NOT{end+1} = vecPrefOri;
					end
				elseif intArea == 3
					if strcmpi(strSubjectType,'DBA')
						cellTuningP_AlbHip{end+1} = vecTuningP_R2_corr;
						cellZetaP_AlbHip{end+1} = vecZetaP_corr;
					elseif strcmpi(strSubjectType,'Bl6')
						cellTuningP_Bl6Hip{end+1} = vecTuningP_R2_corr;
						cellZetaP_Bl6Hip{end+1} = vecZetaP_corr;
					end
				end
			end
			
		end
	end
end

%get all areas
cellAllAreas = arrayfun(@(x) unique({x.sCluster.Area}), sExp,'uniformoutput',false);
cellAllAreas = unique(horzcat(cellAllAreas{:}));

%% plot scatter of zeta-p versus tuning-p in Ctx, NOT, and Hip; colored by Alb/BL6
%remove nans
cellMouseType = {'Alb','Bl6'};
for intMouseType=1:2
	strMouseType = cellMouseType{intMouseType};
	for intAreaType=1:2
		strAreaType = cellAreaGroupsAbbr{intAreaType};
		
		strTuningVar = strcat('cellTuningP_',strMouseType,strAreaType);
		strZetaVar = strcat('cellZetaP_',strMouseType,strAreaType);
		strPrefVar = strcat('cellPrefOri_',strMouseType,strAreaType);
		cellTempTuningVar = eval(strTuningVar);
		cellTempZetaVar = eval(strZetaVar);
		cellTempPrefVar = eval(strPrefVar);
		
		intRecs = numel(cellTempTuningVar);
		for intRec=1:intRecs
			indNan = isnan(cellTempTuningVar{intRec}(:)) | isnan(cellTempZetaVar{intRec}(:)) | isnan(cellTempPrefVar{intRec}(:));
			cellTempTuningVar{intRec}(indNan) = [];
			cellTempZetaVar{intRec}(indNan) = [];
			cellTempPrefVar{intRec}(indNan) = [];
		end
		
		%remove recs with low #
		intCutOff = 3;
		indRemRecs = cellfun(@(x) numel(x)<intCutOff,cellTempTuningVar) | cellfun(@(x) numel(x)<intCutOff,cellTempZetaVar);
		cellTempTuningVar(indRemRecs) = [];
		cellTempZetaVar(indRemRecs) = [];
		cellTempPrefVar(indRemRecs) = [];
		
		
		eval([strTuningVar ' = cellTempTuningVar;'])
		eval([strZetaVar ' = cellTempZetaVar;'])
		eval([strPrefVar ' = cellTempPrefVar;'])
	end
end
%% vectorize
vecTuningP_AlbCtx = cell2vec(cellTuningP_AlbCtx);
vecZetaP_AlbCtx = cell2vec(cellZetaP_AlbCtx);
vecPrefOri_AlbCtx = cell2vec(cellPrefOri_AlbCtx);

vecTuningP_Bl6Ctx = cell2vec(cellTuningP_Bl6Ctx);
vecZetaP_Bl6Ctx = cell2vec(cellZetaP_Bl6Ctx);
vecPrefOri_Bl6Ctx = cell2vec(cellPrefOri_AlbCtx);

vecTuningP_AlbNOT = cell2vec(cellTuningP_AlbNOT);
vecZetaP_AlbNOT = cell2vec(cellZetaP_AlbNOT);
vecPrefOri_AlbNOT = cell2vec(cellPrefOri_AlbCtx);

vecTuningP_Bl6NOT = cell2vec(cellTuningP_Bl6NOT);
vecZetaP_Bl6NOT = cell2vec(cellZetaP_Bl6NOT);
vecPrefOri_Bl6NOT = cell2vec(cellPrefOri_AlbCtx);

%vecTuningP_AlbHip = cell2vec(cellTuningP_AlbHip);
%vecTuningP_Bl6Hip = cell2vec(cellTuningP_Bl6Hip);
%vecZetaP_AlbHip = cell2vec(cellZetaP_AlbHip);
%vecZetaP_Bl6Hip = cell2vec(cellZetaP_Bl6Hip);

%correlations
vecCorrAlbCtx = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbCtx,cellTuningP_AlbCtx);
vecCorrBl6Ctx = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6Ctx,cellTuningP_Bl6Ctx);
vecCorrAlbNOT = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbNOT,cellTuningP_AlbNOT);
vecCorrBl6NOT = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6NOT,cellTuningP_Bl6NOT);
%vecCorrAlbHip = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbHip,cellTuningP_AlbHip);
%vecCorrBl6Hip = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6Hip,cellTuningP_Bl6Hip);

%test
[h,p_Ctx] = ttest2(vecCorrAlbCtx,vecCorrBl6Ctx);
[h,p_Not] = ttest2(vecCorrAlbNOT,vecCorrBl6NOT);
%[h,p_Hip] = ttest2(vecCorrAlbHip,vecCorrBl6Hip);

%% plot
figure;maxfig;
vecRawP = nan(1,6);
for intAreaType=1:2
	strAreaType = cellAreaGroupsAbbr{intAreaType};
	
	strTuningAlbVar = strcat('vecTuningP_Alb',strAreaType);
	strZetaAlbVar = strcat('vecZetaP_Alb',strAreaType);
	strTuningBl6Var = strcat('vecTuningP_Bl6',strAreaType);
	strZetaBl6Var = strcat('vecZetaP_Bl6',strAreaType);
	
	vecTuningP_AlbTemp = eval(strTuningAlbVar);
	vecTuningP_Bl6Temp = eval(strTuningBl6Var);
	vecZetaP_AlbTemp = eval(strZetaAlbVar);
	vecZetaP_Bl6Temp = eval(strZetaBl6Var);
	
	vecTuningZ_Bl6 = -norminv(vecTuningP_Bl6Temp/2);
	vecTuningZ_Alb = -norminv(vecTuningP_AlbTemp/2);
	vecZetaZ_Bl6 = -norminv(vecZetaP_Bl6Temp/2);
	vecZetaZ_Alb = -norminv(vecZetaP_AlbTemp/2);
	
	vecColAlb = [0.9 0 0];
	vecColBl6 = lines(1);
	
	
	dblMax = max(cat(1,vecTuningZ_Bl6,vecTuningZ_Alb,vecZetaZ_Bl6,vecZetaZ_Alb));
	dblStep = 1;
	vecBins=0:dblStep:dblMax;
	vecPlotBins = vecBins(2:end)-dblStep/2;
	[h,p_ZZN] = ttest2(vecZetaZ_Bl6,vecZetaZ_Alb);
	vecRawP(intAreaType) = p_ZZN;
	vecCountsZetaZ_Bl6=histcounts(vecZetaZ_Bl6,vecBins);
	vecCountsZetaZ_Alb=histcounts(vecZetaZ_Alb,vecBins);
	subplot(2,3,intAreaType)
	plot(vecPlotBins,vecCountsZetaZ_Bl6./sum(vecCountsZetaZ_Bl6),'color',vecColBl6);
	hold on
	plot(vecPlotBins,vecCountsZetaZ_Alb./sum(vecCountsZetaZ_Alb),'color',vecColAlb);
	hold off
	legend({'BL6','DBA'},'location','best');
	xlabel('Responsiveness z-score (ZETA)');
	ylabel('Fraction of cells (norm. count)');
	fixfig;
	title(sprintf('%s, ZETA-Z DBA (mu=%.2f) BL6 (mu=%.2f),p=%.5f',strAreaType,mean(vecZetaZ_Alb),mean(vecZetaZ_Bl6),p_ZZN));
	
	[h,p_TZN] = ttest2(vecTuningZ_Bl6,vecTuningZ_Alb);
	vecRawP(intAreaType+3) = p_TZN;
	vecCountsTuningZ_Bl6=histcounts(vecTuningZ_Bl6,vecBins);
	vecCountsTuningZ_Alb=histcounts(vecTuningZ_Alb,vecBins);
	subplot(2,3,intAreaType+3)
	plot(vecPlotBins,vecCountsTuningZ_Bl6./sum(vecCountsTuningZ_Bl6),'color',vecColBl6);
	hold on
	plot(vecPlotBins,vecCountsTuningZ_Alb./sum(vecCountsTuningZ_Alb),'color',vecColAlb);
	hold off
	legend({'BL6','DBA'},'location','best');
	xlabel('Ori tuning z-score (sd)');
	ylabel('Fraction of cells (norm. count)');
	fixfig;
	title(sprintf('%s, Tuning-Z DBA (mu=%.2f) BL6 (mu=%.2f),p=%.5f',strAreaType,mean(vecTuningZ_Alb),mean(vecTuningZ_Bl6),p_TZN));
end
[h,i,vecCorr_P] = fdr_bh(vecRawP(~isnan(vecRawP)));
[vecCorr_P2] = bonf_holm(vecRawP);
vecCorr_P2(vecCorr_P2>1)=1;

%% is BL6/DBA difference larger in NOT than Ctx?
vecZetaZ_AlbCtx = -norminv(cell2vec(cellZetaP_AlbCtx)/2);
vecZetaZ_Bl6Ctx = -norminv(cell2vec(cellZetaP_Bl6Ctx)/2);
vecZetaZ_AlbNOT = -norminv(cell2vec(cellZetaP_AlbNOT)/2);
vecZetaZ_Bl6NOT = -norminv(cell2vec(cellZetaP_Bl6NOT)/2);

vecG1 = cat(1,ones(size(vecZetaZ_Bl6Ctx)),2*ones(size(vecZetaZ_AlbCtx)),ones(size(vecZetaZ_Bl6NOT)),2*ones(size(vecZetaZ_AlbNOT))); %bl6 vs alb
cellG1 = cellSubjectGroups(vecG1)';
vecG2 = cat(1,ones(size(vecZetaZ_Bl6Ctx)),ones(size(vecZetaZ_AlbCtx)),2*ones(size(vecZetaZ_Bl6NOT)),2*ones(size(vecZetaZ_AlbNOT))); %ctx vs not
cellG2 = cellAreaGroupsAbbr(vecG2)';

vecY2 = cat(1,vecZetaZ_Bl6Ctx,vecZetaZ_AlbCtx,vecZetaZ_Bl6NOT,vecZetaZ_AlbNOT);

[p2,tbl2,stats2,terms2] = anovan(vecY2,{cellG1,cellG2},'model','full','display','off');

%plot
vecZetaZ_Mu = [mean(vecZetaZ_Bl6Ctx) mean(vecZetaZ_AlbCtx) mean(vecZetaZ_Bl6NOT) mean(vecZetaZ_AlbNOT)]; 
vecZetaZ_SEM = [std(vecZetaZ_Bl6Ctx)/sqrt(numel(vecZetaZ_Bl6Ctx))...
	std(vecZetaZ_AlbCtx)/sqrt(numel(vecZetaZ_AlbCtx))...
	std(vecZetaZ_Bl6NOT)/sqrt(numel(vecZetaZ_Bl6NOT))...
	std(vecZetaZ_AlbNOT)/sqrt(numel(vecZetaZ_AlbNOT))]; 

subplot(2,3,3)
errorbar([1 2],[mean(vecZetaZ_Bl6Ctx) mean(vecZetaZ_Bl6NOT)],...
	[std(vecZetaZ_Bl6Ctx)/sqrt(numel(vecZetaZ_Bl6Ctx)) std(vecZetaZ_Bl6NOT)/sqrt(numel(vecZetaZ_Bl6Ctx))],...
	'x-','Color',vecColBl6,'CapSize',20);
hold on
errorbar([1 2],[mean(vecZetaZ_AlbCtx) mean(vecZetaZ_AlbNOT)],...
	[std(vecZetaZ_AlbCtx)/sqrt(numel(vecZetaZ_AlbCtx)) std(vecZetaZ_AlbNOT)/sqrt(numel(vecZetaZ_AlbCtx))],...
	'x-','Color',vecColAlb,'CapSize',20);
hold off
ylabel('Mean responsiveness (ZETA)');
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
legend(cellSubjectGroups,'Location','best');
title(sprintf('Interaction Ctx/NOT - Alb/BL6,p=%.6f',p2(3)));
fixfig;grid off;

%with cohens d
vecTuningZ_AlbCtx = -norminv(cell2vec(cellTuningP_AlbCtx)/2);
vecTuningZ_Bl6Ctx = -norminv(cell2vec(cellTuningP_Bl6Ctx)/2);
vecTuningZ_AlbNOT = -norminv(cell2vec(cellTuningP_AlbNOT)/2);
vecTuningZ_Bl6NOT = -norminv(cell2vec(cellTuningP_Bl6NOT)/2);

vecY = cat(1,vecTuningZ_Bl6Ctx,vecTuningZ_AlbCtx,vecTuningZ_Bl6NOT,vecTuningZ_AlbNOT);
vecG1 = cat(1,ones(size(vecTuningZ_Bl6Ctx)),2*ones(size(vecTuningZ_AlbCtx)),ones(size(vecTuningZ_Bl6NOT)),2*ones(size(vecTuningZ_AlbNOT))); %bl6 vs alb
cellG1 = cellSubjectGroups(vecG1)';
vecG2 = cat(1,ones(size(vecTuningZ_Bl6Ctx)),ones(size(vecTuningZ_AlbCtx)),2*ones(size(vecTuningZ_Bl6NOT)),2*ones(size(vecTuningZ_AlbNOT))); %ctx vs not
cellG2 = cellAreaGroupsAbbr(vecG2)';

[p,tbl,stats,terms] = anovan(vecY,{cellG1,cellG2},'model','full','display','off');

subplot(2,3,6)
errorbar([1 2],[mean(vecTuningZ_Bl6Ctx) mean(vecTuningZ_Bl6NOT)],...
	[std(vecTuningZ_Bl6Ctx)/sqrt(numel(vecTuningZ_Bl6Ctx)) std(vecTuningZ_Bl6NOT)/sqrt(numel(vecTuningZ_Bl6Ctx))],...
	'x-','Color',vecColBl6,'CapSize',20);
hold on
errorbar([1 2],[mean(vecTuningZ_AlbCtx) mean(vecTuningZ_AlbNOT)],...
	[std(vecTuningZ_AlbCtx)/sqrt(numel(vecTuningZ_AlbCtx)) std(vecTuningZ_AlbNOT)/sqrt(numel(vecTuningZ_AlbCtx))],...
	'x-','Color',vecColAlb,'CapSize',20);
hold off
ylabel('Mean tuning (z-score)');
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
legend(cellSubjectGroups,'Location','best');
title(sprintf('Interaction Ctx/NOT - Alb/BL6,p=%.6f',p(3)));
fixfig;grid off;drawnow;

export_fig(fullpath(strTargetPath,'TuningAndResponsivenessCtxNOT.tif'));
export_fig(fullpath(strTargetPath,'TuningAndResponsivenessCtxNOT.pdf'));
%% plot
%set scatter vals
dblSize = 200;

figure;maxfig;
for intAreaGroup=1:2
	vecTempAlb = sort(eval(['vecCorrAlb' cellAreaGroupsAbbr{intAreaGroup} ';']));
	vecTempBl6 = sort(eval(['vecCorrBl6' cellAreaGroupsAbbr{intAreaGroup} ';']));
	[h,p_corr] = ttest2(vecTempAlb,vecTempBl6);
	subplot(2,3,intAreaGroup)
	h=errorbar([1 2],[mean(vecTempBl6) mean(vecTempAlb)],[std(vecTempBl6)./sqrt(numel(vecTempBl6)) std(vecTempAlb)./sqrt(numel(vecTempAlb))],...
		'color',0.2*[1 1 1],'Marker','x','MarkerSize',10,'linestyle','none','capsize',50);
	hold on
	%bl6
	vecJitter = (rand(size(vecTempBl6))-0.5)*2*(1/16) + (mod(1:numel(vecTempBl6),2)*2-1)*(1/16);
	scatter(ones(size(vecTempBl6)) + vecJitter,vecTempBl6,dblSize,lines(1),'.');
	%alb
	vecJitter = (rand(size(vecTempAlb))-0.5)*2*(1/16) + (mod(1:numel(vecTempAlb),2)*2-1)*(1/16);
	scatter(2*ones(size(vecTempAlb)) + vecJitter,vecTempAlb,dblSize,[0.9 0 0],'.');
	title(sprintf('%s; r_BL6=%.3f (b);r_DBA=%.3f (r);p=%.6f',cellAreaGroupsAbbr{intAreaGroup},mean(vecTempBl6),mean(vecTempAlb),p_corr),'interpreter','none');
	ylabel('Correlation R(Tuning, Responsive)');
	set(gca,'xtick',[1 2],'xticklabel',{'BL6','DBA'});
	hold off
	xlim([0.5 2.5])
	ylim([-1 1]);
	fixfig;grid off
	h.LineWidth=3;
	
	%plot examples
	%subplot(2,3,intAreaGroup+3)
	vecExampleAlbZ = eval(['cellZetaP_Alb' cellAreaGroupsAbbr{intAreaGroup} '{1}(:);']);
	vecExampleAlbT = eval(['cellTuningP_Alb' cellAreaGroupsAbbr{intAreaGroup} '{1}(:);']);
	
	vecExampleBl6Z = eval(['cellZetaP_Bl6' cellAreaGroupsAbbr{intAreaGroup} '{1}(:);']);
	vecExampleBl6T = eval(['cellTuningP_Bl6' cellAreaGroupsAbbr{intAreaGroup} '{1}(:);']);
	
	%scatter(-norminv(vecExampleBl6Z/2),-norminv(vecExampleBl6T/2),[],vecColBl6)
	%hold on
	%scatter(-norminv(vecExampleAlbZ/2),-norminv(vecExampleAlbT/2),[],vecColAlb)
	%hold off
end


%% plot scatter of zeta-p, tuning-p per recording in Ctx, NOT, and Hip; colored by Alb/BL6
vecZBC = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_Bl6Ctx);
vecTBC = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_Bl6Ctx);
vecZBN = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_Bl6NOT);
vecTBN = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_Bl6NOT);
%vecZBH = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_Bl6Hip);
%vecTBH = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_Bl6Hip);

vecZAC = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_AlbCtx);
vecTAC = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_AlbCtx);
vecZAN = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_AlbNOT);
vecTAN = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_AlbNOT);
%vecZAH = cellfun(@(x) mean(-norminv(x/2)),cellZetaP_AlbHip);
%vecTAH = cellfun(@(x) mean(-norminv(x/2)),cellTuningP_AlbHip);

[h,p2a]=ttest2(vecZBC,vecZAC);
[h,p2b]=ttest2(vecZBN,vecZAN);
%[h,p2c]=ttest2(vecZBH,vecZAH);

[h,p2a]=ttest2(vecTBC,vecTAC);
[h,p2b]=ttest2(vecTBN,vecTAN);
%[h,p2c]=ttest2(vecTBH,vecTAH);

%% plot corr(zeta-p,tuning-p) per recording in Ctx, NOT, and Hip; colored by Alb/BL6
return
%ctx
dblR_AlbCtx = nancorr(-norminv(vecZetaP_AlbCtx/2),-norminv(vecTuningP_AlbCtx/2));
dblR_Bl6Ctx = nancorr(-norminv(vecZetaP_Bl6Ctx/2),-norminv(vecTuningP_Bl6Ctx/2));
figure
subplot(2,3,1)
hold on;
scatter(-norminv(vecZetaP_AlbCtx/2),-norminv(vecTuningP_AlbCtx/2),[],'r','x');
scatter(-norminv(vecZetaP_Bl6Ctx/2),-norminv(vecTuningP_Bl6Ctx/2),[],'b','x');
hold off
fixfig;

%not
[dblR_AlbNot,p_an,rl_an,ru_an] = corrcoef(-norminv(vecZetaP_AlbNOT/2),-norminv(vecTuningP_AlbNOT/2));
[dblR_Bl6Not,p_bn,rl_bn,ru_bn] = corrcoef(-norminv(vecZetaP_Bl6NOT/2),-norminv(vecTuningP_Bl6NOT/2));
x1=-norminv(vecZetaP_AlbNOT/2);
y1=-norminv(vecTuningP_AlbNOT/2);
x2=-norminv(vecZetaP_Bl6NOT/2);
y2=-norminv(vecTuningP_Bl6NOT/2);
[p,z,r1,r2,S] = corrtest2(x1,y1,x2,y2);
p

subplot(2,3,2)
hold on;
scatter(-norminv(vecZetaP_AlbNOT/2),-norminv(vecTuningP_AlbNOT/2),[],'r','x');
scatter(-norminv(vecZetaP_Bl6NOT/2),-norminv(vecTuningP_Bl6NOT/2),[],'b','x');
hold off
fixfig;

%hip
%dblR_AlbHip = nancorr(-norminv(vecZetaP_AlbHip/2),-norminv(vecTuningP_AlbHip/2));
%dblR_Bl6Hip = nancorr(-norminv(vecZetaP_Bl6Hip/2),-norminv(vecTuningP_Bl6Hip/2));

%subplot(2,3,3)
%hold on;
%scatter(-norminv(vecZetaP_AlbHip/2),-norminv(vecTuningP_AlbHip/2),[],'r','x');
%scatter(-norminv(vecZetaP_Bl6Hip/2),-norminv(vecTuningP_Bl6Hip/2),[],'b','x');
%hold off
%fixfig;