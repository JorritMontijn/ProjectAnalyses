%% further analyses
%{
1) how do confusion matrices of ori decoding differ between alb/bl6? 

2) is neural code of nat movs more variable during eye movements in NOT than Ctx? 

3) does info in NOT predict info in V1? 

4) spike shape NOT vs Ctx 

%to do:
A) plot tuning in NOT as function of location in NOT: is a recording closer to the border more likely
to be tuned?

B) plot tuning as second-closest pair decoding

C) plot results as separate recordings (and make selection of recordings based on visual
responsiveness)

D) plot example NOT cells that respond visually but are not orientation tuned

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
cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
cellUseAreas{3} = {'Hippocampal formation','Field CA1','Field CA2','Field CA3','subiculum','dentate gyrus'};
cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};

%% pre-allocate
cellTuningP_AlbCtx = {};
cellTuningP_Bl6Ctx = {};
cellZetaP_AlbCtx = {};
cellZetaP_Bl6Ctx = {};

cellTuningP_AlbNot = {};
cellTuningP_Bl6Not = {};
cellZetaP_AlbNot = {};
cellZetaP_Bl6Not = {};

cellTuningP_AlbHip = {};
cellTuningP_Bl6Hip = {};
cellZetaP_AlbHip = {};
cellZetaP_Bl6Hip = {};
			
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
			
			%get data matrix
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			dblStimDur = median(vecStimOffTime-vecStimOnTime);
			matData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblStimDur);
			
			%include?
			vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
			indUseCells = arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
			
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
				
				%% calc tuning curves & zeta
				%calc tuning curve
				matUseData = matData(vecSelectCells,:);
				sOut = getTuningCurves(matUseData,vecOrientation);
				vecTuningP_A = sOut.vecOriAnova;
				vecTuningP_R2 = sOut.vecFitP;
				[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
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
				[h, crit_p, vecZetaP_corr] = fdr_bh(vecZetaP);
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
					elseif strcmpi(strSubjectType,'Bl6')
						cellTuningP_Bl6Ctx{end+1} = vecTuningP_R2_corr;
						cellZetaP_Bl6Ctx{end+1} = vecZetaP_corr;
					end
				elseif intArea == 2
					if strcmpi(strSubjectType,'DBA')
						cellTuningP_AlbNot{end+1} = vecTuningP_R2_corr;
						cellZetaP_AlbNot{end+1} = vecZetaP_corr;
					elseif strcmpi(strSubjectType,'Bl6')
						cellTuningP_Bl6Not{end+1} = vecTuningP_R2_corr;
						cellZetaP_Bl6Not{end+1} = vecZetaP_corr;
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


%% plot scatter of zeta-p versus tuning-p in Ctx, NOT, and Hip; colored by Alb/BL6
%vectorize
vecTuningP_AlbCtx = cell2vec(cellTuningP_AlbCtx);
vecTuningP_Bl6Ctx = cell2vec(cellTuningP_Bl6Ctx);
vecZetaP_AlbCtx = cell2vec(cellZetaP_AlbCtx);
vecZetaP_Bl6Ctx = cell2vec(cellZetaP_Bl6Ctx);

vecTuningP_AlbNot = cell2vec(cellTuningP_AlbNot);
vecTuningP_Bl6Not = cell2vec(cellTuningP_Bl6Not);
vecZetaP_AlbNot = cell2vec(cellZetaP_AlbNot);
vecZetaP_Bl6Not = cell2vec(cellZetaP_Bl6Not);

vecTuningP_AlbHip = cell2vec(cellTuningP_AlbHip);
vecTuningP_Bl6Hip = cell2vec(cellTuningP_Bl6Hip);
vecZetaP_AlbHip = cell2vec(cellZetaP_AlbHip);
vecZetaP_Bl6Hip = cell2vec(cellZetaP_Bl6Hip);

%correlations
vecCorrAlbCtx = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbCtx,cellTuningP_AlbCtx);
vecCorrBl6Ctx = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6Ctx,cellTuningP_Bl6Ctx);
vecCorrAlbNot = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbNot,cellTuningP_AlbNot);
vecCorrBl6Not = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6Not,cellTuningP_Bl6Not);
vecCorrAlbHip = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_AlbHip,cellTuningP_AlbHip);
vecCorrBl6Hip = cellfun(@(x,y) nancorr(x(:),y(:)),cellZetaP_Bl6Hip,cellTuningP_Bl6Hip);

%% plot scatter of corr(zeta-p,tuning-p) per recording in Ctx, NOT, and Hip; colored by Alb/BL6
error('remove nans & recordings with few cells, plot recording-based corrs, test significance')
%ctx
dblR_AlbCtx = nancorr(-norminv(vecZetaP_AlbCtx/2),-norminv(vecTuningP_AlbCtx/2));
dblR_Bl6Ctx = nancorr(-norminv(vecZetaP_Bl6Ctx/2),-norminv(vecTuningP_Bl6Ctx/2));

subplot(2,3,1)
hold on;
scatter(-norminv(vecZetaP_AlbCtx/2),-norminv(vecTuningP_AlbCtx/2),[],'r','x');
scatter(-norminv(vecZetaP_Bl6Ctx/2),-norminv(vecTuningP_Bl6Ctx/2),[],'b','x');
hold off
fixfig;

%not
[dblR_AlbNot,p_an,rl_an,ru_an] = corrcoef(-norminv(vecZetaP_AlbNot/2),-norminv(vecTuningP_AlbNot/2));
[dblR_Bl6Not,p_bn,rl_bn,ru_bn] = corrcoef(-norminv(vecZetaP_Bl6Not/2),-norminv(vecTuningP_Bl6Not/2));

subplot(2,3,2)
hold on;
scatter(-norminv(vecZetaP_AlbNot/2),-norminv(vecTuningP_AlbNot/2),[],'r','x');
scatter(-norminv(vecZetaP_Bl6Not/2),-norminv(vecTuningP_Bl6Not/2),[],'b','x');
hold off
fixfig;

%hip
dblR_AlbHip = nancorr(-norminv(vecZetaP_AlbHip/2),-norminv(vecTuningP_AlbHip/2));
dblR_Bl6Hip = nancorr(-norminv(vecZetaP_Bl6Hip/2),-norminv(vecTuningP_Bl6Hip/2));

subplot(2,3,3)
hold on;
scatter(-norminv(vecZetaP_AlbHip/2),-norminv(vecTuningP_AlbHip/2),[],'r','x');
scatter(-norminv(vecZetaP_Bl6Hip/2),-norminv(vecTuningP_Bl6Hip/2),[],'b','x');
hold off
fixfig;