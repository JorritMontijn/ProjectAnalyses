%% further analyses
%2) is neural code of nat movs more variable during eye movements in NOT than Ctx?

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
cellAreaGroups = {'Vis. ctx','NOT'};
cellAreaGroupsAbbr = {'Ctx','NOT'};
cellSubjectGroups = {'BL6','DBA'};

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
cellAreasPerExp = cell(1,numel(sExp));
cellStimsPerExp = cell(1,numel(sExp));
for intExp=1:numel(sExp)
	sRec = sExp(intExp);
	strName=[sRec.sJson.subject '_' sRec.sJson.date];
	cellAreasPerExp{intExp} = unique({sExp(intExp).sCluster.Area})';
	cellStimsPerExp{intExp} = [strName ':' char(cell2vec(cellfun(@(x) x.strExpType,sExp(intExp).cellBlock,'UniformOutput',false)))'];
end

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
		vecBlocksDG = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
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
			indRemTrials = vecBlinkFractionPerTrial > 0;
			
			%get movement
			vecEyeMovement = sqrt(diff(vecPupilLocY).^2 + diff(vecPupilLocX).^2);
			
			%% prep data
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
				vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
				
				%% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				%get data matrix
				vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
				vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
				
				if numel(vecStimOnTime) <= 10,close;continue;end
				intPopCounter = intPopCounter + 1;
				cellSpikeT = {sRec.sCluster(:).SpikeTimes};
				
				%include?
				vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
				%indUseCells = vecZetaP(:)<0.05 & arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
				indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
				
				%% split cells into areas
				%build cell vectors
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				
				%% split data into bins
				intTypeCV = 2;
				dblLambda = 100;
				dblStimDur = median(vecStimOffTime-vecStimOnTime);
				if dblStimDur < 9
					dblMovieDur = 8; %new; actually 8.33
					intBinNr = 8;
				else
					dblMovieDur = 20; %old
					intBinNr = 20;
				end
				dblBinDur = dblMovieDur/intBinNr;
				vecBinOnset = linspace(0,dblMovieDur-dblBinDur,intBinNr);
				
				%build artificial stim vector
				intRepNum = numel(vecStimOnTime);
				vecBinOnT = flat(vecStimOnTime(:)' + vecBinOnset(:));
				vecBinIdx = flat(repmat((1:intBinNr)',[1 intRepNum]));
				matData = getSpikeCounts(cellSpikeT,vecBinOnT,dblBinDur);
				
				%% decode movie
				matUseData = matData(vecSelectCells,:);
				vecPriorDistribution = intRepNum*ones(1,intBinNr);
				[dblPerformanceLR2,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionLR2] = ...
					doCrossValidatedDecodingLR(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution,dblLambda);
				
				%get movement & probability of correct time bin
				intMovieBins = size(matPostProbability,2);
				vecProbCorrect = nan(1,intMovieBins);
				vecMovePupil = nan(1,intMovieBins);
				for intMovieBin=1:intMovieBins
					%get probability
					vecProbCorrect(intMovieBin) = matPostProbability(vecBinIdx(intMovieBin),intMovieBin);
					
					%eye movement
					vecMovePupil(intMovieBin) = [];
				end
				
				
				subplot(2,4,4+(intType-1)*4)
				imagesc(matConfusionLR2)
				title(['LR2: ' strName '_' num2str(intBlock)],'interpreter','none');
				return
			end
		end
	end
	%% plot
	
end
%% save
