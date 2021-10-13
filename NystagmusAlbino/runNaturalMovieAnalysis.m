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
strTargetPath = 'F:\Data\Results\AlbinoProject';

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% define area categories
%cortex
cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
cellUseAreas{3} = {'Hippocampal formation','CA1','CA2','CA3','subiculum','dentate gyrus'};
cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};

%% pre-allocate
			
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
			
			%% prep data
			%get data matrix
			for intType=1:2
				if intType == 1
					figure;
					vecStimOnTime = sRec.sSources.cellBlock{intBlock}.structEP.ActOnNI(~indRemTrials) - dblFirstSamp/dblSampNi;
					strTit = 'NI-time';
				else
					vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
					strTit = 'Stim-time';
				end
				if numel(vecStimOnTime) <= 10,close;continue;end
			intPopCounter = intPopCounter + 1;
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			
			%include?
			vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
			indUseCells = vecZetaP(:)<0.05 & arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
			
			%% split cells into areas
			%build cell vectors
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			end
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			%% split data into bins
			intTypeCV = 2;
			dblLambda = 100;
			intBinNr = 20;
			dblMovieDur = roundi(median(vecStimOffTime-vecStimOnTime)*2,0);
			dblBinDur = dblMovieDur/intBinNr;
			vecBinOnset = linspace(0,dblMovieDur-dblBinDur,intBinNr);
			
			%build artificial stim vector
			intRepNum = numel(vecStimOnTime);
			vecBinOnT = flat(vecStimOnTime(:)' + vecBinOnset(:));
			vecBinIdx = flat(repmat((1:intBinNr)',[1 intRepNum]));
			matData = getSpikeCounts(cellSpikeT,vecBinOnT,dblBinDur);
			
			% decode movie
			for intArea = 1%:numel(cellUseAreas)
			%% select cells
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				%% plot
				%[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(cellSpikeT{vecSelectCells(5)},vecStimOnTime,20,[],3)
			
				%% decode movie
				matUseData = matData(vecSelectCells,:);
				vecPriorDistribution = intRepNum*ones(1,intBinNr);
				[dblPerformanceCV,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionML] = doCrossValidatedDecodingML(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution);
				[dblPerformanceTM,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusionTM] = doCrossValidatedDecodingTM(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution);
				[dblPerformanceLR,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionLR] = ...
					doCrossValidatedDecodingLR(matUseData,vecBinIdx,intTypeCV,[],dblLambda);
				[dblPerformanceLR2,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionLR2] = ...
					doCrossValidatedDecodingLR(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution,dblLambda);
				
				
				subplot(2,4,1+(intType-1)*4)
				imagesc(matConfusionML)
				title([strTit '; ML: ' strName '_' num2str(intBlock)],'interpreter','none'); 
				
				subplot(2,4,2+(intType-1)*4)
				imagesc(matConfusionTM)
				title(['TM: ' strName '_' num2str(intBlock)],'interpreter','none'); 
			
				subplot(2,4,3+(intType-1)*4)
				imagesc(matConfusionLR)
				title(['LR: ' strName '_' num2str(intBlock)],'interpreter','none'); 
				
				subplot(2,4,4+(intType-1)*4)
				imagesc(matConfusionLR2)
				title(['LR2: ' strName '_' num2str(intBlock)],'interpreter','none'); 
				
				if intType == 2
					maxfig;drawnow;
					export_fig(fullpath(strTargetPath,['NatMovDecoding' getDate '.jpg']));
					export_fig(fullpath(strTargetPath,['NatMovDecoding' getDate '.pdf']));
				end
			end
		end
	end
	end
	%% plot
	
end
%% save
return
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.jpg']));
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.pdf']));