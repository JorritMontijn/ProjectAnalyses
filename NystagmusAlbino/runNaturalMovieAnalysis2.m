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
matAggCorrWt = nan(0,numel(cellUseAreas),1);
matAggCorrAlb = nan(0,numel(cellUseAreas),1);
matAggMoveWt = nan(0,numel(cellUseAreas),2);
matAggMoveAlb = nan(0,numel(cellUseAreas),2);

%% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);

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
			%vecEyeMovement = vecPupilSize;
			
			%% prep data
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			vecCorr = nan(1,numel(cellUseAreas));
			matMeanMov = nan(2,numel(cellUseAreas));
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
				vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
				
				%% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				%get data matrix
				vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
				vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
				if ~isfield(sBlock,'vecPupilStimOn'),continue;end
				vecPupilStimOn = sBlock.vecPupilStimOn(~indRemTrials);
				vecPupilStimOff = sBlock.vecPupilStimOff(~indRemTrials);
		
				if numel(vecStimOnTime) <= 10,close;continue;end
				intPopCounter = intPopCounter + 1;
				cellSpikeT = {sRec.sCluster(:).SpikeTimes};
				
				%include?
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
				[phat,pci] = binofit(dblPerformanceLR2*numel(vecDecodedIndexCV),numel(vecDecodedIndexCV));
				if pci(1) < (1/intBinNr),continue;end
				
				%get pupil starts
				vecBinOnPupilT = flat(vecPupilStimOn(:)' + vecBinOnset(:));
				
				%get movement & probability of correct time bin
				intMovieBins = size(matPostProbability,2);
				vecProbCorrect = nan(1,intMovieBins);
				vecIsCorrect = nan(1,intMovieBins);
				vecMovePupil = nan(1,intMovieBins);
				for intMovieBin=1:intMovieBins
					%get probability
					vecProbCorrect(intMovieBin) = matPostProbability(vecBinIdx(intMovieBin),intMovieBin);
					vecIsCorrect(intMovieBin) = vecDecodedIndexCV(intMovieBin)==vecBinIdx(intMovieBin);
					
					%get eye movement
					dblStart = vecBinOnPupilT(intMovieBin);
					dblStop = dblStart + dblBinDur;
					intStartT = max([1 find(vecPupilTime > dblStart,1) - 1]);
					intStopT = min([numel(vecPupilTime) find(vecPupilTime > dblStop,1) + 1]);
					vecSelectFrames = intStartT:intStopT;
					vecMovePupil(intMovieBin) = mean(vecEyeMovement(vecSelectFrames));
				end
				
				%test correct/incorrect difference
				vecMoveCorrect = vecMovePupil(vecIsCorrect==1);
				vecMoveIncorrect = vecMovePupil(vecIsCorrect==0);
				[hT,pT]=ttest2(vecMoveCorrect,vecMoveIncorrect);
				matMeanMov(intArea,1) = mean(vecMoveCorrect);
				matMeanMov(intArea,2) = mean(vecMoveIncorrect);
				
				%calculate correlation
				[R,P,RL,RU] = corrcoef(vecMovePupil',vecProbCorrect');
				vecCorr(intArea) = R(1,2);
				
				%figure
				%scatter(vecMovePupil,vecProbCorrect,'.')
				%xlabel('Pupil movement');
				%ylabel('P(correct)');
				%title(sprintf('%s,%s; %sB%d, R=%.3f, p=%.3f; T-p=%.3f',strSubjectType,strAreaGroup,strName,intBlock,R(1,2),P(1,2),pT),'interpreter','none');
				%fixfig;
				
			end
			
			if strcmp(strSubjectType,'BL6')
				matAggCorrWt(end+1,:,:) = vecCorr;
				matAggMoveWt(end+1,:,:) = matMeanMov;
				%matAggPerfWt(end+1,:,:) = matPerf;
				%cellAggPrefWt(end+1,:) = cellPref;
			else
				matAggCorrAlb(end+1,:,:) = vecCorr;
				matAggMoveAlb(end+1,:,:) = matMeanMov;
				%matAggPerfAlb(end+1,:,:) = matPerf;
				%cellAggPrefAlb(end+1,:) = cellPref;
			end
		end
	end
	%% plot
	
end
%% data
vecCorrCtxBL6 = matAggCorrWt(:,1);
vecCorrCtxBL6(isnan(vecCorrCtxBL6)) = [];
vecCorrNOTBL6 = matAggCorrWt(:,2);
vecCorrNOTBL6(isnan(vecCorrNOTBL6)) = [];

vecCorrCtxAlb = matAggCorrAlb(:,1);
vecCorrCtxAlb(isnan(vecCorrCtxAlb)) = [];
vecCorrNOTAlb = matAggCorrAlb(:,2);
vecCorrNOTAlb(isnan(vecCorrNOTAlb)) = [];

matDiffMoveAlb = matAggMoveAlb(:,:,1) - matAggMoveAlb(:,:,2)
vecDiffMoveAlbCtx = matDiffMoveAlb(:,1);
vecDiffMoveAlbNOT = matDiffMoveAlb(:,2);

matDiffMoveBL6 = matAggMoveWt(:,:,1) - matAggMoveWt(:,:,2)
vecDiffMoveBL6Ctx = matDiffMoveBL6(:,1);
vecDiffMoveBL6NOT = matDiffMoveBL6(:,2);

[h,p_Ctx2] = ttest2(vecDiffMoveBL6Ctx,vecDiffMoveAlbCtx);
[h,p_NOT2] = ttest2(vecDiffMoveBL6NOT,vecDiffMoveAlbNOT);

%% test
[h,p_Ctx] = ttest2(vecCorrCtxBL6,vecCorrCtxAlb);
[h,p_NOT] = ttest2(vecCorrNOTBL6,vecCorrNOTAlb);

%interaction
vecY = cat(1,vecCorrCtxBL6,vecCorrCtxAlb,vecCorrNOTBL6,vecCorrNOTAlb);
vecG1 = cat(1,ones(size(vecCorrCtxBL6)),2*ones(size(vecCorrCtxAlb)),ones(size(vecCorrNOTBL6)),2*ones(size(vecCorrNOTAlb))); %bl6 vs alb
cellG1 = cellSubjectGroups(vecG1)';
vecG2 = cat(1,ones(size(vecCorrCtxBL6)),ones(size(vecCorrCtxAlb)),2*ones(size(vecCorrNOTBL6)),2*ones(size(vecCorrNOTAlb))); %ctx vs not
cellG2 = cellAreaGroupsAbbr(vecG2)';

[p,tbl,stats,terms] = anovan(vecY,{cellG1,cellG2},'model','full','display','off');

%% plot
figure
errorbar([1 2],[mean(vecCorrCtxBL6) mean(vecCorrNOTBL6)],...
	[std(vecCorrCtxBL6)/sqrt(numel(vecCorrCtxBL6)) std(vecCorrNOTBL6)/sqrt(numel(vecCorrNOTBL6))],...
	'x-','Color',vecColBl6,'CapSize',20);
hold on
errorbar([1 2],[mean(vecCorrCtxAlb) mean(vecCorrNOTAlb)],...
	[std(vecCorrCtxAlb)/sqrt(numel(vecCorrCtxAlb)) std(vecCorrNOTAlb)/sqrt(numel(vecCorrNOTAlb))],...
	'x-','Color',vecColAlb,'CapSize',20);
hold off
ylabel(sprintf('Pearson R(Pupil movement,Decoding)'));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
legend(cellSubjectGroups,'Location','best');
title(sprintf('Ctx,p=%.4f,NOT,p=%.4f, Interaction,p=%.4f',p_Ctx,p_NOT,p(3)));
fixfig;grid off;drawnow;

export_fig([strTargetPath filesep sprintf('NatMovDecodingMovementCorr.tif')]);
export_fig([strTargetPath filesep sprintf('NatMovDecodingMovementCorr.pdf')]);
	