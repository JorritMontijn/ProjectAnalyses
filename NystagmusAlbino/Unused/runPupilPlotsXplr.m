%% load data
strDataPath = 'F:\Data\Processed\Neuropixels';
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

%% plot
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
figure;
for intSubType=1:2
	if intSubType == 1
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx);
			cellStim = sRec.cellBlock{intBlock};
			
			vecStimOn = cellStim.vecPupilStimOnTime;
			if isfield(sRec.sPupil,'vecPupilTimeNI')
				vecPupilT = sRec.sPupil.vecPupilTimeNI+dblOffsetT;
			else
				vecPupilT = sRec.sPupil.vecPupilTime+dblOffsetT;
			end
			
			%split by ori
			sTrialObjects = cellStim.sStimObject(cellStim.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24,continue,end
			
			%% pupil loc
			vecPupilLocX = sRec.sPupil.vecPupilFixedCenterX;
			vecPupilLocY = sRec.sPupil.vecPupilFixedCenterY;
			vecPupilSize = sRec.sPupil.vecPupilFixedRadius;
			vecPupilSync = sRec.sPupil.vecPupilSyncLum;
			%if ~isfield(sRec.sPupil,'vecPupilTimeNI')
			%doPEP(vecPupilT,sRec.sPupil.vecPupilSyncLum,vecStimOn);pause
			%end
			
			%% are they tracking the stimulus? i.e., is more time spent moving in the direction of the stimulus?
			%[matTE,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilLocX,vecStimOn,[-0.2 1.3]);
			vecWindowBinEdges = 0.05:0.02:1;
			[matTE_LocX] = getRespMat(vecPupilT,vecPupilLocX,vecStimOn,vecWindowBinEdges);
			[matTE_LocY] = getRespMat(vecPupilT,vecPupilLocY,vecStimOn,vecWindowBinEdges);
			vecWindowBinEdgesSize = -0.2:0.02:1.3;
			[matTE_Size,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSize,vecStimOn,vecWindowBinEdgesSize);
			[matTE_Sync,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSync,vecStimOn,vecWindowBinEdgesSize);
			matAggTE_LocX = cat(1,matAggTE_LocX,matTE_LocX);
			matAggTE_LocY = cat(1,matAggTE_LocY,matTE_LocY);
			matAggTE_Size = cat(1,matAggTE_Size,matTE_Size);
			matAggTE_Sync = cat(1,matAggTE_Sync,matTE_Sync);
			vecAggOriIdx = cat(1,vecAggOriIdx,vecOriIdx);
			vecAggCounts = vecAggCounts+vecCounts;
		end
	end
	%% plot
	intType=2;
	if intType == 1
	%get mean per ori, X
	subplot(2,4,1+(intSubType-1)*4)
	matTEX_dt = diff(matAggTE_LocX,[],2);
	matMove_MuX=splitapply(@nanmean,matTEX_dt,vecAggOriIdx);
	matNansX=splitapply(@(x) sum(isnan(x)),matTEX_dt,vecAggOriIdx);
	vecNanFrac=mean(matNansX,2)./vecAggCounts;
	vecEligibleX = (1 - vecNanFrac).*size(matMove_MuX,2);
	vecFracMoveRight = sum(matMove_MuX>0,2)./vecEligibleX;
	plot(vecUnique,vecFracMoveRight)
	
	%get mean per ori, Y
	subplot(2,4,2+(intSubType-1)*4)
	matTEY_dt = diff(matAggTE_LocY,[],2);
	matMove_MuY=splitapply(@nanmean,matTEY_dt,vecAggOriIdx);
	matNansY=splitapply(@(x) sum(isnan(x)),matTEY_dt,vecAggOriIdx);
	vecNanFrac=mean(matNansY,2)./vecAggCounts;
	vecEligibleY = (1 - vecNanFrac).*size(matMove_MuY,2);
	vecFracMoveDown = sum(matMove_MuY>0,2)./vecEligibleY;
	plot(vecUnique,vecFracMoveDown)
	
	
	%get mean per time
	subplot(2,4,3+(intSubType-1)*4)
	matAggTE_SizeZ = matAggTE_Size;%(matAggTE_Size-nanmean(matAggTE_Size(:)))/nanstd(matAggTE_Size(:));
	matSizeMu=nanmean(matAggTE_SizeZ,1);
	matSizeSd=nanstd(matAggTE_SizeZ,[],1);
	vecNans = sum(~isnan(matAggTE_SizeZ),1);
	errorbar(vecWindowBinCenters,matSizeMu,matSizeSd./sqrt(vecNans))
	hold on
	subplot(2,4,4+(intSubType-1)*4)
	matAggTE_SyncZ = (matAggTE_Sync-nanmean(matAggTE_Sync(:)))/nanstd(matAggTE_Sync(:));
	matSyncMu=nanmean(matAggTE_SyncZ,1);
	matSyncSd=nanstd(matAggTE_SyncZ,[],1);
	vecNans = sum(~isnan(matAggTE_SyncZ),1);
	errorbar(vecWindowBinCenters,matSyncMu,matSyncSd./sqrt(vecNans))
	hold off
	
	
	else
		%% plot v2
		%split ori in quadrants
		vecBoundQ = 45:90:360;
		vecOris = mod(vecUnique + 0,360);
		vecQ1 = vecOris(vecAggOriIdx) <= vecBoundQ(1) | vecOris(vecAggOriIdx) > vecBoundQ(end);
		vecQ2 = vecOris(vecAggOriIdx) <= vecBoundQ(2) & vecOris(vecAggOriIdx) > vecBoundQ(1);
		vecQ3 = vecOris(vecAggOriIdx) <= vecBoundQ(3) & vecOris(vecAggOriIdx) > vecBoundQ(2);
		vecQ4 = vecOris(vecAggOriIdx) <= vecBoundQ(4) & vecOris(vecAggOriIdx) > vecBoundQ(3);
		
		%get mean per ori, X
	subplot(2,4,1+(intSubType-1)*4)
	matTEX_dt = diff(matAggTE_LocX,[],2);
	matFracMoveRight = sum(matTEX_dt>0,2)./size(matTEX_dt,2);
	matFracNormMoveRight = matFracMoveRight / mean(matFracMoveRight);
	vecFMR_Q1 = matFracNormMoveRight(vecQ1);
	vecFMR_Q2 = matFracNormMoveRight(vecQ2);
	vecFMR_Q3 = matFracNormMoveRight(vecQ3);
	vecFMR_Q4 = matFracNormMoveRight(vecQ4);
	vecMeans = [mean(vecFMR_Q1) mean(vecFMR_Q2) mean(vecFMR_Q3) mean(vecFMR_Q4)];
	vecSems = [std(vecFMR_Q1)/sqrt(numel(vecFMR_Q1)) ...
		std(vecFMR_Q2)/sqrt(numel(vecFMR_Q2)) ...
		std(vecFMR_Q3)/sqrt(numel(vecFMR_Q3)) ...
		std(vecFMR_Q4)/sqrt(numel(vecFMR_Q4))];
	subplot(2,4,3+(intSubType-1)*4)
	errorbar(0:90:270,vecMeans,vecSems,'xb');
	
	
		%get mean per ori, Y
	subplot(2,4,1+(intSubType-1)*4)
	matTEY_dt = diff(matAggTE_LocY,[],2);
	matFracMoveUp = sum(matTEY_dt>0,2)./size(matTEY_dt,2);
	matFracNormMoveUp = matFracMoveUp / mean(matFracMoveUp);
	vecFMU_Q1 = matFracNormMoveUp(vecQ1);
	vecFMU_Q2 = matFracNormMoveUp(vecQ2);
	vecFMU_Q3 = matFracNormMoveUp(vecQ3);
	vecFMU_Q4 = matFracNormMoveUp(vecQ4);
	vecMeans = [mean(vecFMU_Q1) mean(vecFMU_Q2) mean(vecFMU_Q3) mean(vecFMU_Q4)];
	vecSems = [std(vecFMU_Q1)/sqrt(numel(vecFMU_Q1)) ...
		std(vecFMU_Q2)/sqrt(numel(vecFMU_Q2)) ...
		std(vecFMU_Q3)/sqrt(numel(vecFMU_Q3)) ...
		std(vecFMU_Q4)/sqrt(numel(vecFMU_Q4))];
	subplot(2,4,4+(intSubType-1)*4)
	errorbar(0:90:270,vecMeans,vecSems,'xb');
	
	return
	end
end