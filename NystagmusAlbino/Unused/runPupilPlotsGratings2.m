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
		dblAverageMouseHeadTiltInSetup = -15; %average tilt due to head-bar placement; procedures changed between BL6 and DBA experiments
		boolInvertX = 1; %other eye was recorded, so temporonasal is nasotemporal
	elseif intSubType == 2
		strSubjectType = 'DBA';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = -15;
		boolInvertX = 0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecTrialCounts = [];
	matFracMoveRight = zeros(24,0);
	matFracMoveUp = zeros(24,0);
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		
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
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			vecPupilStimOn = sBlock.vecStimOnTime+median(vecPupilLatency);
			vecPupilStimOff = sBlock.vecStimOffTime+median(vecPupilLatency);
			vecPupilT = sRec.sPupil.vecTime;
			
			%split by ori
			sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24,continue,end
			
			%% pupil loc
			vecPupilLocX = sRec.sPupil.vecCenterX;
			vecPupilLocY = sRec.sPupil.vecCenterY;
			vecPupilSize = sRec.sPupil.vecRadius;
			vecPupilSync = sRec.sPupil.vecSyncLum;
			if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
				vecPupilBlinks = sRec.sPupil.vecBlinks;
			else
				%filter absvidlum
				dblLowPass = 0.01/(1/median(diff(vecPupilT)));
				[fb,fa] = butter(2,dblLowPass,'high');
				vecPupilAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum));
				vecPupilBlinks = vecPupilAbsVidLum > 5;
			end
			dblFs = 1/median(diff(vecPupilT));
			
			%% are they tracking the stimulus? i.e., is more time spent moving in the direction of the stimulus?
			%[matTE,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilLocX,vecStimOn,[-0.2 1.3]);
			vecWindowBinEdges = 0.05:0.02:1;
			[matTE_LocX] = getRespMat(vecPupilT,vecPupilLocX,vecPupilStimOn,vecWindowBinEdges);
			[matTE_LocY] = getRespMat(vecPupilT,vecPupilLocY,vecPupilStimOn,vecWindowBinEdges);
			[matTE_Blink] = getRespMat(vecPupilT,double(vecPupilBlinks),vecPupilStimOn,vecWindowBinEdges);
			matTE_Blink(isnan(matTE_Blink))=0;
			vecBlinkFrac = sum(matTE_Blink==1,2)./sum(~isnan(matTE_Blink),2);
			vecWindowBinEdgesSize = -0.2:0.02:1.3;
			[matTE_Size,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSize,vecPupilStimOn,vecWindowBinEdgesSize);
			[matTE_Sync,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSync,vecPupilStimOn,vecWindowBinEdgesSize);
			%remove trials with >1/5 of blinking
			indRemTrials = vecBlinkFrac>(1/5);
			matTE_LocX(indRemTrials,:) = [];
			matTE_LocY(indRemTrials,:) = [];
			matTE_Size(indRemTrials,:) = [];
			matTE_Sync(indRemTrials,:) = [];
			vecTheseOris = vecOrientation(~indRemTrials);
			vecTheseOris = mod(vecTheseOris+dblAverageMouseHeadTiltInSetup,360);
			%add data
			matAggTE_LocX = cat(1,matAggTE_LocX,matTE_LocX);
			matAggTE_LocY = cat(1,matAggTE_LocY,matTE_LocY);
			matAggTE_Size = cat(1,matAggTE_Size,matTE_Size);
			matAggTE_Sync = cat(1,matAggTE_Sync,matTE_Sync);
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTheseOris);
			vecAggOriIdx = cat(1,vecAggOriIdx,vecOriIdx);
			
			vecAggCounts = vecAggCounts+vecCounts;
			
			%calc right movement
			matTEX_dt = diff(matTE_LocX,[],2);
			if boolInvertX,matTEX_dt=-matTEX_dt;end
			matMoveR_MuX=splitapply(@(x) nansum(x>0),matTEX_dt,vecOriIdx); %0 or mean?
			matMoveL_MuX=splitapply(@(x) nansum(x<0),matTEX_dt,vecOriIdx); %0 or mean?
			matMoveR_Frac = matMoveR_MuX ./ (matMoveL_MuX+matMoveR_MuX);
			vecFracMoveRight = mean(matMoveR_Frac,2);
			matFracMoveRight(:,end+1) = vecFracMoveRight;
			if any(vecFracMoveRight>1),error;end
			
			%calc left movement
			matTEY_dt = diff(matTE_LocY,[],2);
			matMoveU_MuY=splitapply(@(x) nansum(x>0),matTEY_dt,vecOriIdx); %0 or mean?
			matMoveD_MuY=splitapply(@(x) nansum(x<0),matTEY_dt,vecOriIdx); %0 or mean?
			matMoveU_Frac = matMoveU_MuY ./ (matMoveD_MuY+matMoveU_MuY);
			vecFracMoveUp = mean(matMoveU_Frac,2);
			matFracMoveUp(:,end+1) = vecFracMoveUp;
			
			vecTrialCounts(end+1) = size(matTE_LocX,1);
		end
	end
	
	%% plot
	%get mean per ori, X
	subplot(2,4,1+(intSubType-1)*4)
	matTEX_dt = diff(matAggTE_LocX,[],2);
	if boolInvertX,matTEX_dt=-matTEX_dt;end
	matMoveR_MuX=splitapply(@(x) nansum(x>0),matTEX_dt,vecAggOriIdx); %0 or mean?
	matMoveL_MuX=splitapply(@(x) nansum(x<0),matTEX_dt,vecAggOriIdx); %0 or mean?
	matMoveR_Frac = matMoveR_MuX ./ (matMoveL_MuX+matMoveR_MuX);
	vecFracMoveRight = mean(matMoveR_Frac,2);
	vecFracMoveRight = (vecFracMoveRight + circshift(vecFracMoveRight,-1) + circshift(vecFracMoveRight,1))/3;
	vecFracMoveRightMinMax = [min(vecFracMoveRight) max(vecFracMoveRight)];
	vecNormFracRight = vecFracMoveRight-min(vecFracMoveRight);
	vecNormFracRight = vecNormFracRight./max(vecNormFracRight(:));
	vecOris = (0:15:359)';
	polarplot([deg2rad(vecOris); deg2rad(vecOris(1))],[vecFracMoveRight; vecFracMoveRight(1)],'r');
	rlim(vecFracMoveRightMinMax);%vecFracMoveUpMinMax)
	title(sprintf('%s; Fraction of time eye is moving right',strSubjectType))
	%set(gca,'ytick',[0 1],'yticklabel',vecFracMoveRightMinMax)
	fixfig;
	
	%get mean per ori, Y
	subplot(2,4,2+(intSubType-1)*4)
	matTEY_dt = -diff(matAggTE_LocY,[],2);
	matMoveU_MuY=splitapply(@(x) nansum(x>0),matTEY_dt,vecAggOriIdx); %0 or mean?
	matMoveD_MuY=splitapply(@(x) nansum(x<0),matTEY_dt,vecAggOriIdx);%0 or mean?
	matMoveU_Frac = matMoveU_MuY ./ (matMoveD_MuY+matMoveU_MuY);
	vecFracMoveUp = mean(matMoveU_Frac,2);
	
	vecFracMoveUp = (vecFracMoveUp + circshift(vecFracMoveUp,-1) + circshift(vecFracMoveUp,1))/3;
	vecFracMoveUpMinMax = [min(vecFracMoveUp) max(vecFracMoveUp)];
	vecNormFracUp = vecFracMoveUp-min(vecFracMoveUp);
	vecNormFracUp = vecNormFracUp./max(vecNormFracUp(:));
	h=polarplot([deg2rad(vecOris); deg2rad(vecOris(1))],[vecFracMoveUp; vecFracMoveUp(1)],'b');
	rlim(vecFracMoveUpMinMax);%vecFracMoveUpMinMax)
	%set(gca,'ytick',[0 1],'yticklabel',vecFracMoveUpMinMax)
	%xlim([-15 360]);
	title('Fraction of time eye is moving up')
	fixfig;
	
	%{
	subplot(2,4,4+(intSubType-1)*4)
	matAggTE_SyncZ = (matAggTE_Sync-nanmean(matAggTE_Sync(:)))/nanstd(matAggTE_Sync(:));
	matSyncMu=nanmean(matAggTE_SyncZ,1);
	matSyncSd=nanstd(matAggTE_SyncZ,[],1);
	vecNans = sum(~isnan(matAggTE_SyncZ),1);
	errorbar(vecWindowBinCenters,matSyncMu,matSyncSd./sqrt(vecNans))
	hold off
	%}
	
	% plot v2
	%{
	%split ori in quadrants
	vecBoundQ = 45:90:360; %0 is rightward, 90 is up
	vecOris = mod(vecOris + 0,360);
	vecQ1 = vecOris(vecAggOriIdx) < vecBoundQ(1) | vecOris(vecAggOriIdx) > vecBoundQ(end); %right
	vecQ2 = vecOris(vecAggOriIdx) < vecBoundQ(2) & vecOris(vecAggOriIdx) > vecBoundQ(1); %up
	vecQ3 = vecOris(vecAggOriIdx) < vecBoundQ(3) & vecOris(vecAggOriIdx) > vecBoundQ(2); %left
	vecQ4 = vecOris(vecAggOriIdx) < vecBoundQ(4) & vecOris(vecAggOriIdx) > vecBoundQ(3); %down
	
	%get mean per ori, X
	matTEX_dt = diff(matAggTE_LocX,[],2);
	dblCutOff = 0;%nanmean(matTEX_dt(:));
	matFracMoveRight = sum(matTEX_dt>dblCutOff,2)./size(matTEX_dt,2);
	matFracNormMoveRight = 1-(matFracMoveRight / mean(matFracMoveRight));
	vecFMR_Q1 = matFracNormMoveRight(vecQ1);
	vecFMR_Q2 = matFracNormMoveRight(vecQ2);
	vecFMR_Q3 = matFracNormMoveRight(vecQ3);
	vecFMR_Q4 = matFracNormMoveRight(vecQ4);
	%}
	%plot v3
	vecFracNormMoveRight = 1-(vecFracMoveRight / mean(vecFracMoveRight));
	vecFMR_Q1 = vecFracNormMoveRight((vecOris < 90) | vecOris > 270);
	vecFMR_Q3 = vecFracNormMoveRight((vecOris > 90) & vecOris < 270);
	vecMeans = [mean(vecFMR_Q1) mean(vecFMR_Q3)];
	vecSems = [std(vecFMR_Q1)/sqrt(numel(vecFMR_Q1)) ...
		std(vecFMR_Q3)/sqrt(numel(vecFMR_Q3)) ...
		];
	
	subplot(2,4,3+(intSubType-1)*4)
	errorbar([0.2 0.8],vecMeans,vecSems,'xr');
	[h,pLR] = ttest2(vecFMR_Q1,vecFMR_Q3);
	title(sprintf('L vs R, p=%.3e',pLR));
	ylim([-0.06 0.06]);
	xlim([0 1]);
	set(gca,'xtick',[0.2 0.8],'xticklabel',{'Left','Right'});
	ylabel('Norm. frac. eye moving right');
	xlabel('Stim. dir');
	fixfig;
	
	%get mean per ori, Y
	%{
	matTEY_dt = -diff(matAggTE_LocY,[],2);
	dblCutOff = 0;%nanmean(matTEY_dt(:));
	matFracMoveUp = sum(matTEY_dt>dblCutOff,2)./size(matTEY_dt,2);
	matFracNormMoveUp = 1-(matFracMoveUp / mean(matFracMoveUp));
	vecFMU_Q1 = matFracNormMoveUp(vecQ1);
	vecFMU_Q2 = matFracNormMoveUp(vecQ2);
	vecFMU_Q3 = matFracNormMoveUp(vecQ3);
	vecFMU_Q4 = matFracNormMoveUp(vecQ4);
	vecMeans = [mean(vecFMU_Q2)  mean(vecFMU_Q4)];
	vecSems = [...
		std(vecFMU_Q2)/sqrt(numel(vecFMU_Q2)) ...
		std(vecFMU_Q4)/sqrt(numel(vecFMU_Q4))];
	%}
	%plot v3
	vecFracNormMoveUp = 1-(vecFracMoveUp / mean(vecFracMoveUp));
	vecFMU_Q2 = vecFracNormMoveUp(vecOris > 0 & vecOris < 180);
	vecFMU_Q4 = vecFracNormMoveUp(vecOris > 0 & vecOris > 180);
	vecMeans = [mean(vecFMU_Q2) mean(vecFMU_Q4)];
	vecSems = [std(vecFMU_Q2)/sqrt(numel(vecFMU_Q2)) ...
		std(vecFMU_Q4)/sqrt(numel(vecFMU_Q4)) ...
		];
	
	%subplot(2,16,[14 15 16]+(intSubType-1)*16)
	subplot(2,4,4+(intSubType-1)*4)
	errorbar([0.2 0.8],vecMeans,vecSems,'xb');
	[h,pUD] = ttest2(vecFMU_Q2,vecFMU_Q4);
	title(sprintf('U vs D, p=%.3e',pUD));
	ylim([-0.06 0.06]);
	xlim([0 1]);
	set(gca,'xtick',[0.2 0.8],'xticklabel',{'Down','Up'});
	ylabel('Norm. frac. eye moving up');
	xlabel('Stim. dir');
	fixfig;
	
	%% save data
	indRem=all(matFracMoveRight==0,1) | all(matFracMoveUp==0,1) | any(isnan(matFracMoveRight),1) | any(isnan(matFracMoveRight),1);
	cellFracMoveRight{intSubType} = matFracMoveRight(:,~indRem);
	cellFracMoveUp{intSubType} = matFracMoveUp(:,~indRem);
end

%% save plot
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['GratingTracking.tif']));
export_fig(fullpath(strTargetPath,['GratingTracking.pdf']));

return
%% plot per rec
cellFracMoveRightFilt{1} = (cellFracMoveRight{1} + circshift(cellFracMoveRight{1},[-1 0]) + circshift(cellFracMoveRight{1},[1 0]))/3;
cellFracMoveRightFilt{2} = (cellFracMoveRight{2} + circshift(cellFracMoveRight{2},[-1 0]) + circshift(cellFracMoveRight{2},[1 0]))/3;
cellFracMoveUpFilt{1} = (cellFracMoveUp{1} + circshift(cellFracMoveUp{1},[-1 0]) + circshift(cellFracMoveUp{1},[1 0]))/3;
cellFracMoveUpFilt{2} = (cellFracMoveUp{2} + circshift(cellFracMoveUp{2},[-1 0]) + circshift(cellFracMoveUp{2},[1 0]))/3;
figure
subplot(2,3,1)
hold on
errorbar(vecOris,mean(cellFracMoveRightFilt{1},2),std(cellFracMoveRightFilt{1},[],2)/sqrt(size(cellFracMoveRightFilt{1},2)),'r');
errorbar(vecOris,mean(cellFracMoveUpFilt{1},2),std(cellFracMoveUpFilt{1},[],2)/sqrt(size(cellFracMoveUpFilt{1},2)),'b');
hold off

subplot(2,3,4)
hold on
errorbar(vecOris,mean(cellFracMoveRightFilt{2},2),std(cellFracMoveRightFilt{2},[],2)/sqrt(size(cellFracMoveRightFilt{2},2)),'r');
errorbar(vecOris,mean(cellFracMoveUpFilt{2},2),std(cellFracMoveUpFilt{2},[],2)/sqrt(size(cellFracMoveUpFilt{2},2)),'b');
hold off

vecBiasLR_BL6 = mean(cellFracMoveRight{1},1);
[h,pHB] = ttest(vecBiasLR_BL6,0.5);
vecBiasLR_DBA = mean(cellFracMoveRight{2},1);
[h,pHA] = ttest(vecBiasLR_DBA,0.5);

vecBiasUD_BL6 = mean(cellFracMoveUp{1},1);
[h,pVB] = ttest(vecBiasUD_BL6,0.5);
vecBiasUD_DBA = mean(cellFracMoveUp{2},1);
[h,pVA] = ttest(vecBiasUD_DBA,0.5);
