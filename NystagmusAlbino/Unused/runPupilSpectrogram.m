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
boolSavePlots = true;

%% plot
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
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
		if ~isfield(sRec,'sPupil'),continue;end
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
			
			indRemVals = [false diff(vecPupilT)<=0];
			dblStartT = vecPupilStimOn(1)-30;
			dblStopT = vecPupilStimOn(end)+30;
			indUseSamples = vecPupilT>dblStartT & vecPupilT<dblStopT & ~indRemVals;
			if all(indUseSamples==0),continue;end
			
			%% pupil loc
			vecPupilTime = vecPupilT(indUseSamples);
			vecPupilLocX = sRec.sPupil.vecCenterX(indUseSamples);
			vecPupilLocY = sRec.sPupil.vecCenterY(indUseSamples);
			vecPupilSize = sRec.sPupil.vecRadius(indUseSamples);
			vecPupilSync = sRec.sPupil.vecSyncLum(indUseSamples);
			if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
				vecPupilBlinks = sRec.sPupil.vecBlinks(indUseSamples);
			else
				%filter absvidlum
				dblLowPass = 0.01/(1/median(diff(vecPupilT)));
				[fb,fa] = butter(2,dblLowPass,'high');
				vecPupilAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum(indUseSamples)));
				vecPupilBlinks = vecPupilAbsVidLum > 5;
			end
			dblFs = 1/median(diff(vecPupilT));
			
			%% spectrogram
			[s,f,t,m] = spectrogram(vecPupilLocY,floor(dblStopT-dblStartT),[],[],dblFs);
			mLog = log(m);
			%mLog = mLog - mean(mLog,2);
			
			figure
			subplot(2,3,[1 2 4 5])
			imagesc(t,f(2:end),mLog(2:end,:))
			title([sRec.sJson.subject '_' sRec.sJson.date 'B' num2str(intBlock) ],'interpreter','none')
			xlabel('Time (s)');
			ylabel('Frequency (Hz)');
			h=colorbar;
			clabel(h,'Power (dB)')
			fixfig;
			grid off;
			axis xy
			
			subplot(2,3,3)
			vecP = imfilt(log(mean(m,2)),normpdf(-2:2,0,1)');
			plot(f,vecP - median(vecP))
			ylim([-5 10])
			xlabel('Frequency (Hz)');
			ylabel('Normalized power (dB)');
			fixfig;
			drawnow;maxfig;
			
			%% save
			if boolSavePlots
				drawnow;
				export_fig(fullpath(strTargetPath,['NatMovSpectrum_' sRec.sJson.subject '_' sRec.sJson.date  'B' num2str(intBlock) '.jpg']));
				export_fig(fullpath(strTargetPath,['NatMovSpectrum_' sRec.sJson.subject '_' sRec.sJson.date  'B' num2str(intBlock) '.pdf']));
			end
			%{
			%% v3
			%split vector into trials, concat as matrix
			%remove trials with blinks
			%pspectrum per trial
			%average over trials
			dblFs = 1/median(diff(vecPupilT));
			intTrials = numel(vecStimOn);
			dblDur = median(vecStimOff - vecStimOn);
			intSampN = round(dblDur*dblFs);
			matEyeX = nan(intTrials,intSampN);
			matEyeY = nan(intTrials,intSampN);
			matEyeB = nan(intTrials,intSampN);
			for intTrial=1:intTrials
				intStart = find(vecPupilTime > vecStimOn(intTrial),1);
				matEyeX(intTrial,:) = vecPupilLocX(intStart:(intStart+intSampN-1));
				matEyeY(intTrial,:) = vecPupilLocY(intStart:(intStart+intSampN-1));
				matEyeB(intTrial,:) = vecPupilBlinks(intStart:(intStart+intSampN-1));
			end
			indRemTrials = (sum(matEyeB,2)./intSampN) > 0;
			matEyeX(indRemTrials,:) = [];
			matEyeY(indRemTrials,:) = [];
			%pspec
			intTrials = size(matEyeX,1);
			for intTrial=1:intTrials
				[p,vecF] = pspectrum(matEyeX(intTrial,:),dblFs);
				if intTrial==1
					matSpecX = nan(intTrials,numel(vecF));
					matSpecY = nan(intTrials,numel(vecF));
				end
				matSpecX(intTrial,:) = log(p);
				[p,vecF] = pspectrum(matEyeY(intTrial,:),dblFs);
				matSpecY(intTrial,:) = log(p);
			end
			vecSpec = (mean(matSpecY,1) + mean(matSpecY,1))/2;
			figure
			%plot(vecF,matSpecY)
			plot(vecF,vecSpec)
			title([sRec.sJson.subject '_' sRec.sJson.date],'interpreter','none')
			%}
		end
	end
	
end