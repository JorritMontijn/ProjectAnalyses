function cellStim = SyncET_NoEphysOld(sPupil,cellStim,cellRecNames,vecStartT)
	
	%% extract data
	vecPupilSyncTime = sPupil.vecPupilFullSyncLumT;
	vecFiltSyncLum = sPupil.vecPupilSyncLumFilt;
	vecPupilStimOn = sPupil.vecPupilStimOn;
	if ~exist('cellRecNames','var')
		cellRecNames = cellfill('N/A',size(cellStim));
	end
	
	%% list of recordings to ignore warnings
	cellIgnore = {'20191212_MP3_RunDriftingGratingsR01_16_34_41.mat','20191212_MP3_RunNaturalMovieR01_15_14_31.mat','20210212_MA7_16_35_17_RunDriftingGratings.mat','20210215_MA7_17_08_25_RunDriftingGratings.mat',...
		'20191213_MP3_RunNaturalMovieR01_15_49_31.mat'};
	
	%% run
	intLastPupilStop = 0;
	intLogs = numel(cellStim);
	for intLogFile = 1:intLogs
		%% calculate stimulus times
		fprintf('>Log file "%s" [%s]\n',cellRecNames{intLogFile},getTime);
		%return
		intThisNumTrials = numel(~isnan(cellStim{intLogFile}.structEP.ActOffSecs));
		
		%% align eye-tracking data
		%if intLogFile == 1,continue;end
		%get pupil on/offsets
		if ~exist('intLastPupilStop','var') || isempty(intLastPupilStop)
			intLastPupilStop = 0;
		end
		
		%% get pupil signals
		dblM = (1+4.7e-05);
		vecThisPupilOn = dblM*vecPupilStimOn(vecPupilStimOn>intLastPupilStop);
		
		%get stim signals
		vecThisStimOn = cellStim{intLogFile}.structEP.ActOnSecs;
		
		vecRefT = vecThisStimOn;
		vecSignalT = vecThisPupilOn;
		if ~exist('vecStartT','var') || numel(vecStartT) < intLogFile
			[dblStartT,intFlagOut] = SyncEvents(vecThisStimOn,vecThisPupilOn,true);
		else
			dblStartT = vecStartT(intLogFile);
		end
		%get aligment over trials
		dblOffset = 0.5*median(diff(vecRefT));
		vecCorrPupilOn = (vecThisStimOn - vecThisStimOn(1))/dblM + dblStartT;
		if isnan(dblStartT),continue;end
		[vecPlotT,matTracePerTrial] = getTraceInTrial(vecPupilSyncTime,vecFiltSyncLum,vecCorrPupilOn - dblOffset,median(diff(vecPupilSyncTime)),1.5*median(diff(vecRefT)));
		if 0
		figure
		hAx=gca;
		plot(hAx,vecPlotT,matTracePerTrial(1:100,:));
		hAx.ColorOrder = redbluepurple;
		end
		
		%realign
		[cellUpDown,cellCrit]=dimfun(1,@DP_GetUpDown,matTracePerTrial,0.6,1);
		matUD = cell2mat(cellUpDown);
		[cellOnset]=dimfun(1,@find,matUD,1,'first');
		vecOnset = cell2mat(cellOnset);
		vecOnsetT = vecPlotT(vecOnset) - dblOffset;
		cellStim{intLogFile}.structEP.ActOnPupil = vecCorrPupilOn + vecOnsetT;
		
		%%
		%calculate signal
		matR = corr(matTracePerTrial');
		vecMeanR = mean(matR);
		mdl = fitlm(1:numel(vecMeanR),vecMeanR);
		dblP = table2array(mdl.Coefficients(2,4));
		dblR2 = mdl.Rsquared.Adjusted;
		if dblP < 0.05 && dblR2 > 0.08 && ~ismember(cellRecNames{intLogFile},cellIgnore)
			error([mfilename ':SyncError'],'Inter-trial synchronization pulse correlation significantly depends on time: this indicates an incorrect alignment, signal decay, non-stationarity, or another problem');
		end
		
		%increment last stop
		intLastPupilStop = vecCorrPupilOn(end);
	end
