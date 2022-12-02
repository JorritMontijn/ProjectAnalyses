%% load data
strDataPath = 'E:\DataPreProcessed';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		if ~isfield(sLoad.sAP,'sPupil') || isempty(sLoad.sAP.sPupil),continue;end
		sLoad.sAP.name = sFiles(intFile).name;
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

%% plot
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
indRemoveNoEye = ~contains(strrep(cellExperiment,'-',''),cellUseForEyeTracking);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
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
	indUseRecs = contains(cellSubjectType,strSubjectType) & ~indRemoveNoEye
	vecRunRecs = find(indUseRecs);
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		sRec.name = sFiles(intFile).name;
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false)
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true))
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			%{
			if isfield(cellStim,'vecPupilStimOnTime')
				vecPupilTime = sRec.sPupil.vecPupilTimeime+dblOffsetT;
				vecStimOn = cellStim.vecPupilStimOnTime;
				vecStimOff = cellStim.vecPupilStimOffTime;
			elseif isfield(sRec.sPupil,'vecPupilTimeimeNI')
				vecPupilTime = sRec.sPupil.vecPupilTimeimeNI+dblOffsetT;
				vecStimOn = cellStim.vecStimOnTime;
				vecStimOff = cellStim.vecStimOffTime;
			else
				continue;
			end
			%}
			
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
			
			%% get rf resps
			%{
			indNOT = contains({sRec.sCluster.Area},'Nucleus of the optic tract','IgnoreCase',true);
			sNeuronsNOT = sRec.sCluster(indNOT);
			
			intNeuron = 1;
			vecSpikes = sNeuronsNOT(intNeuron).SpikeTimes;
			%}
			%Rfs
			if isfield(sBlock,'vecPupilStimOn')
			vecPupilStimOn = sBlock.vecPupilStimOn;
			vecPupilStimOff = sBlock.vecPupilStimOff;
			else
			vecPupilStimOn = sBlock.vecStimOnTime;
			vecPupilStimOff = sBlock.vecStimOffTime;
			end
			sStimObject = sBlock.sStimObject;
			intStims = numel(sStimObject);
			
			[boolArray,dblCritVal] = DP_GetUpDown(vecPupilSync,0.01,0.99);
			vecSyncLum = vecPupilSync;%double(boolArray);
			
			%test
			
			%
			dblOffsetT = -median(diff(vecPupilStimOn))/4;
			%[vecRefT,matTracePerTrial] = getTraceInTrial(vecPupilTime,vecSyncLum,vecPupilStimOn,1/dblFs,median(diff(vecPupilStimOn)));
			[vecRefT,matTracePerTrial] = getTraceInTrial(vecPupilTime,vecSyncLum,vecPupilStimOn+dblOffsetT,1/dblFs,median(diff(vecPupilStimOn)));
			vecRefT = vecRefT+dblOffsetT;
			figure
			subplot(2,3,1)
			h=plot(vecRefT,matTracePerTrial(1:10:end,:));
			colororder(gca,redbluepurple(numel(h)))
				subplot(2,3,2)
			imagesc(vecRefT,[],matTracePerTrial)
			title(sRec.name,'interpreter','none')
			
			%{
			vecBlinkFractionPerTrial = true(size(vecStimOn));
			for intTrial=1:numel(vecBlinkFractionPerTrial)
				vecStimOn(intTrial)
				vecStimOff(intTrial)
				
			end
			%}
			%% get PD
			%load
			sMetaNI = sRec.sSources.sMeta;
			dblSampRateReportedNI = DP_SampRate(sMetaNI);
			intFirstSample = str2double(sMetaNI.firstSample);
			dblT0_NI = intFirstSample/dblSampRateReportedNI;
			
			strPathNidq = sRec.sSources.sEphysNidq.folder;
			strFileNidq = sRec.sSources.sEphysNidq.name;
			fprintf('   Loading raw data ... [%s]\n',getTime);
			matDataNI = -DP_ReadBin(-inf, inf, sMetaNI, strrep(strFileNidq,'.meta','.bin'), strPathNidq); %1=PD,2=sync pulse
			fprintf('   Done! [%s]\n',getTime);
			
			sNiCh = PP_GetNiCh(sRec.sSources.sMetaVar,sMetaNI);
			vecDiodeSignal = matDataNI(sNiCh.intStimOnsetCh,:);
			dblSubsampleBy = 100;
			vecSubSampledV = double(vecDiodeSignal(1:dblSubsampleBy:end));
			vecT = (1:numel(vecDiodeSignal))/dblSampRateReportedNI;
			vecSubSampledT = double(vecT(1:dblSubsampleBy:end));
			[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecSubSampledT,vecSubSampledV,sBlock.vecStimOnTime+dblOffsetT,(1/dblSampRateReportedNI)*dblSubsampleBy,median(diff(sBlock.vecStimOnTime)));
			vecRefT2 = vecRefT2+dblOffsetT;
			
			%plot
			subplot(2,3,4)
			h=plot(vecRefT2,matTracePerTrial2(1:10:end,:));
			colororder(gca,redbluepurple(numel(h)))
				
			subplot(2,3,5)
			imagesc(vecRefT2,[],matTracePerTrial2)
			
			%% plot spikes
			%%{
			%cellCortex = {'posteromedial visual area','primary visual'};
			%indCtx = contains({sRec.sCluster.Area},cellCortex,'IgnoreCase',true);
			%sNeuronsCtx = sRec.sCluster(indCtx);
			vecStimOn = sBlock.vecStimOnTime;
			vecStimOff = sBlock.vecStimOnTime;
			
			vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
			[C,vecN] = findmin(vecZetaP,1);
			for intN=vecN
				if ~any(cell2vec({sRec.sCluster.ZetaP})<0.05),continue;end
			vecAllSpikes = cell2vec({sRec.sCluster(intN).SpikeTimes});
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecAllSpikes,vecStimOn-0.5,2);
			vecBins = 0:0.01:2;
			intTrials = numel(vecStimOn);
			matSpikes = zeros(intTrials,numel(vecBins)-1);
			for intTrial = 1:intTrials
				matSpikes(intTrial,:)=histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
			end
			subplot(2,3,6)
			imagesc(vecBins(2:end)-0.5,[],matSpikes+1);
			title(sprintf('neuron %d',intN))
			end
			
			%% go through patches and concatenate repetitions
			%cellCellSpikeT = cell(size(sStimObject(1).LinLoc)); %{X,Y}{Rep} = [T_s1 T_s2 ... T_sN]
			
			%% go through stims
			%{
			matRespOn = repmat(nan(size(sStimObject(1).LinLoc)),[1 1 intStims]);
			matRespOff = repmat(nan(size(sStimObject(1).LinLoc)),[1 1 intStims]); %#ok<*REPMAT>
			
			for intStim=1:intStims
				%get spikes this stim
				dblRate = [];
				
				%add to on matrix & off matrix
				matThisOn = nan(size(sStimObject(1).LinLoc));
				matThisOff = nan(size(sStimObject(1).LinLoc));
				matThisOn(sStimObject.LinLocOn) = dblRate;
				matThisOff(sStimObject.LinLocOff) = dblRate;
				matRespOn(:,:,intStim) = matThisOn;
				matRespOff(:,:,intStim) = matThisOff;
			end
			matRespOnMu = mean(matRespOn,3);
			matRespOnSd = std(matRespOn,[],3);
			matRespOffMu = mean(matRespOff,3);
			matRespOffSd = std(matRespOff,[],3);
			%}
			return
		end
	end
	
end
