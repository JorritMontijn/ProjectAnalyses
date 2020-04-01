%% load spiking data & plot tuning curves

%% set recording
strMouse = 'MB2';
strDate = '20190508';
intBlock = 4;
strBlock = num2str(intBlock);
strDataRoot = 'D:\Data\Raw\ePhys';
strStimLog = [strDataRoot filesep 'StimLogs' filesep strMouse '_' strDate];

% set paths
ops.root = [strDataRoot filesep 'KiloSortBinaries']; % 'openEphys' only: where raw files are
ops.rec  = [strMouse '_' strDate '_B' strBlock]; %which recording to process

%% load stimulus info
%load logging file
strPathLogs = strcat(strDataRoot,filesep,'StimLogs');
strSubDir = [strPathLogs filesep strMouse '_' strDate filesep];
sFiles = dir([strSubDir strDate '*_B' strBlock '_*.mat']);
if numel(sFiles) == 1
	sLog = load([strSubDir sFiles(1).name]);
end
structEP = sLog.structEP;

%load triggers

%% TDT data
%get meta data & trigger times
sMetaData = struct;
sMetaData.Mytank = [strDataRoot filesep 'DataTanksTDT' filesep strMouse '_' strDate];
sMetaData.Myblock = ['Block-' strBlock];
sMetaData = getMetaDataTDT(sMetaData);
vecStimOnTime = sMetaData.Trials.stim_onset;
matWord = sMetaData.Trials.word;
[vecStimOnTime,matWord] = checkTriggersTDT(vecStimOnTime,matWord);
vecTrialStartTime = matWord(:,1);
if any(vecTrialStartTime > vecStimOnTime),vecTrialStartTime = vecStimOnTime - 0.5;end
vecStimType = matWord(:,2);
if isfield(sMetaData.Trials,'stim_offset')
	vecStimOffTime = checkTriggersTDT(sMetaData.Trials.stim_offset,matWord);
elseif isfield(sMetaData.Trials,'target_onset')
	vecStimOffTime = checkTriggersTDT(sMetaData.Trials.target_onset,matWord);
else
	vecStimOffTime = vecStimOnTime + 0.5; %use 500 ms as default duration
end

%% load clustered data into matlab using https://github.com/cortex-lab/spikes
% load some of the useful pieces of information from the kilosort and manual sorting results into a struct
strLoadDir = [ops.root filesep ops.rec];
sSpikes = loadKSdir(strLoadDir);

%% load the information from the cluster_groups.csv file with cluster labels
% cids is length nClusters, the cluster ID numbers
% cgs is length nClusters, the "cluster group":
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
[vecClusterIdx, vecClusterType] = readClusterGroupsCSV([strLoadDir filesep 'cluster_groups.csv']);

%% transform data to spikes per cluster
vecAllSpikeTimes = sSpikes.st;
vecAllSpikeClust = sSpikes.clu;
vecSingleUnits = vecClusterIdx(vecClusterType==2);
intNumSU = numel(vecSingleUnits);
vecMultiUnits = vecClusterIdx(vecClusterType==1);
intNumMU = numel(vecMultiUnits);
%assign single unit spikes
SU_st = cell(1,intNumSU); %single unit spike times
for intClustSUA=1:intNumSU
	SU_st{intClustSUA} = vecAllSpikeTimes(vecAllSpikeClust==vecSingleUnits(intClustSUA));
end
%assign multi unit spikes
MU_st = cell(1,intNumMU); %multi unit spike times
for intClustMUA=1:intNumMU
	MU_st{intClustMUA} = vecAllSpikeTimes(vecAllSpikeClust==vecMultiUnits(intClustMUA));
end

%% transform to spikes to response matrix
vecStimOnTime = vecStimOnTime + 0.05;
vecStimOffTime = vecStimOnTime + 0.2;
vecTrialStartTime = vecStimOnTime - 0.2;
matStimCounts = getSpikeCounts(SU_st,vecStimOnTime,vecStimOffTime);
matStimResp = bsxfun(@rdivide,matStimCounts,(vecStimOffTime-vecStimOnTime)'); %transform to Hz
matBaseCounts = getSpikeCounts(SU_st,vecTrialStartTime,vecStimOnTime);
matBaseResp = bsxfun(@rdivide,matBaseCounts,(vecStimOnTime-vecTrialStartTime)'); %transform to Hz
matResp = matStimResp - matBaseResp;

%% get tuning curves
vecStimOriDegrees = structEP.Orientation;
[vecOriIdx,vecOriTypes,vecReps] = label2idx(vecStimOriDegrees);
vecAngles = deg2rad(vecStimOriDegrees);
[sOut] = getTuningCurves(matResp,vecStimOriDegrees);
vecDeltaPrime = getDeltaPrime(matResp,vecAngles);

%% make PSTH
%gather data
sOptions.handleFig = 1;
for intSU = 1:intNumSU
	vecSpikeTimes = SU_st{intSU};
	vecWindow = -0.5:0.05:1.5;
	figure
	for intStimType=1:numel(vecOriTypes)
		vecEvents = vecStimOnTime(vecOriIdx==intStimType);
		subplot(3,3,intStimType)
		[vecMean,vecSEM] = doPEP(vecSpikeTimes,vecWindow,vecEvents,sOptions);
		xlabel('Time from stim on (s)');
		ylabel('Spiking rate (Hz)');
		title(sprintf('Single Unit %d; ori %d (%.1f degs)',intSU,intStimType,vecOriTypes(intStimType)));
	end
end
