%% message
fprintf('Loading raw data for %s_%s, B%s [%s]\n',strMouse,strDate,strBlock,getTime);

%% load stimulus info
%check if recording is single or multi-block
cellBlocks=strsplit(strBlock,'-');
structEP = struct;
dblLastTime = 0;
intStimTypeCounter = 0;
%pre-allocate time vectors
vecStimOnTime = [];
vecStimOffTime = [];
vecTrialStartTime = [];
vecStimType = [];
vecEyeTimestamps = [];
matEyeData = [];
vecEnvTimestamps = [];
matEnvData = [];
vecRawTimestamps = [];
vecLFPTimestamps = [];
matLFPData = [];
	
boolLoadOnlyDG = true;

for intBlock=1:numel(cellBlocks)
	%% load logging file
	strThisBlock = cellBlocks{intBlock};
	strPathLogs = strcat(strDataRoot,filesep,'StimLogs');
	strSubDir = [strPathLogs filesep strMouse '_' strDate filesep];
	sFiles = dir([strSubDir strDate '*_B' strThisBlock '_*.mat']);
	if numel(sFiles) == 1
		sLog = load([strSubDir sFiles(1).name]);
	end
	
	%% build structEP
	if dblLastTime == 0
		%assign
		structEP = sLog.structEP;
		%set log for which trials belong to which stim script; 
		%{'Script',[StartTrial EndTrial],(# of stim types)}
		structEP.cellStimBlock = {structEP.strFile,[1 numel(structEP.ActStartSecs)],numel(unique(structEP.ActStimType))};
		
		%concatenate sStimObject as cell array
		intObjNr = numel(structEP.sStimObject);
		cellStimObject = cell(1,intObjNr);
		for intObject=1:intObjNr
			cellStimObject(intObject) = {structEP.sStimObject(intObject)};
		end
		structEP.cellStimObject = cellStimObject;
	else
		%get data
		structEP_Old = structEP;
		structEP_New = sLog.structEP;
		
		%check if it's the same stimulus type
		if strcmpi(structEP_New.strFile,structEP_Old.strFile)
			boolUseSame = true;
		else
			boolUseSame = false;
		end
		
		%concatenate fields
		cellNewData = {structEP_New.strFile,[1 numel(structEP_New.ActStartSecs)],numel(unique(structEP_New.ActStimType))};
		structEP.cellStimBlock = cat(1,structEP_Old.cellStimBlock,cellNewData);
		cellFieldsOld = fieldnames(structEP);
		cellFieldsNew = fieldnames(structEP_New);
		cellFields = unique(cat(1,cellFieldsOld,cellFieldsNew));
		for intField=1:numel(cellFields)
			strField = cellFields{intField};
			if strcmpi(strField,'cellStimBlock')
				%ignore
			elseif strcmpi(strField,'sStimObject')
				%concatenate sStimObject as cell array
				intObjNr = numel(structEP_New.(strField));
				cellNewStimObject = cell(1,intObjNr);
				for intObject=1:intObjNr
					cellNewStimObject(intObject) = {structEP_New.(strField)(intObject)};
				end
				structEP.cellStimObject = cat(2,structEP_Old.cellStimObject,cellNewStimObject);
			elseif ~isfield(structEP_New,strField) && numel(structEP_Old.(strField)) > 1
				%add nans
				structEP.(strField) = cat(2,structEP_Old.(strField), nan(size(structEP_New.ActStimType)));
			elseif ~isfield(structEP_Old,strField) && numel(structEP_New.(strField)) > 1
				%add retroactive nans
				structEP.(strField) = cat(2,nan(size(structEP_Old.ActStimType)), structEP_New.(strField));
			elseif isfield(structEP_New,strField) && isfield(structEP_Old,strField) && numel(structEP_Old.(strField)) > 1
				if ~isempty(strfind(strField,'StimType')) %#ok<STREMP>
					if boolUseSame || (boolLoadOnlyDG && isfield(structEP_New,'Orientation'))
						%do not add previous # of stim types
						structEP.(strField) = cat(2,structEP_Old.(strField), structEP_New.(strField));
					else
						%add previous # of stim types
						structEP.(strField) = cat(2,structEP_Old.(strField), structEP_New.(strField) + max(structEP_Old.(strField)));
					end
				else
					if ~isempty(strfind(strField,'Secs')) %#ok<STREMP>
						%add time offset
						structEP.(strField) = cat(2,structEP_Old.(strField), structEP_New.(strField) + dblLastTime);
					else
						%simple concatenate
						structEP.(strField) = cat(2,structEP_Old.(strField), structEP_New.(strField));
					end
				end
			else
				%scalar, set to nan
				structEP.(strField) = nan;
			end
		end
	end
	
	%% TDT data
	%get meta data & trigger times
	sMetaData = struct;
	sMetaData.Mytank = [strDataRoot filesep 'DataTanksTDT' filesep strMouse '_' strDate];
	sMetaData.Myblock = ['Block-' strThisBlock];
	sMetaData = getMetaDataTDT(sMetaData);
	
	vecThisStimOnTime = sMetaData.Trials.stim_onset;
	matThisWord = sMetaData.Trials.word;
	[vecThisStimOnTime,matThisWord] = checkTriggersTDT(vecThisStimOnTime,matThisWord);
	vecThisTrialStartTime = matThisWord(:,1);
	if any(vecThisTrialStartTime > vecThisStimOnTime),vecThisTrialStartTime = vecThisStimOnTime - 0.5;end
	vecThisStimType = matThisWord(:,2);
	if isfield(sMetaData.Trials,'stim_offset')
		vecThisStimOffTime = checkTriggersTDT(sMetaData.Trials.stim_offset,matThisWord);
	elseif isfield(sMetaData.Trials,'target_onset')
		vecThisStimOffTime = checkTriggersTDT(sMetaData.Trials.target_onset,matThisWord);
	else
		vecThisStimOffTime = vecThisStimOnTime + 0.5; %use 500 ms as default duration
	end
	
	%% get raw 
	strRawEvent = 'dRAW';
	intRaw = ismember({sMetaData.strms.name},strRawEvent);
	sMetaData.Myevent = strRawEvent;
	[vecTimestampsRAW,matBlockData,vecChannels] = getRawDataTDT(sMetaData);
	vecRawTimestamps = cat(2,vecRawTimestamps,vecTimestampsRAW + dblLastTime);
	
	%% build LFP
	%Clean raw data, there is a 1-sample mismatch in top 2 channels
	for intCh = 31:32
		matBlockData(intCh,:) = circshift(matBlockData(intCh,:),[0 -1]);
		matBlockData(intCh,end) = matBlockData(intCh,end-1);
	end
	
	%re-reference odd by average of all odd channels, and even by even
	matBlockData(1:2:end,:) = bsxfun(@minus,matBlockData(1:2:end,:),cast(mean(matBlockData(1:2:end,:),1),'like',matBlockData)); %odd
	matBlockData(2:2:end,:) = bsxfun(@minus,matBlockData(2:2:end,:),cast(mean(matBlockData(2:2:end,:),1),'like',matBlockData)); %even
				
	%resample to 1.5KHz
	dblResampFreq = 1500;
	vecTimestampsLFP = min(vecTimestampsRAW):(1/dblResampFreq):max(vecTimestampsRAW);
	matBlockData = squeeze(get(resample(timeseries(matBlockData,vecTimestampsRAW),vecTimestampsLFP),'Data'));
	matFiltered = nan(size(matBlockData));
	%filter each channel
	for intCh=1:size(matBlockData,1)
		%get data
		vecResampled = matBlockData(intCh,:);
		
		
		%filter 50Hz
		vecWindow = [45 55]./(dblResampFreq./2);
		[fb,fa] = butter(2,vecWindow,'stop');
		vecFiltered = filtfilt(fb,fa,double(vecResampled));
		
		%filter to 1-300Hz
		vecWindow2 = [1 300]./(dblResampFreq./2);
		[fb,fa] = butter(2,vecWindow2,'bandpass');
		matFiltered(intCh,:) = filtfilt(fb,fa,vecFiltered);
		
		%calc power
		%[vecFreq,vecPower] = getPowerSpectrum(matFiltered(intCh,:),dblSampFreq,2);
		%plot(vecFreq,vecPower);
	end
	
	%assign data
	vecLFPTimestamps = cat(2,vecLFPTimestamps,vecTimestampsLFP + dblLastTime);
	matLFPData = cat(2,matLFPData,matFiltered);
	
	%% get envelope
	strEnvEvent = 'dENV';
	intEnv = ismember({sMetaData.strms.name},strEnvEvent);
	sMetaData.Myevent = strEnvEvent;
	[vecTimestampsENV,matBlockData,vecChannels] = getRawDataTDT(sMetaData);
	vecEnvTimestamps = cat(2,vecEnvTimestamps,vecTimestampsENV + dblLastTime);
	matEnvData = cat(2,matEnvData,matBlockData);
	
	%% get eye-tracking and diode
	strEyeEvent = 'Eye_';
	intEyeTracker = ismember({sMetaData.strms.name},strEyeEvent);
	sMetaData.Myevent = strEyeEvent;
	[vecTimestamps,matBlockData,vecChannels] = getRawDataTDT(sMetaData);
	vecEyeTimestamps = cat(2,vecEyeTimestamps,vecTimestamps + dblLastTime);
	matEyeData = cat(2,matEyeData,matBlockData);
	
	%% combine TDT
	vecStimOnTime = cat(1,vecStimOnTime,vecThisStimOnTime + dblLastTime);
	vecStimOffTime = cat(1,vecStimOffTime,vecThisStimOffTime + dblLastTime);
	vecTrialStartTime = cat(1,vecTrialStartTime,vecThisTrialStartTime + dblLastTime);
	vecStimType = cat(1,vecStimType,vecThisStimType + intStimTypeCounter);
	
	%% increment time
	dblLastTime = dblLastTime + vecTimestampsRAW(end);
	intStimTypeCounter = max(vecStimType);
end
%% message
fprintf('Retrieving clustered spiking data... \n');

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

%% get spike depths
[spikeAmps, vecAllSpikeDepth] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps); 
dblRecDepth = sSession(intSelectSession).Recording(intSelectRec).Depth;
vecChannelDepth = sSession(intSelectSession).DepthPerCh + dblRecDepth;
vecChannelDepth = vecChannelDepth(end:-1:1) - 31*25;

%% transform data to spikes per cluster
vecAllSpikeTimes = sSpikes.st;
vecAllSpikeClust = sSpikes.clu;
vecSingleUnitClusters = vecClusterIdx(vecClusterType==2 | vecClusterType==3);
intNumSU = numel(vecSingleUnitClusters);
vecMultiUnitClusters = vecClusterIdx(vecClusterType==1);
intNumMU = numel(vecMultiUnitClusters);
%assign single unit spikes & depth
SU_st = cell(1,intNumSU); %single unit spike times
SU_depth_ch = nan(1,intNumSU); %single unit spike times
SU_depth_micron = nan(1,intNumSU); %single unit spike times
for intClustSUA=1:intNumSU
	intClustIdx = vecSingleUnitClusters(intClustSUA);
	vecSpikeIDs = intClustIdx==vecAllSpikeClust;
	
	SU_depth_ch(intClustSUA) = mean(vecAllSpikeDepth(vecSpikeIDs));
	SU_depth_micron(intClustSUA) = getFractionalEntry(vecChannelDepth,SU_depth_ch(intClustSUA));
	SU_st{intClustSUA} = vecAllSpikeTimes(vecSpikeIDs);
end
%assign multi unit spikes
MU_st = cell(1,intNumMU); %multi unit spike times
for intClustMUA=1:intNumMU
	MU_st{intClustMUA} = vecAllSpikeTimes(vecAllSpikeClust==vecMultiUnitClusters(intClustMUA));
end
dblLastSpike = max(vecAllSpikeTimes);

%% remove trials  that end before last spike
indRemTrials = vecStimOffTime > dblLastSpike;
vecStimOnTime(indRemTrials) = [];
vecStimOffTime(indRemTrials) = [];
vecTrialStartTime(indRemTrials) = [];
vecStimType(indRemTrials) = [];

%% message
fprintf('\b   Done! [%s]\n',getTime);
