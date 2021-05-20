%% find data
close all;
clear all;
strDir = 'F:\Drive\EyeTrackingProcessed';
if ~strcmp(strDir(end),filesep),strDir(end+1) = filesep;end
%strRec = 'EyeTrackingProcessed2021-02-12_R.mat';
sFiles = dir([strDir '*Processed*2019*.mat']);
sFiles = cat(1,sFiles,dir([strDir '*2019*Processed*.mat']));
sFiles = cat(1,sFiles,dir([strDir '*2020*Processed*.mat']));
sFiles = cat(1,sFiles,dir([strDir '*Processed*2020*.mat']));
%intRec = 3;


for intRec=1:numel(sFiles)
	%% find stim files
	strRec = sFiles(intRec).name;
	strSearchFormat = '\d{4}[-_]?\d{2}[-_]?\d{2}';
	[intB,intE]=regexp(strRec,strSearchFormat);
	strRecFmt0 = strRec(intB:intE);
	strRecFmt1 = strrep(strrep(strRecFmt0,'-',''),'_','');
	strRecFmt2 = strcat(strRecFmt1(1:4),'-',strRecFmt1(5:6),'-',strRecFmt1(7:8));
	%get files
	strTrackFile = strRec;
	sStimFiles = dir([strDir '*'  strRecFmt1 '*.mat']);
	%get preprocessed data
	strDataDir = 'F:\Data\Processed\Neuropixels\';
	sDataFiles = dir([strDataDir '*'  strRecFmt2 '*.mat']);
	if numel(sDataFiles) ~= 1
		error help
	end
	sLoad = load(fullfile(sDataFiles(1).folder,sDataFiles(1).name));
	cellStimOld = sLoad.sAP.cellStim;
	
	% determine temporal order
	cellFiles = {sStimFiles(:).name};
	intLogs = numel(sStimFiles);
	vecTimes = nan(1,intLogs);
	for intLogFile = 1:intLogs
		cellSplit = strsplit(cellFiles{intLogFile}(1:(end-4)),'_');
		vecTimes(intLogFile) = str2double(cat(2,cellSplit{end-2:end}));
	end
	[dummy,vecReorderStimFiles] = sort(vecTimes);
	sStimFiles = sStimFiles(vecReorderStimFiles);
	
	%% prepro pupil
	%load eye-tracking data
	sLoad = load([strDir strRec]);
	sPupil = sLoad.sPupil;
	
	%filter sync lum
	vecPupilSyncTime = sPupil.vecPupilFullSyncLumT;
	dblSampRatePupil = 1/median(diff(vecPupilSyncTime));
	
	%filter to 0.1-30Hz
	vecWindow2 = [0.01 30]./(dblSampRatePupil./2);
	[fb,fa] = butter(2,vecWindow2,'bandpass');
	vecFiltSyncLum = filtfilt(fb,fa, double(sPupil.vecPupilFullSyncLum));
	boolPupilSync1 = vecFiltSyncLum>(-std(vecFiltSyncLum)/2);
	boolPupilSync2 = vecFiltSyncLum>(std(vecFiltSyncLum)/3);
	
	%get on/off
	[boolPupilSync0,dblCritValStim]=DP_GetUpDown(vecFiltSyncLum,0.5,0.9);
	vecChangePupilSync0 = diff(boolPupilSync0);
	vecChangePupilSync1 = diff(boolPupilSync1);
	vecChangePupilSync2 = diff(boolPupilSync2);
	vecPupilSyncOn = [(find(vecChangePupilSync0 == 1 | vecChangePupilSync1 == 1 | vecChangePupilSync2 == 1)+1)];
	
	%find transitions
	vecPupilStimOn = vecPupilSyncTime(vecPupilSyncOn);
	sPupil.vecPupilStimOn = vecPupilStimOn;
	sPupil.vecPupilSyncLumFilt = vecFiltSyncLum;
	
	%% load stim data
	cellStim = {};
	intC = 0;
	for intStimFile = 1:intLogs
		sLoad = load([strDir sStimFiles(intStimFile).name]);
		sEP = sLoad.structEP;
		%get data
		vecOnSecs = sEP.ActOnSecs;
		if all(isnan(vecOnSecs))
			continue;
		else
			sStims=sEP.sStimObject(sEP.ActStimType(~isnan(sEP.ActStimType)));
			intC = intC + 1;
			cellStim{intC}.structEP = sEP;
		end
	end
	
	%% synchronize pupil time to stim time
	cellRecNames = {sStimFiles(:).name};
	vecStartT = zeros(size(cellRecNames));
	for intRec=1:numel(cellRecNames)
		if isfield(cellStimOld{intRec}.structEP,'vecPupilStimOnTime')
			vecStartT(intRec) = cellStimOld{intRec}.structEP.vecPupilStimOnTime(1);
		else
			vecStartT(intRec) = nan;
		end
	%vecStartT = cellfun(@(x) x.structEP.vecPupilStimOnTime(1),cellStimOld);
	end
	cellStim = SyncET_NoEphysOld(sPupil,cellStim,cellRecNames,vecStartT);
	
	%% run analysis
	getPupilPlotsNoEphys;%(sPupil,cellStim,cellRecNames,sFiles)
end