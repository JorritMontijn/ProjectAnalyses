%% find data
close all;
clear all;
strDir = 'F:\Drive\EyeTrackingProcessed';
if ~strcmp(strDir(end),filesep),strDir(end+1) = filesep;end
%strRec = 'EyeTrackingProcessed2021-02-12_R.mat';
sFiles = dir([strDir '*Processed*2021*.mat']);
sFiles = cat(1,sFiles,dir([strDir '*2021*Processed*.mat']));
%intRec = 3;
for intRec=1:numel(sFiles)
	strRec = sFiles(intRec).name;
	strSearchFormat = '\d{4}[-_]?\d{2}[-_]?\d{2}';
	[intB,intE]=regexp(strRec,strSearchFormat);
	strRecFmt0 = strRec(intB:intE);
	strRecFmt1 = strrep(strrep(strRecFmt0,'-',''),'_','');
	strRecFmt2 = strcat(strRecFmt1(1:4),'-',strRecFmt1(5:6),'-',strRecFmt1(7:8));
	%get files
	strTrackFile = strRec;
	strDir = sFiles(intRec).folder;
	if ~strcmp(strDir(end),filesep),strDir(end+1) = filesep;end

	sStimFiles1 = dir([strDir '*'  strRecFmt1 '*.mat']);
	sStimFiles2 = dir([strDir '*'  strRecFmt2 '*.mat']);
	sStimFiles = cat(1,sStimFiles1,sStimFiles2);
	sStimFiles(contains({sStimFiles.name},'Processed')) = [];
	
	%% load data
	%load eye-tracking data
	sLoad = load([strDir strRec]);
	sPupil = sLoad.sPupil;
	vecFrTimePupil = sPupil.sSyncData.matSyncData(2,:);
	[dblStartF,intIdx] = min(vecFrTimePupil);
	vecPupilTime = sPupil.vecPupilTime; %time corresponding to data with same clock as vecVidTimePupil
	dblFR = 1/median(diff(vecPupilTime));
	dblOffsetT = dblStartF*dblFR;
	
	vecVidTimePupil = sPupil.sSyncData.matSyncData(1,1:(end-1)); %timestamps of video
	vecVidTimePupil = vecVidTimePupil - (vecVidTimePupil(intIdx) - dblStartF*dblFR);
	vecNITimePupil = sPupil.sSyncData.matSyncData(3,1:(end-1)); %synchronous timestamps of NI time
	
	
	%load stims
	cellStim = {};
	intC = 0;
	for intStimFile = 1:numel(sStimFiles)
		sLoad = load([strDir sStimFiles(intStimFile).name]);
		sEP = sLoad.structEP;
		%get data
		vecOnNI = sEP.ActOnNI;
		if all(isnan(vecOnNI))
			error([mfilename ':AllNaN'],sprintf('Data is all nan, please remove %s',sStimFiles(intStimFile).name));
		else
			vecOffNI = sEP.ActOffNI;
			sStims=sEP.sStimObject(sEP.ActStimType(~isnan(sEP.ActStimType)));
			intC = intC + 1;
			cellStim{intC}.structEP = sEP;
		end
	end
	
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
	
	%% find starting times from NI timestamps
	vecStartT = zeros(size(cellStim));
	for intLog=1:numel(cellStim)
		vecActOnNI=cellStim{intLog}.structEP.ActOnNI;
		vecPupilT=sPupil.sSyncData.matSyncData(1,:) - sPupil.sSyncData.matSyncData(1,1) + sPupil.vecPupilTime(1);
		vecPupilNI=sPupil.sSyncData.matSyncData(3,:);
		intStart = find(vecPupilNI>vecActOnNI(1),1)-1;
		vecStartT(intLog) = vecPupilT(intStart);
	end
	%% synchronize pupil time to stim time
	cellRecNames = {sStimFiles(:).name};
	cellStim = SyncET_NoEphysOld(sPupil,cellStim,cellRecNames,vecStartT);
	
	%% run analysis
	getPupilPlotsNoEphys;%(sPupil,cellStim,cellRecNames,sFiles)
end