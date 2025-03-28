%% stimulus descriptions
%{
%http://observatory.brain-map.org/visualcoding/stimulus/natural_movies#
%help.brain-map.org/download/attachments/10616846/VisualCoding_VisualStimuli.pdf

Drifting Gratings
Stimulus. This stimulus was used to measure the direction tuning, orientation tuning and temporal frequency
tuning of the cells. The total stimulus duration was 31.5 minutes.
The stimulus consisted of a full field drifting sinusoidal grating at a single spatial frequency (0.04 cycles/degree)
and contrast (80%). The grating was presented at 8 different directions (separated by 45?) and at 5 temporal
frequencies (1, 2, 4, 8, 15 Hz). Each grating was presented for 2 seconds, followed by 1 second of mean
luminance gray before the next grating. Each grating condition (direction & temporal frequency combination)
was presented 15 times, in a random order. There were blank sweeps (i.e. mean luminance gray instead of
grating) presented approximately once every 20 gratings.

Natural Scenes
Stimulus. Receptive fields measured using synthetic stimuli, such as gratings and noise, often fail to predict
responses to natural stimuli (David et al., 2004). This stimulus allows us to measure cells? responses to natural
scenes to explore non-linear processing in their receptive fields.
The stimulus consisted of 118 natural images. Images 1-58 were from the Berkeley Segmentation Dataset
(Martin et al., 2001), images 59-101 from the van Hateren Natural Image Dataset (van Hateren and van der
Schaaf, 1998), and images 102-118 are from the McGill Calibrated Colour Image Database (Olmos and
Kingdom, 2004). The images were presented in grayscale and were contrast normalized and resized to 1174 x
918 pixels. Images were luminance-matched, but a bug in our code resulted in an error in this matching such
that the mean luminance of the images were not accurately matched, and were within 6% of each other.
The images were presented for 0.25 seconds each, with no inter-image gray period. Each image was presented
~50 times, in random order, and there were blank sweeps (i.e. mean luminance gray instead of an image)
roughly once every 100 images.

Natural Movies
Stimulus. This stimulus was used to measure responses to natural movies.
Three different clips were used from the opening scene of the movie Touch of Evil (Welles, 1958). This scene
was chosen as it is a long take with a range of different features and scales of motion. Natural Movie 1 and
Natural Movie 2 were both 30 second clips while Natural Movie 3 was a 120 second clip. All clips had been
contrast normalized and were presented in grayscale at 30 fps. Each movie was presented 10 times in a row
with no inter-trial gray period.
%}

%% define locations
strDisk = 'F:';
strDataSource = '\Data\Processed\AllenBrainVisualEphys\nwb_files\visual-behavior-neuropixels-0.5.0\';
strDataTarget = strcat(strDisk,strDataSource,'Aggregates',filesep);
%load metadata
sCSV_behav = loadcsv(strcat(strDisk,strDataSource,'project_metadata\behavior_sessions.csv'));
sCSV_ephys = loadcsv(strcat(strDisk,strDataSource,'project_metadata\ecephys_sessions.csv'));
%filter sessions
boolKeep = contains(sCSV_ephys.structure_acronyms,'VISp') & cellfun(@strcmp,sCSV_ephys.abnormal_histology,cellfill(',',size(sCSV_ephys.abnormal_histology)));
vecSessions = sCSV_ephys.ecephys_session_id(boolKeep);
%vecSessions = 756029989;
strFormat = 'ecephys_session_%d';

%% loop through sessions
for intSes=1:numel(vecSessions)
	
	%% load data
	strDataFile = sprintf(strFormat,vecSessions(intSes));
	strDataPath = sprintf('%s%s%s%d%s',strDataSource,'behavior_ecephys_sessions',filesep,vecSessions(intSes),filesep);
	strTargetFile = strcat(strDisk,fullpath(strDataPath,[strDataFile,'.nwb']));
	fclose('all');
	fprintf('Processing session %d/%d: %s [%s]\n',intSes,numel(vecSessions),strDataFile,getTime);
	nwb = nwbRead(strTargetFile);
	sNWB = ExpandNWB(nwb,[],false);
	%sNWB2 = PruneStruct(sNWB);
	
	%% get spike times & unit locations
	unit_ids = sNWB.units.id; % array of unit ids represented within this
	spike_times = sNWB.units.spike_times; % array of unit ids represented within this
	spike_cluster_ids = sNWB.units.spike_times_index; % array of unit ids represented within this
	
	%get data
	vecElectrodeIDs = sNWB.general_extracellular_ephys_electrodes.id;
	cellElectrodeLocations = sNWB.general_extracellular_ephys_electrodes.vectordata.location;
	vecUnitPeakElectrodeIDs = sNWB.units.vectordata.peak_channel_id;
	
	%loop
	intUnits = numel(unit_ids);
	cellSpikes = cell(1,intUnits);
	cellAreas = cell(1,intUnits);
	intLastEnd = 0;
	for intUnit=1:intUnits
		cellAreas{intUnit} = cellElectrodeLocations{vecElectrodeIDs==vecUnitPeakElectrodeIDs(intUnit)};
		cellSpikes{intUnit} = spike_times((intLastEnd+1):spike_cluster_ids(intUnit));
		intLastEnd = spike_cluster_ids(intUnit);
	end
	
	%% get cluster properties
	% cut-offs
	dblThresh_isi_violations = 0.5;% < 0.5
	dblThresh_amplitude_cutoff = 0.1;%< 0.1
	dblThresh_presence_ratio = 0.9;%> 0.9
	
	%"Set" can use: .keys .values .Count .size .set .remove .clear .get .export
	vecCluster_ID=sNWB.units.vectordata.cluster_id;
	indCluster_Good=strcmp(sNWB.units.vectordata.quality,'good');
	vecCluster_Viol=sNWB.units.vectordata.isi_violations;
	vecCluster_AmpC=sNWB.units.vectordata.amplitude_cutoff;
	vecCluster_Pres=sNWB.units.vectordata.presence_ratio;
	
	indInclude = indCluster_Good & ...
		(vecCluster_Viol < dblThresh_isi_violations) & ...
		(vecCluster_AmpC < dblThresh_amplitude_cutoff) & ...
		(vecCluster_Pres > dblThresh_presence_ratio);
	
	%create neuron structure
	vecUseCells = find(indInclude & contains(cellAreas(:),{'VIS','LP','LG'}));
	intNeurons = numel(vecUseCells);
	
	matWaveforms = sNWB.units.waveform_mean(:,sNWB.units.vectordata.waveform_mean_index);
	clear sNeuron
	cellUnitFields = fieldnames(sNWB.units.vectordata);
	for i=1:numel(vecUseCells)
		intNeuron = vecUseCells(i);
		
		sUnit = struct;
		sUnit.OrigIdx = intNeuron;
		sUnit.Area = cellAreas{intNeuron};
		sUnit.SpikeTimes = cellSpikes{intNeuron};
		sUnit.Waveform = matWaveforms(:,intNeuron);
		
		for j=1:numel(cellUnitFields)
			strField = cellUnitFields{j};
			sUnit.(strField) = sNWB.units.vectordata.(strField)(i);
		end
		sNeuron(i) = sUnit;
	end
	
	%% subject
	sSubject = sNWB.general_subject;
	structStimAgg = struct;
	
	%% stimuli, trials
	if isfield(sNWB,'intervals_trials')
		sTrials = sNWB.intervals_trials.vectordata;
		sTrials.start_time = sNWB.intervals_trials.start_time;
		sTrials.stop_time = sNWB.intervals_trials.stop_time;
		sTrials.id = sNWB.intervals_trials.id;
		structStimAgg.sTrials=sTrials;
	end
	
	%% stimuli, DG
	vecAllStimTimes = sNWB.processing.stimulus.nwbdatainterface.timestamps.data;
	if isfield(sNWB.intervals,'drifting_gratings_presentations')
		strSource = 'drifting_gratings_presentations';
	elseif isfield(sNWB.intervals,'drifting_gratings_75_repeats_presentations')
		strSource = 'drifting_gratings_75_repeats_presentations';
	else
		%no gratings
		strSource = '';
	end
	if ~isempty(strSource)
		sSource = sNWB.intervals.(strSource);
		sDG = sSource.vectordata;%orientation, SF, TF, etc
		sDG.Name = strSource;
		sDG.stimIdx = table2array(sSource.timeseries(:,1))+1;
		sDG.startT = sSource.start_time;
		sDG.stopT = sSource.stop_time;
		%vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
		structStimAgg.sDG = sDG;
		vecGood = sDG.orientation(~isnan(sDG.orientation));
		intStims = numel(unique(vecGood));
		intReps = numel(vecGood)/intStims;
		fprintf('Ses %d (%s) has %.1f repetitions of %d drifting gratings\n',intSes,sprintf(strFormat,vecSessions(intSes)),intReps,intStims);
	end
	
	%% stimuli, nat scenes
	cellNatScenes = {'natural_scenes_presentations','Natural_Images_Lum_Matched'};
	for i=1:numel(cellNatScenes)
		strCheckVar = cellNatScenes{i};
		cellFields = fieldnames(sNWB.intervals);
		intField = find(contains(cellFields,strCheckVar));
		if isempty(intField),continue;end
		strField = cellFields{intField};
		sNS = sNWB.intervals.(strField).vectordata;%orientation, SF, TF, etc
		sNS.start_time = sNWB.intervals.(strField).start_time;%orientation, SF, TF, etc
		sNS.stop_time = sNWB.intervals.(strField).stop_time;%orientation, SF, TF, etc
		%remove non-scenes
		if isfield(sNS,'frame')
			indBlanks = sNS.frame==-1;
			vecRealScenes = sNS.frame(~indBlanks);
		else
			indBlanks = sNS.omitted==1;
			vecRealScenes = val2idx(sNS.image_name(~indBlanks));
		end
		intScenes = numel(unique(vecRealScenes));
		vecNatSceneReps(intSes) = numel(vecRealScenes)/intScenes;
		fprintf('Ses %d (%s) has %.1f repetitions of %d unique natural scenes\n',intSes,sprintf(strFormat,vecSessions(intSes)),vecNatSceneReps(intSes),intScenes)
		structStimAgg.sNS = sNS;
	end
	
	%% stimuli, nat movies
	if isfield(sNWB.intervals,'natural_movie_one_presentations')
		sNM = sNWB.intervals.natural_movie_one_presentations.vectordata;%orientation, SF, TF, etc
		%remove non-scenes; 30fps, 30s long(?)
		indBlanks = sNM.frame==-1;
		vecMovieFrames = sNM.frame(~indBlanks);
		intReps = sum(diff(vecMovieFrames)<0)+1;
		vecNatSceneReps(intSes) = numel(vecRealScenes)/intScenes;
		fprintf('Ses %d (%s) has %d repetitions of natural movie one\n',intSes,sprintf(strFormat,vecSessions(intSes)),intReps);
	
		sNM.Name = 'natural_movie_one_presentations';
		sNM.stimIdx = sNM.frame+1;
		sNM.startT = sSource.start_time;
		sNM.stopT = sSource.stop_time;
		%vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
		structStimAgg.sNM = sNM;
	end
	
	%% remove trials during invalid times
	sRemoveIntervals = sNWB.intervals_invalid_times;
	
	%% eye-tracking
	%not in all v1 sessions!
	vecPupilSize = [];
	vecPupilT = [];
	vecPupilBlink = [];
	matPupilXY = [];
	if isfield(sNWB.processing,'filtered_gaze_mapping')
		%% v1
		vecPupilT = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.pupil_area.timestamps; %sec
		vecPupilSize = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.pupil_area.data; %pixel^2
		matPupilXY = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.screen_coordinates.data; %cm
	elseif isfield(sNWB.acquisition,'EyeTracking')
		%% v2
		vecPupilT = sNWB.acquisition.EyeTracking.eye_tracking.timestamps; %sec
		vecPupilSize = sNWB.acquisition.EyeTracking.pupil_tracking.area; %pixel^2
		matPupilXY = sNWB.acquisition.EyeTracking.pupil_tracking.data; %m
		vecPupilBlink = sNWB.acquisition.EyeTracking.likely_blink.data;
	end
	
	%% running
	if isfield(sNWB.processing.running.nwbdatainterface,'running_speed')
		%% running, v1
		vecRunningSpeed = sNWB.processing.running.nwbdatainterface.running_speed.data; %speed (cm/s)
		vecRunningTime = sNWB.processing.running.nwbdatainterface.running_speed.timestamps; %time (s)
	else
		%% behavior, v2
		vecRunningSpeed = sNWB.processing.running.nwbdatainterface.speed.data; %speed (cm/s)
		vecRunningTime = sNWB.processing.running.nwbdatainterface.speed.timestamps; %time (s)
	end
	
	%% licking & rewards
	%already part of sTrials
% 	vecLickT = [];
% 	try
% 		vecLickT = sNWB.processing.licking.nwbdatainterface.licks.timestamps;
% 	end
% 	vecRewardT = [];
% 	try
% 		vecRewardT = sNWB.processing.rewards.nwbdatainterface.autorewarded.timestamps;
% 	end
% 	
% 	%add to struct stim
% 	structStimAgg.vecLickT = vecLickT;
% 	structStimAgg.vecRewardT = vecRewardT;
% 	
	%% optotagging
	vecOptoT = [];
	vecOptoVal = [];
	try
		vecOptoT = sNWB.processing.optotagging.nwbdatainterface.optotagging.timestamps;
		vecOptoVal = sNWB.processing.optotagging.nwbdatainterface.optotagging.data;
	end
	
	%% get optogenetic modulation
	
	%gather data
	sSes = struct;
	sSes.SesId = intSes;
	sSes.Exp = strDataFile;
	sSes.vecUseCells = vecUseCells;
	sSes.cellAreas = cellAreas(vecUseCells); %redundant, but good as a sanity check
	
	%units
	sSes.sNeuron = sNeuron;
	
	%stims & licking
	sSes.structStimAgg = structStimAgg;
	
	%other behavior
	sSes.vecPupilT = vecPupilT;
	sSes.vecPupilBlink = vecPupilBlink;
	sSes.vecPupilSize = vecPupilSize;
	sSes.matPupilXY = matPupilXY;
	sSes.vecRunningTime = vecRunningTime;
	sSes.vecRunningSpeed = vecRunningSpeed;
	
	%opto
	sSes.vecOptoT = vecOptoT;
	sSes.vecOptoVal = vecOptoVal;
	
	%to-be-removed-interval structure
	sSes.sRemoveIntervals = sRemoveIntervals;
	
	%% save data
	strOutputFile = strcat(strDataTarget,'AggSes_',strDataFile,'_ProcDate',getDate,'.mat');
	fprintf('Saving data to "%s"\n',strOutputFile);
	save(strOutputFile,'sSes','-v7.3');
	fprintf('Done.\n');
end
