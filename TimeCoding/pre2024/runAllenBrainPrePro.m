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
strDataSource = '\Data\Processed\AllenBrainVisualEphys\';
strDataTarget = strcat(strDisk,strDataSource,'Aggregates',filesep,'AggSes',getDate,'.mat');
%load metadata
sCSV = loadcsv(strcat(strDisk,strDataSource,'sessions.csv'));
%select VIP
vecSessions = sCSV.id;
%vecSessions = 756029989;
strFormat = 'ecephys_session_%d';

%% loop through sessions
for intSes=1:numel(vecSessions)
	
	%% load data
	strDataFile = sprintf(strFormat,vecSessions(intSes));
	strTargetFile = strcat(strDisk,fullpath(strDataSource,'nwb_files',[strDataFile,'.nwb']));
	fclose('all');
	fprintf('Processing session %d/%d: %s [%s]\n',intSes,numel(vecSessions),strDataFile,getTime);
	nwb = nwbRead(strTargetFile);
	sNWB = ExpandNWB(nwb,[],false);
	%sNWB2 = PruneStruct(sNWB);
	
	%% stimuli, nat scenes
	if isfield(sNWB.intervals,'natural_scenes_presentations')
		sNS = sNWB.intervals.natural_scenes_presentations.vectordata;%orientation, SF, TF, etc
		%remove non-scenes
		indBlanks = sNS.frame==-1;
		vecRealScenes = sNS.frame(~indBlanks);
		intScenes = numel(unique(vecRealScenes));
		vecNatSceneReps(intSes) = numel(vecRealScenes)/intScenes;
		fprintf('Ses %d (%s) has %d repetitions of %d unique natural scenes\n',intSes,sprintf(strFormat,vecSessions(intSes)),vecNatSceneReps(intSes),intScenes)
	end
	%% stimuli, nat scenes
	if isfield(sNWB.intervals,'natural_scenes_presentations')
		sNM = sNWB.intervals.natural_movie_one_presentations.vectordata;%orientation, SF, TF, etc
		%remove non-scenes; 30fps, 30s long(?)
		indBlanks = sNM.frame==-1;
		vecMovieFrames = sNM.frame(~indBlanks);
		intReps = sum(diff(vecMovieFrames)<0)+1;
		vecNatSceneReps(intSes) = numel(vecRealScenes)/intScenes;
		fprintf('Ses %d (%s) has %d repetitions of natural movie one\n',intSes,sprintf(strFormat,vecSessions(intSes)),intReps);
	end
	
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
	
	%% subject
	sSubject = sNWB.general_subject;
	
	%% running
	vecRunningSpeed = sNWB.processing.running.nwbdatainterface.running_speed.data; %speed (cm/s)
	vecRunningTime = sNWB.processing.running.nwbdatainterface.running_speed.timestamps; %time (s)
	
	%% stimuli, DG
	structStimAgg = struct;
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
	if isfield(sNWB.intervals,'natural_scenes_presentations')
		strSource = 'natural_movie_one_presentations';
	else
		%no gratings
		strSource = '';
	end
	if ~isempty(strSource)
		sSource = sNWB.intervals.(strSource);
		sNS = sSource.vectordata;%orientation, SF, TF, etc
		sNS.Name = strSource;
		sNS.stimIdx = sNS.frame+1;
		sNS.startT = sSource.start_time;
		sNS.stopT = sSource.stop_time;
		%vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
		structStimAgg.sNS = sNS;
	end
	
	%% stimuli, nat movies
	if isfield(sNWB.intervals,'natural_movie_one_presentations')
		strSource = 'natural_movie_one_presentations';
	elseif isfield(sNWB.intervals,'natural_movie_two_presentations')
		strSource = 'natural_movie_one_presentations';
	else
		%no gratings
		strSource = '';
	end
	if ~isempty(strSource)
		sSource = sNWB.intervals.(strSource);
		sNM = sSource.vectordata;%orientation, SF, TF, etc
		sNM.Name = strSource;
		sNM.stimIdx = sNM.frame+1;
		sNM.startT = sSource.start_time;
		sNM.stopT = sSource.stop_time;
		%vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
		structStimAgg.sNM = sNM;
	end
	
	%% remove trials during invalid times
	sRemoveIntervals = sNWB.intervals_invalid_times;
	
	%% eye-tracking
	%not in all sessions!
	vecPupilSize = [];
	vecPupilT = [];
	matPupilXY = [];
	if isfield(sNWB.processing,'filtered_gaze_mapping')
		vecPupilT = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.pupil_area.timestamps; %sec
		vecPupilSize = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.pupil_area.data; %pixel^2
		matPupilXY = sNWB.processing.filtered_gaze_mapping.nwbdatainterface.screen_coordinates.data; %cm
	end
	
	%% get optogenetic modulation
	%cull cells
	vecUseCells = find(indInclude & contains(cellAreas(:),'VIS'));
	intNeurons = numel(vecUseCells);
	
	%gather data
	sSes(intSes).Exp = strDataFile;
	sSes(intSes).vecUseCells = vecUseCells;
	sSes(intSes).cellSpikes = cellSpikes(vecUseCells);
	sSes(intSes).cellAreas = cellAreas(vecUseCells);
	
	%stims
	sSes(intSes).structStimAgg = structStimAgg;
	
	%behavior
	sSes(intSes).vecPupilT = vecPupilT;
	sSes(intSes).vecPupilSize = vecPupilSize;
	sSes(intSes).matPupilXY = matPupilXY;
	sSes(intSes).vecRunningTime = vecRunningTime;
	sSes(intSes).vecRunningSpeed = vecRunningSpeed;
	
	%to-be-removed-interval structure
	sSes(intSes).sRemoveIntervals = sRemoveIntervals;
end


%% saving data
fprintf('Saving data to "%s"\n',strDataTarget);
save(strDataTarget,'sSes','-v7.3');
fprintf('Done.\n');