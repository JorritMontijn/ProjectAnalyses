
%% define locations
strDisk = 'F:';
strDataSource = '\Data\Processed\AllenBrainVisualEphys\';
strDataTarget = strcat(strDisk,strDataSource,'Aggregates',filesep,'AggSes',getDate,'.mat');
%load metadata
sCSV = loadcsv(strcat(strDisk,strDataSource,'sessions.csv'));
%select VIP
vecSessions = sCSV.id(contains(sCSV.genotype,'Vip'));
strFormat = 'session_%d';
vecNatSceneReps = zeros(size(vecSessions));

%% loop through sessions
for intSes=1:numel(vecSessions)
	%% skip if not laser opto
	if vecSessions(intSes) < 789848216
		continue;
	end
	
	%% load data
	strDataFile = sprintf(strFormat,vecSessions(intSes));
	strTargetFile = strcat(strDisk,strDataSource,strDataFile,filesep,strDataFile,'.nwb');
	fclose('all');
	nwb = nwbRead(strTargetFile);
	sNWB = ExpandNWB(nwb);
	%sNWB2 = PruneStruct(sNWB);
	
	%% stimuli, nat scenes
	if isfield(sNWB.intervals,'natural_scenes_presentations')
		sNS = sNWB.intervals.natural_scenes_presentations.vectordata;%orientation, SF, TF, etc
		%remove non-scenes
		indBlanks = sNS.frame==-1;
		vecRealScenes = sNS.frame(~indBlanks);
		intScenes = numel(unique(vecRealScenes));
		vecNatSceneReps(intSes) = numel(vecRealScenes)/intScenes;
		fprintf('Ses %d (%s) has %d repetitions of %d unique natural scenes\n',intSes,vecSessions(intSes),intScenes/numel(vecRealScenes),intScenes)
	end
	continue;
	
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
	
	%% opto stim
	%{
	Identifying Cre+ units
	Now that we know how to align spikes, we can start assessing which units are reliably driven by the optotagging stimulus and are likely to be Cre+.
	There are a variety of ways to do this, but these are the most important things to keep in mind:
		Spikes that occur precisely at the start or end of a light pulse are likely artifactual, and need to be ignored.
		The bright blue light required for optotagging can be seen by the mouse, so any spikes that occur more than 40 ms after the stimulus onset may result from retinal input, as opposed to direct optogenetic drive.
		The rate of false negatives (Cre+ cells that are not light-driven) will vary across areas, across depths, and across sessions. We've tried our best to evenly illuminate the entire visual cortex, and to use light powers that can drive spikes throughout all cortical layers, but some variation is inevitable.
	For these reasons, we've found that the 10 ms pulses are the most useful stimulus for finding true light-evoked activity. These pulses provide a long enough artifact-free window to observe light-evoked spikes, but do not last long enough to be contaminated by visually driven activity.
	Using the DataArray we created previously, we can search for units that increase their firing rate during the 10 ms pulse:
	baseline = da.sel(time_relative_to_stimulus_onset=slice(-0.01,-0.002))
	baseline_rate = baseline.sum(dim='time_relative_to_stimulus_onset').mean(dim='trial_id') / 0.008
	evoked = da.sel(time_relative_to_stimulus_onset=slice(0.001,0.009))
	evoked_rate = evoked.sum(dim='time_relative_to_stimulus_onset').mean(dim='trial_id') / 0.008
	%}
	
	sOptoStim = sNWB.processing.optotagging.dynamictable.optogenetic_stimulation; %speed (cm/s)
	sOptoParams = sOptoStim.vectordata;
	sOptoMeta =sNWB.processing.optotagging.nwbdatainterface.optotagging; %time (s)
	vecOptoDur = sOptoParams.duration;
	vecOptoLvl = sOptoParams.level;
	indUse10ms = vecOptoDur > 0.007 & vecOptoDur < 0.02;
	vecOptoEventsT = sOptoStim.start_time(indUse10ms);
	
	%% running
	vecRunningSpeed = sNWB.processing.running.nwbdatainterface.running_speed.data; %speed (cm/s)
	vecRunningTime = sNWB.processing.running.nwbdatainterface.running_speed.timestamps; %time (s)
	
	%% stimuli, DG
	vecAllStimTimes = sNWB.processing.stimulus.nwbdatainterface.timestamps.data;
	if isfield(sNWB.intervals,'drifting_gratings_presentations')
		sDG = sNWB.intervals.drifting_gratings_presentations.vectordata;%orientation, SF, TF, etc
		vecDG_stimIdx = table2array(sNWB.intervals.drifting_gratings_presentations.timeseries(:,1))+1;
		vecDG_startT = sNWB.intervals.drifting_gratings_presentations.start_time;
		vecDG_stopT = sNWB.intervals.drifting_gratings_presentations.stop_time;
		vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
	elseif isfield(sNWB.intervals,'drifting_gratings_75_repeats_presentations')
		sDG = sNWB.intervals.drifting_gratings_75_repeats_presentations.vectordata;%orientation, SF, TF, etc
		vecDG_stimIdx = table2array(sNWB.intervals.drifting_gratings_75_repeats_presentations.timeseries(:,1))+1;
		vecDG_startT = sNWB.intervals.drifting_gratings_75_repeats_presentations.start_time;
		vecDG_stopT = sNWB.intervals.drifting_gratings_75_repeats_presentations.stop_time;
		vecDG_stimT = vecAllStimTimes(vecDG_stimIdx); %same as vecDG_startT
	else
		%no gratings
	end
	
	%% remove trials during invalid times
	sRemoveIntervals = sNWB.intervals_invalid_times;
	if ~isempty(sRemoveIntervals)
		error
	end
	
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
	dblUseMaxDur = 1;%100/1000;
	vecPreOnsets = vecOptoEventsT - dblUseMaxDur/2;
	vecZetaP = nan(1,intNeurons);
	hTic = tic;
	for intNeuron=1:intNeurons
		if toc(hTic) > 5,hTic=tic;fprintf('Proc %d.1: Neuron %d/%d [%s]\n',intSes,intNeuron,intNeurons,getTime);end
		intUnit = vecUseCells(intNeuron);
		vecZetaP(intNeuron) = getZeta(cellSpikes{intUnit},vecPreOnsets,dblUseMaxDur);
	end
	
	intResampNum = 100;
	intPlot = 0;
	intLatencyPeaks = 4;
	vecZetaP2 = nan(1,intNeurons);
	matLatencies = nan(4,intNeurons);
	for intNeuron = find(vecZetaP<0.05)
		if toc(hTic) > 5,hTic=tic;fprintf('Proc %d.2: Neuron %d/%d [%s]\n',intSes,intNeuron,intNeurons,getTime);end
		intUnit = vecUseCells(intNeuron);
		[vecZetaP2(intNeuron),vecLatencies] = getZeta(cellSpikes{intUnit},vecPreOnsets,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks);
		matLatencies(:,intNeuron) = vecLatencies;
	end
	
	%gather data
	sSes(intSes).vecUseCells = vecUseCells;
	sSes(intSes).vecZetaP = vecZetaP;
	sSes(intSes).vecZetaP2 = vecZetaP2;
	sSes(intSes).matLatencies = matLatencies;
	sSes(intSes).cellSpikes = cellSpikes(vecUseCells);
	sSes(intSes).cellAreas = cellAreas(vecUseCells);
	sSes(intSes).vecOptoEventsT = vecOptoEventsT;
	sSes(intSes).dblUseMaxDur = dblUseMaxDur;
end
vecNatSceneReps
return
%% saving data
fprintf('Saving data to "%s"\n',strDataTarget);
save(strDataTarget,'sSes');
fprintf('Done.\n');