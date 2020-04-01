%cellMultiSes cellAggregate vecBlock

%% define variables
clear all;
close all;

boolUseNeuropilSubtraction = true;
strTargetPath = 'D:\Data\Processed\StimDetectionAgg';
strMasterPath = 'D:\Data\Processed\imagingdata';
strEyeTrackingPath = 'D:\Data\Processed\imagingvideo';
for intMouse = 1:8
	clearvars -except strTargetPath strMasterPath strEyeTrackingPath intMouse boolUseNeuropilSubtraction
	if intMouse == 1
		strSes = '20140207';
		vecRecordings = 1:8;
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intMouse == 2
		strSes = '20140314';
		vecRecordings = 1:7;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 3
		strSes = '20140425';
		vecRecordings = 1:8;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == -1
		%return; %exclude?: bad behavior, weird signals
		strSes = '20140430';
		vecRecordings = 2:8;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 4
		strSes = '20140507';
		vecRecordings = 1:3;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 5
		strSes = '20140530';
		vecRecordings = 1:7;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 6
		strSes = '20140604';
		vecRecordings = 1:4;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 7
		strSes = '20140711';
		vecRecordings = 1:8;
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intMouse == 8
		strSes = '20140715';
		vecRecordings = 1:4;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	end
	
	%put in structure
	sIn.strSes = strSes;
	sIn.strMasterPath = strMasterPath;
	sIn.vecRecordings = vecRecordings;
	
	%define parameters
	cellFields = {'neuron','astrocyte','PV'};
	intObjectTypes = numel(cellFields);
	intRemTrialType = 0; %0=all trials; 1=only short RT; 2=only long RT
	dblCutOffRT = 1.2; %seconds
	intAgg = 0; %will calculate everything over aggregate session if ~=0
	sParams = struct;
	sParams.boolSaveData = true;
	
	%calculate block definitions
	intStartSes = 1;
	intStopSes = length(vecRecordings);
	vecBlockTypes = unique(vecBlock);
	vecBiDirSes = [nan 9];
	intNumBlocks = length(vecBlockTypes);
	matObjectNumPreRem = zeros(intNumBlocks,intObjectTypes);
	matObjectNumPostRem = zeros(intNumBlocks,intObjectTypes);
	cellKeepList = cell(intNumBlocks,intObjectTypes);
	vecFirstBlock = length(vecBlockTypes);
	vecLastBlock = length(vecBlockTypes);
	for intBlockType = vecBlockTypes
		vecFirstBlock(intBlockType) = find(vecBlock == intBlockType,1,'first');
		vecLastBlock(intBlockType) = find(vecBlock == intBlockType,1,'last');
	end
	
	%% get stim detection data per session
	%create links to files
	intCounter = 0;
	cellAggregate = cell(1,length(vecRecordings));
	cellEyeTracking = cell(1,length(vecRecordings));
	for intRec=vecRecordings
		intCounter = intCounter + 1;
		cellAggregate{intCounter} = sprintf('%s%s%s%sxyt%02d%s%sxyt%02d_ses',strMasterPath,filesep,strSes,filesep,intRec,filesep,strSes,intRec);
		cellEyeTracking{intCounter} = sprintf('%s%s%s%sEyeTrackData_%sxyt%02d',strEyeTrackingPath,filesep,strSes,filesep,strSes,intRec);
	end
	
	%check if NPS needs to be added
	if boolUseNeuropilSubtraction
		strSes = [strSes 'NPS']; %#ok<AGROW>
	end
	
	%aggregate all behavioral sessions
	fprintf('Loading all data and creating inclusion list for [%s]..\n',cell2str(cellFields));
	for intThisBlock=1:intNumBlocks
		cellPresenceList{intThisBlock,numel(cellFields)} = [];
		cellResponseList{intThisBlock,numel(cellFields)} = [];
	end
	cellSes = cell(1,numel(cellAggregate));
	for intFile=1:numel(cellAggregate)
		%% load data
		sLoad = load([cellAggregate{intFile} '.mat']);
		ses = sLoad.ses;
		if exist([cellEyeTracking{intFile} '.mat'],'file')
			sLoad = load([cellEyeTracking{intFile} '.mat']);
			sEyeTracking = sLoad.sEyeTracking;
		else
			sEyeTracking = [];
		end
		
		%% remove trials
		if intRemTrialType == 1 %only short
			indRemTrials = ses.structStim.vecTrialRespSecs-ses.structStim.SecsOn > dblCutOffRT;
		elseif intRemTrialType == 2 %only long
			indRemTrials = ses.structStim.vecTrialRespSecs-ses.structStim.SecsOn < dblCutOffRT;
		else indRemTrials = [];
		end
		if ~isempty(indRemTrials)
			cellNames=fieldnames(ses.structStim);
			for intField=1:length(cellNames)
				strField = cellNames{intField};
				ses.structStim.(strField) = ses.structStim.(strField)(~indRemTrials);
			end
		end
		
		%% create inclusion list
		for intField=1:numel(cellFields)
			strField = cellFields{intField};
			intThisBlock = vecBlock(intFile);
			if ~isfield(ses,strField),continue;end
			matObjectNumPreRem(intThisBlock,intField) = numel(ses.(strField));
			indPresenceList = cellPresenceList{intThisBlock,intField};
			indResponseList = cellResponseList{intThisBlock,intField};
			if isempty(indPresenceList),indPresenceList = true(1,numel(ses.(strField)));end
			if isempty(indResponseList),indResponseList = false(1,numel(ses.(strField)));end
			for intObject=1:numel(ses.(strField))
				if strcmp(ses.(strField)(intObject).strPresence,'absent')
					indPresenceList(intObject) = false;
				end
				if ~strcmp(ses.(strField)(intObject).strRespType,'silent')
					indResponseList(intObject) = true;
				end
			end
			cellPresenceList{intThisBlock,intField} = indPresenceList;
			cellResponseList{intThisBlock,intField} = indResponseList;
		end
		
		%% load and add eyetracking data
		%{
			offsets present:
			20140207xyt01,xyt03
			20140530xyt04
			20140604xyt01-05
			20140711xyt02-08
			20140715xyt01-04

			offsets missing:
			20140207xyt02,xyt04-08
			20140314xyt02-07
			20140425xyt01-08
			20140507xyt01-03
			20140530xyt01-03,xyt05-07

			pupil detection
			20140207
			xyt01-3: good
			xyt04-8: bad

			20140314
			xyt02-07: good

			20140425
			xyt01-08: good

			20140507
			xyt01-03: good
		%}
		if isempty(sEyeTracking)
			%create dummy data
			sEyeTracking.strSes = strSes;
			sEyeTracking.vecWeight = nan(1,2);
			sEyeTracking.vecRoundness = nan(1,2);
			sEyeTracking.vecPupilLuminance = nan(1,2);
			sEyeTracking.vecPosX = nan(1,2);
			sEyeTracking.vecPosY = nan(1,2);
			sEyeTracking.vecArea = nan(1,2);
			sEyeTracking.vecEvents = nan(1,2);
			sEyeTracking.indEvents = nan(1,2);
		end
		if strcmp(ses.session,'20140207') && ismember(ses.recording,[1 3])
			boolEndPresent = true;
		elseif strcmp(ses.session,'20140530') && ismember(ses.recording,4)
			boolEndPresent = true;
		elseif strcmp(ses.session,'20140604') && ismember(ses.recording,1:5)
			boolEndPresent = true;
		elseif strcmp(ses.session,'20140711') && ismember(ses.recording,2:8)
			boolEndPresent = true;
		elseif strcmp(ses.session,'20140715') && ismember(ses.recording,1:4)
			boolEndPresent = true;
		else
			boolEndPresent = false;
		end
		
		%{
		if strcmp(ses.session,'20140207') && ismember(ses.recording,[2 7])
			boolEndPresent = false;
		elseif strcmp(ses.session,'20140711') && ismember(ses.recording,1)
			boolEndPresent = false;
		else
			boolEndPresent = true;
		end
		%}
		
		%put in output
		ses.sEyeTracking = getResampledEyeTrackingData(sEyeTracking,ses,boolEndPresent);
		
		%% assign to output
		cellSes{intFile} = ses;
		
		%msg
		fprintf('Processed [%s]\n',cellAggregate{intFile});
	end
	%%
	for intField=1:intObjectTypes
		strType = cellFields{intField};
		if ~isfield(ses,strType),continue;end
		for intThisBlock=1:intNumBlocks
			if strcmpi(strType,'astrocyte')
				cellKeepList{intThisBlock,intField} = cellPresenceList{intThisBlock,intField};
			else
				cellKeepList{intThisBlock,intField} = cellPresenceList{intThisBlock,intField};% & cellResponseList{intThisBlock,intField};
			end
		end
	
		
		%add list to parameter structure
		for intThisBlock=vecBlockTypes
			intBlockStart = find(vecBlock == intThisBlock,1,'first');
			intBlockStop = find(vecBlock == intThisBlock,1,'last');
			indKeepList = cellKeepList{intThisBlock,intField};
			%indKeepList = true(size(cellKeepList{intThisBlock}));
			intObjectNum = matObjectNumPreRem(intThisBlock,intField);
			intObjectsIncluded = sum(indKeepList);
			matObjectNumPostRem(intThisBlock,intField) = intObjectsIncluded;
			fprintf('>> Will include [%d/%d] %ss for block %d [%d - %d]\n',intObjectsIncluded,intObjectNum,strType,intThisBlock,intBlockStart,intBlockStop)
		end
	end
	
	%% build multi ses aggregate structure combination
	fprintf('\nBuilding multi-recording aggregate per block...\n')
	cellMultiSes = cell(1,intNumBlocks);
	for intBlockType = vecBlockTypes
		sSesAggregate = [];
		for intSes=vecFirstBlock(intBlockType):vecLastBlock(intBlockType)
			%load session file
			ses = cellSes{intSes};
			fprintf('Processing session %d/%d [recording block %d].. [%s]\n',intSes,length(vecBlock),intBlockType,getTime);
			
			%remove all excluded neurons
			for intField=1:intObjectTypes
				strField = cellFields{intField};
				if ~isfield(ses,strField),continue;end
				indKeepList = cellKeepList{intBlockType,intField};
				ses = doRecalcdFoF(ses,6,indKeepList,strField);
				
				%recalculate with NPS if requested
				if boolUseNeuropilSubtraction
					%recalculate dF/F0 with neuropil subtraction
					dblFraction=0.5;
					dblSecWindowSize=30;
					dblNeuropilSubtractionFactor=0.7;
					ses = doRecalcdFoF(ses,5,[],strField,dblFraction,dblSecWindowSize,dblNeuropilSubtractionFactor); %originally pre (5)
					
				else
					%recalculate dF/F0 without neuropil subtraction
					ses = doRecalcdFoF(ses,3,[],strField);
					
				end
				
				%recalculate spikes
				if strcmpi(strField,'astrocyte')
					dblTau = 10;
				else
					dblTau = 0.5;
				end
				dblSamplingFreq = ses.samplingFreq;
				dblTotDurSecs = length(ses.(strField)(1).dFoF)/dblSamplingFreq;
				for intObject=1:numel(ses.(strField))
					ptrTime = tic;
					[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(ses.(strField)(intObject).dFoF,dblSamplingFreq,dblTau);
					ses.(strField)(intObject).apFrames = apFrames;
					ses.(strField)(intObject).apSpikes = apSpikes;
					ses.(strField)(intObject).vecSpikes = vecSpikes;
					ses.(strField)(intObject).expFit = expFit;
					fprintf('%s %d/%d [duration %.1f seconds]; %d transients dectected; mean spiking rate is %.2f Hz\n',strField,intObject,numel(ses.(strField)),toc(ptrTime),sum(apSpikes),sum(apSpikes)/dblTotDurSecs)
					pause(0);
				end
			end
			%create aggregate
			sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate);
		end
		for intField=1:intObjectTypes
			strField = cellFields{intField};
			if ~isfield(ses,strField),continue;end
			sSesAggregate = doRecalcdFoF(sSesAggregate,6,[],strField); %remove neurons with extreme dF/F0 values
		end
		cellMultiSes{intBlockType} = sSesAggregate;
	end
	
	%% save data
	strDataFile = sprintf('dataPreProAggregate%s',strSes);
	strFullTarget = [strTargetPath filesep strDataFile];
	fprintf('Saving to: %s\n',strFullTarget);
	save(strFullTarget,'cellMultiSes','cellAggregate','vecBlock','-v7.3');
	fprintf('\b    Done! Time is %s\n\n',getTime)
end