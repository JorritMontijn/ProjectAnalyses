%run data with NOT
%S1L3
%S2L1
%S2L4 => eye tracking wrong?
%S3L1 => probably not NOT
%S3L3
cellRun = {'S1L3','S2L1','S3L3'};
cellRun = {};
cellRun(end+1:end+2) = {'S2L4','S3L1'};

%get data
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
strDataTarget = 'F:\Data\Processed\PupilNOT\';
clear sFiles;
for intFile=1:numel(cellRun)
	sFile = dir([strDataSourcePath '*' cellRun{intFile} '*.mat']);
	sFiles(intFile) = sFile(1);
end
cellFiles = {sFiles(:).name}';
close all;

%% pre-allocate
cellRunAreas = {...
	'lateral geniculate',...Area 1
	'Primary visual',...Area
	'Lateral posterior nucleus',...Area
	'Anterior pretectal nucleus',...Area
	'Nucleus of the optic tract',...Area
	'Superior colliculus',...Area
	'Anteromedial visual',...Area
	'posteromedial visual',...Area
	'Anterolateral visual',...Area
	'Lateral visual',...Area
	'Rostrolateral area',...Area
	'Anterior area',...Area
	'Subiculum',...Area
	'Field CA1',...Area
	'Field CA2',...Area
	'Field CA3',...Area
	'Dentate gyrus',...Area
	'Retrosplenial'...Area
	};

dblFrameDur = 0.0099;
vecPlotWindow = [-1 1.5];
vecBinsT = vecPlotWindow(1):dblFrameDur:vecPlotWindow(2);

%% go through files
for intFile=1:numel(cellFiles)
	strTarget = cellFiles{intFile};
	strRec = strTarget((end-10):(end-7));
	sLoad = load(fullfile(strDataSourcePath,strTarget));
	sAP = sLoad.sAP;
	fprintf('Processing pupil STA for "%s" [%s]\n',strTarget,getTime);
	%% run pupil prep header
	boolSkipSpikes = true;
	runPupilHeader;
	%outputs:
	vecFixedPupilTime;
	vecOrientation;
	vecPupilStimOnTime;
	vecGlitchesPerTrials;
	vecStimOnTime;
	vecMoveX;
	vecMoveY;
	vecMoveR;
	
	%% spike-triggered pupil movement
	%z-score
	varTypeCV = 1;
	dblLambda = 200;
	% regroup
	vecRunAreaIdx = groupInto(cellClustAreas,cellRunAreas);
	[varDataOut1,vecUnique,vecCounts,cellSelect] = label2idx(vecRunAreaIdx);
	
	%assign neurons that belong to an area with fewer than 5 neurons to "other"
	indRemGroups = vecCounts < 2;
	indRemove = any(cell2mat(cellSelect(indRemGroups)),2);
	vecRunAreaIdx(indRemove) = 0;
	vecUnique(indRemGroups) = [];
	
	%analysis
	vecOrientation;
	vecPupilStimOnTime;
	vecGlitchesPerTrials;
	vecStimOnTime;
	
	vecMoveX;
	vecMoveY;
	vecMoveR;
	
	cellArea = {sAP.sCluster.Area};
	cellSpikes = {sAP.sCluster.SpikeTimes};
	vecWindow = round(-0.5/dblFrameDur):round(1/dblFrameDur);
	vecWindowT = vecWindow*dblFrameDur;
	vecN = nan(1,intNeurons);
	matM = nan(numel(vecWindow),intNeurons,3);
	matS = nan(numel(vecWindow),intNeurons,3);
	hTic=tic;
	dblMsgInterval = 300;
	for intNeuron=1:intNeurons
		if toc(hTic) > dblMsgInterval
			hTic=tic;
			fprintf('  Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
		vecSpikes = cellSpikes{intNeuron};
		intSpikes = numel(vecSpikes);
		matPupil = nan(intSpikes,numel(vecWindow),3);
		intPupN = numel(vecFixedPupilTime);
		parfor intSpike=1:numel(vecSpikes)
			intFrameIdx = 1+sum(vecFixedPupilTime<vecSpikes(intSpike));
			vecSelect = vecWindow + intFrameIdx;
			if vecSelect(1) < 1 || vecSelect(end) > intPupN,continue;end
			for intVar=1:3
				if intVar==1
					matY = vecMoveR;
				elseif intVar==2
					matY = vecMoveX;
				else
					matY = vecMoveY;
				end
				matPupil(intSpike,:,intVar) = matY(vecSelect);
			end
		end
		indRem = any(isnan(matPupil(:,:,1)),2);
		matPupil(indRem,:,:) = [];
		
		%plot
		vecN(intNeuron) = sum(~indRem);
		for intVar=1:3
			matM(:,intNeuron,intVar) = mean(matPupil(:,:,intVar),1);
			matS(:,intNeuron,intVar) = std(matPupil(:,:,intVar),[],1);
		end
		%errorbar(vecWindow*dblFrameDur,vecM,vecS/sqrt(sum(~indRem)));
	end
	toc(hTic)
	%{
	vecN(vecN<1000)=0;
	for intVar=1:3
		subplot(2,3,intVar)
		plot(vecWindow*dblFrameDur,matM(:,:,intVar).*sqrt(vecN))
	end
	%}
	%% aggregate data
	sPupilSTA = struct;
	sPupilSTA.File = strTarget;
	sPupilSTA.Rec = strRec;
	sPupilSTA.Info = '[(time) x (neuron) x (radius-change,x-change,y-change)]';
	
	sPupilSTA.matM = matM;
	sPupilSTA.matS = matS;
	sPupilSTA.vecN = vecN;
	sPupilSTA.vecWindowT = vecWindowT(:);
	sPupilSTA.cellArea = cellArea;
	
	
	%% save data
	strFile = sprintf('PupilSTA_%s_%s.mat',strRec,getDate);
	strTargetFile = fullfile(strDataTarget,strFile);
	save(strTargetFile,'sPupilSTA');
	fprintf('Pupil STA data saved to "%s" [%s]\n',strTargetFile,getTime);
end

