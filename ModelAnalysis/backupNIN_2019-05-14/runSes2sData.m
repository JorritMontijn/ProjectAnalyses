%{
The variable 'output' contains both the V1 (output.nev) and V2 (output.plx) data.
The spike time (*.TimeStamps) are matrices with many rows and 4 columns:

col 1: the electrode number (1-96), electrode 2000 is the digital stimulus
	marker which indicates the onset of the stimulus;

col 2: sort code, where 0 and/or 255 are junk, but any other number is a
	valid SUA or MUA

col 3: the time at which the event occurred relative to stimulus onset

col 4: the trial number

The stimulus are presented for 1280 ms and separated by 1500 ms of blank screen (ISI).

The stimulus order is in the variable 'stim', where 1-8 are the 8 gratings
	(0:22.5:157.5); stimulus 0 is the blank (gray screen) between stimuli.
%}
%% set recording
for intLoadExp=[-21:-1:-27]
clearvars -except intLoadExp;
close all;
if intLoadExp == -21
	%S008
	cellMulti = {'20140410xyt01_ses.mat', '20140414xyt01_ses.mat', '20140415xyt01_ses.mat', '20140416xyt01_ses.mat', '20140424xyt01_ses.mat', '20140429xyt01_ses.mat', '20140501xyt01_ses.mat'};
	strExperiment = 'S008';
elseif intLoadExp == -22
	%S006
	cellMulti = {'20140514xyt01_ses.mat', '20140519xyt01_ses.mat', '20140520xyt03_ses.mat', '20140521xyt01_ses.mat', '20140522xyt01_ses.mat', '20140527xyt01_ses.mat', '20140528xyt01_ses.mat'};
	strExperiment = 'S006';
elseif intLoadExp == -23
	%S010
	cellMulti = {'20140527xyt04_ses.mat', '20140602xyt01_ses.mat', '20140605xyt01_ses.mat', '20140610xyt01_ses.mat', '20140611xyt02_ses.mat'};
	strExperiment = 'S010';
elseif intLoadExp == -24
	%P011
	cellMulti = {'20140804xyt02_ses.mat', '20140808xyt01_ses.mat', '20140811xyt01_ses.mat', '20140811xyt02_ses.mat', '20140815xyt01_ses.mat', '20140819xyt01_ses.mat', '20140821xyt01_ses.mat', '20140826xyt01_ses.mat', '20140902xyt01_ses.mat', '20140902xyt02_ses.mat'};
	strExperiment = 'P011';
elseif intLoadExp == -25
	%P012
	cellMulti = {'20140714xyt01_ses.mat', '20140716xyt01_ses.mat', '20140722xyt01_ses.mat', '20140728xyt01_ses.mat', '20140728xyt02_ses.mat', '20140801xyt01_ses.mat'};
	strExperiment = 'P012';
elseif intLoadExp == -26
	%P016
	cellMulti = {'20141121xyt01_ses.mat', '20141126xyt01_ses.mat', '20141202xyt01_ses.mat', '20141205xyt04_ses.mat', '20141208xyt02_ses.mat', '20141212xyt02_ses.mat', '20141215xyt01_ses.mat', '20141219xyt01_ses.mat', '20141219xyt02_ses.mat'};
	strExperiment = 'P016';
elseif intLoadExp == -27
	%P014
	cellMulti = {'20141202xyt04_ses.mat', '20141205xyt01_ses.mat', '20141208xyt03_ses.mat', '20141212xyt03_ses.mat', '20141215xyt03_ses.mat', '20141219xyt04_ses.mat', '20141219xyt05_ses.mat', '20141224xyt02_ses.mat'};
	strExperiment = 'P014';
end

%% set paths
strSourceDir = 'D:\Data\Processed\JoGu\';
strTargetDir = 'A:\SimAggregates\';
boolPlot = false;
boolTestTuning = false;

%% load aggregate data
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
%get data
[cellSes, cellSesInclude, indInclude] = getIncludeNeurons(strSourceDir, cellMulti, boolPlot, boolTestTuning);
drawnow;

%build aggregate
sSesAggregate = [];
for intSes=1:length(cellSes)
	%load session file
	ses = cellSes{intSes};
	
	%check for SecsOn/SecsOff
	if ~isfield(ses.structStim,'SecsOn')
		ses.structStim.SecsOff = ses.structStim.FrameOff / ses.samplingFreq;
		ses.structStim.SecsOn = ses.structStim.FrameOn / ses.samplingFreq;
	end
	
	%create aggregate
	sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate);
	clear ses;
end
dblFraction = 0.5;
dblSecWindowSize = 30;
dblNeuropilSubtractionFactor = 0.5;
[sSesAggregate,indKeepList] = doRecalcdFoF(sSesAggregate,8,[],'neuron',dblFraction,dblSecWindowSize,dblNeuropilSubtractionFactor);


%% meta data
%identifier
dblTransformToDeltaT = 1/sSesAggregate.samplingFreq;
structStim = sSesAggregate.structStim;

%% transform data
intCellsV1 = numel(sSesAggregate.neuron);
intTrials = numel(structStim.Orientation);
vecOrientations = unique(structStim.Orientation);
intStimTypes = numel(vecOrientations);
vecTrialOris = structStim.Orientation;
vecTrialStimType = nan(1,intTrials);
vecTrialStimRep = nan(1,intTrials);
matModelResp = nan(intCellsV1,intTrials);
for intStimType=1:intStimTypes
	dblThisOri = vecOrientations(intStimType);
	indTheseTrials = vecTrialOris==dblThisOri;
	vecTrialStimType(indTheseTrials) = intStimType;
	vecTrialStimRep(indTheseTrials) = 1:sum(indTheseTrials);
	vecTheseTrials = find(indTheseTrials);
	for intTrialIdx=1:numel(vecTheseTrials)
		intTrial = vecTheseTrials(intTrialIdx);
		intStart = structStim.FrameOn(intTrial);
		intStop = structStim.FrameOff(intTrial);
		for intNeuron=1:intCellsV1
			matModelResp(intNeuron,intTrial) = mean(sSesAggregate.neuron(intNeuron).dFoF(intStart:intStop));
		end
	end
end

%% assign to variables
dblStimDur=structStim.SecsOff(1) - structStim.SecsOn(1);
dblTrialDur=structStim.SecsOn(2) - structStim.SecsOn(1);

vecTrialStartSecs = structStim.SecsOn(:)';
vecStimStartSecs = structStim.SecsOn(:)';
vecStimStopSecs = structStim.SecsOff(:)';
vecTrialEndSecs = structStim.SecsOff(:)';
vecTrialStimRep = vecTrialStimRep(:)';
vecTrialStimType = vecTrialStimType(:)';

%% V2 placeholders
intCellsV2 = 0;

%build V1/V2 combos
cellSpikeTimesV2 = {};
intCortexCells = intCellsV1 + intCellsV2;
vecCellArea = ones(1,intCellsV1);

%% build sData
sData = struct;

%identifier
sData.strConnFile = ['JoGuExp' strExperiment];
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:(vecTrialStartSecs(end)+dblTrialDur);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = 'OriDrift8';
sData.matModelResp = matModelResp;

sData.vecTrialStartSecs = vecTrialStartSecs(:)';
sData.vecStimStartSecs = vecStimStartSecs(:)';
sData.vecStimStopSecs = vecStimStopSecs(:)';
sData.vecTrialEndSecs = vecTrialEndSecs(:)';
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intCellsV1;
sData.intCellsV2 = 0;
sData.intCortexCells = intCortexCells;
sData.vecCellArea = vecCellArea;

%sStimInputs
sData.dblTrialDur=dblTrialDur;
sData.vecTrialOris=vecOrientations(vecTrialStimType);
sData.vecTrialOriIdx=vecTrialStimType;
sData.vecTrialStimRep = vecTrialStimRep;
sData.vecStimTypeOris=vecOrientations;
sData.vecStimTypeSFs=0.05*ones(size(vecOrientations));
sData.vecStimTypeTFs=1.5*ones(size(vecOrientations));
sData.vecStimTypeContrasts = 100*ones(size(vecOrientations));
sData.vecStimTypeLuminance = 100*ones(size(vecOrientations));

%% build stimulus properties
cellParamTypes{1} = sData.vecStimTypeOris;
cellParamTypes{2} = sData.vecStimTypeSFs;
cellParamTypes{3} = sData.vecStimTypeTFs;
cellParamTypes{4} = sData.vecStimTypeContrasts ;
cellParamTypes{5} = sData.vecStimTypeLuminance;
matStimTypeCombos = buildStimCombos(cellParamTypes);

sData.matStimTypeCombos=matStimTypeCombos;

%% save
strOutputFile = strcat('Simulation_xAreaExperiment_JoGu',strExperiment,'_',getDate,'.mat');
fprintf('%s transformation complete, saving data to <%s%s> [%s]\n',strExperiment,strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
fprintf('Done! [%s]\n\n',getTime);
end


