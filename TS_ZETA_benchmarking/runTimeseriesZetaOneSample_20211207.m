%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
intResampNum = 251;%[100 200 500 1000 2000];
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = 'SimOneSample1TsZeta';
intNeurons = 100;
dblFracDiffSpikes = 1/2;
dblTau = 2;
dblTau0 = (63/1000);
dblNoise = 0.025;
dblSamplingFreq = 25;
dblSamplingInterval = 1/dblSamplingFreq;
boolQuick = false;
boolDirectQuantile = false;

%set indicator properties
sIndicatorProps = struct;
sIndicatorProps.dblTimescale = dblTau;
sIndicatorProps.dblNoise = dblNoise;

%% pre-allocate output variables
matTtest = nan(intNeurons,2);
matTsZeta = nan(intNeurons,2);
matAnova = nan(intNeurons,2);

%% generate data
hTicN = tic;
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
	end
	
	%% tuning params
	dblBaseRate = exprnd(1);
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	dblJitter = 5; %in ms'
	
	%% stimulus data
	dblStimDur = 2;
	dblPreBaseDur = 1;
	dblPostBaseDur = 1;
	dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
	dblUseTrialDur = dblTrialDur;
	intOris = 24;
	dblStepDeg = 360/intOris;
	vecOris = linspace(0,360-dblStepDeg,intOris);
	intReps = 10;
	intTrials = intOris*intReps;
	vecTrialAngles = nan(1,intTrials);
	for intRep=1:intReps
		vecTrialAngles(((intRep-1)*intOris+1):(intRep*intOris)) = vecOris(randperm(intOris));
	end
	vecTrialStart = 10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*intTrials+eps));
	matTrialT1 = cat(2,vecTrialStart',vecTrialStart'+dblStimDur);
	matTrialT1 = matTrialT1 + 0.05*rand(size(matTrialT1));
	
	%% generate data
	% generate bursts
	intAddSpikes1 = intTrials/5;
	
	% generate peak
	dblStartDelay = 0.1;
	[vecSpikeTimes1,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT1,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes1,dblStartDelay);
	[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes1,dblSamplingFreq,sIndicatorProps,boolQuick);
	%add empty end
	dblEndDur = vecTrialStart(end) + dblTrialDur*5;
	vecAddT = vecTimestamps(end):(1/dblSamplingFreq):dblEndDur;
	vecTimestamps = cat(2,vecTimestamps,vecAddT(2:end));
	vecdFoF = cat(2,vecdFoF,zeros(size(vecAddT(2:end))));
	
	%real+rand
	for intRunType=vecRunTypes
		%% get visual responsiveness
		%set derivative params
		if intRunType ==2
			vecJitterPerTrial = 2*linspace(dblUseTrialDur/intTrials,dblUseTrialDur,intTrials)';
			matEventTimes = bsxfun(@plus,matTrialT1,vecJitterPerTrial(randperm(numel(vecJitterPerTrial))));
		else
			matEventTimes = matTrialT1;
		end
		
		%ANOVA
		hTicA = tic;
		[vecRefT2,matTracePerTrial] = getTraceInTrial(vecTimestamps,vecdFoF,matEventTimes(:,1),dblSamplingInterval,dblUseTrialDur);
		dblBinWidth = median(diff(vecRefT2));
		dblAnovaP=anova1(matTracePerTrial,[],'off');
		dblAnovaDur = toc(hTicA);
		matAnova(intNeuron,intRunType) = dblAnovaP;
		
		%TS-ZETA
		hTicZ = tic;
		intPlot = 0;
		%continue;
		[dblZetaP,sZETA] = zetatstest(vecTimestamps,vecdFoF,matEventTimes,dblUseTrialDur,intResampNum,intPlot,boolDirectQuantile);
		%pause
		% assign data
		dblMeanP = sZETA.dblMeanP;
		dblMeanZ = -norminv(dblMeanP/2);
		dblZetaZ = sZETA.dblZETA;
		matTsZeta(intNeuron,intRunType) = dblZetaP;
		matTtest(intNeuron,intRunType) = dblMeanP;
	end
end

%% save
if boolSave
	save([strDataTargetPath strRec 'Q' num2str(boolDirectQuantile) '.mat' ],...
		'matAnova','matTtest','matTsZeta','strRec');
end
