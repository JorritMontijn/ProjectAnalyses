%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta\Data\';
vecRunTypes = [1 2];
intResampNum = 250;
boolSave = true;
strFigPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta\Data\';

%% set variables
intNeurons = 100;
vecSampFreqs = 2.^(3:6);
intSampFreqNum = numel(vecSampFreqs);
vecTrialDur = [2.^(-5:0)];
intTrialDurNum = numel(vecTrialDur);
vecTau = [2.^(-1:6)];
dblTau0 = (63/1000);
dblNoise = 0.025;
intTauNum = numel(vecTau);
boolQuick = false;

%% pre-allocate output variables
matTtest = nan(intNeurons,2,intTrialDurNum);
matZeta = nan(intNeurons,2,intTrialDurNum);
matTsZeta = nan(intNeurons,2,intTrialDurNum,intTauNum,intSampFreqNum);
%% run
for intTrialDurIdx=1:intTrialDurNum
	dblTrialDur=vecTrialDur(intTrialDurIdx);
	
%% generate data
for intNeuron=1:intNeurons
	%% message
	fprintf('Processing trial dur %d/%d; neuron %d/%d [%s]\n',intTrialDurIdx,intTrialDurNum,intNeuron,intNeurons,getTime);
	hTic=tic;
	
	%% stimulus data
	dblStimDur = 0.5*dblTrialDur;
	dblPreBaseDur = 0.25*dblTrialDur;
	dblPostBaseDur = 0.25*dblTrialDur;
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
	matTrialT = cat(2,vecTrialStart',vecTrialStart'+dblStimDur);
	
	%% generate data
	dblBaseRate = exprnd(1);
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	
	% generate peak
	dblJitter = 5; %in ms'
	intAddSpikes = round(intTrials*0.1);
	[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes,10/1000);
	
	%% go through sampling freqs
	for intTauIdx=1:intTauNum
		dblTau = vecTau(intTauIdx);
	for intSampFreqIdx=1:intSampFreqNum
		dblSamplingFreq = vecSampFreqs(intSampFreqIdx);
		if dblTrialDur < (1/dblSamplingFreq)
			continue;
		end
		
		%fprintf('Running tau %d/%d, samp freq %d/%d [%s]\n',intTauIdx,intTauNum,intSampFreqIdx,intSampFreqNum,getTime);
		
		%set indicator properties
		sIndicatorProps = struct;
		sIndicatorProps.dblTimescale = dblTau;
		sIndicatorProps.dblNoise = dblNoise;
		[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick);
		%add empty end
		dblEndDur = vecTrialStart(end) + dblTrialDur*5;
		vecAddT = vecTimestamps(end):(1/dblSamplingFreq):dblEndDur;
		vecTimestamps = cat(2,vecTimestamps,vecAddT(2:end));
		vecdFoF = cat(2,vecdFoF,zeros(size(vecAddT(2:end))));
		
		%real+rand
		for intRunType=vecRunTypes
			
			%randomize
			if intRunType ==2
				vecJitter = 2*dblTrialDur*((rand(size(vecTrialStart))-1/2)*2);
				vecUseTrialStart = vecTrialStart + vecJitter;
				strRand = 'Rand';
			else
				vecUseTrialStart = vecTrialStart;
				strRand = 'Real';
			end
			matUseTrialT = cat(2,vecUseTrialStart',vecUseTrialStart'+dblStimDur);
			
			%run zeta & ts-zeta
			intPlot = 0;
			[dblTsZetaP,sTsZETA] = zetatstest(vecTimestamps,vecdFoF,vecUseTrialStart,dblTrialDur,intResampNum,intPlot);
			dblTsZeta = sTsZETA.dblZETA;
			matTsZeta(intNeuron,intRunType,intTrialDurIdx,intTauIdx,intSampFreqIdx) = dblTsZeta;
			%save
			if intTauIdx == intTauNum && intSampFreqIdx == intSampFreqNum
				[dblZetaP,sZETA] = zetatest(vecSpikeTimes,matUseTrialT,dblTrialDur,intResampNum,0);
				dblZeta = sZETA.dblZETA;
				dblMeanZ = sZETA.dblMeanZ;
				matZeta(intNeuron,intRunType,intTrialDurIdx) = dblZeta;
				matTtest(intNeuron,intRunType,intTrialDurIdx) = dblMeanZ;
			end

			end
		end
	end
end
end
%% save
if boolSave
	save([strDataTargetPath 'TsZetaFeatureConjunctionsResamp' num2str(intResampNum) '.mat' ],...
		'matZeta','matTtest','matTsZeta','vecTrialDur','vecTau','vecSampFreqs');
end

