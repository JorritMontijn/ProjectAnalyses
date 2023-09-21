%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;

if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta\';
	strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta\';
	strDataSourcePath = '';
	
end
strFigPath = [strPath 'Figs\'];
strDataTargetPath = [strPath 'Data\'];

vecRunTypes = [1 2];
intResampNum = 250;
boolSave = true;

%% set variables
strRec = 'SimSampFreqs';
intNeurons = 100;
vecSampFreqs = [1 5 25 100];
intSampFreqs = numel(vecSampFreqs);

%% pre-allocate output variables
matZeta = nan(intNeurons,2,intSampFreqs);
matTsZeta = nan(intNeurons,2,intSampFreqs);
hTic = tic;
%% generate data
for intNeuron=1:intNeurons
	%% message
	if toc(hTic) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTic=tic;
	end
	
	%% spiking params
	dblBaseRate = exprnd(1);
	dblBurstEventRate = exprnd(0.2);
	dblBurstDuration = normrnd(20,4);
	dblBurstISI = (0.5+exprnd(2.9))/1000;
	sSpikingParams.dblBaseRate = dblBaseRate;
	sSpikingParams.dblBurstEventRate = dblBurstEventRate;
	sSpikingParams.dblBurstDuration = min(max(10,dblBurstDuration),40);
	sSpikingParams.dblBurstISI = dblBurstISI;
	
	%% tuning params
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate+exprnd(1); %mean single-spike rate during stimulus (exponential ISI)
	dblPrefBurstEventRate = dblBurstEventRate+exprnd(1); %mean evoked rate of burst events (Hz) (exponential inter-event times)
	
	sTuningParams.boolDoublePeaked = boolDoublePeaked;
	sTuningParams.dblPrefOri = dblPrefOri;
	sTuningParams.dblKappa = dblKappa;
	sTuningParams.dblPrefRate = dblPrefRate;
	sTuningParams.dblPrefBurstEventRate = dblPrefBurstEventRate;
	
	%% stimulus data
	dblStimDur = 3;
	dblPreBaseDur = 3;
	dblPostBaseDur = 2;
	dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
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
	% generate bursts
	[vecSpikeTimes1,dblPrefOri] = getGeneratedBurstingData(vecTrialAngles,matTrialT,sSpikingParams,sTuningParams);
	
	% generate peak
	dblJitter = 5; %in ms'
	[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,0,0,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri);
	
	%combine
	vecSpikeTimes = sort(cat(1,vecSpikeTimes1,vecSpikeTimes2));
	
	%% go through sampling freqs
	for intSampFreqIdx=1:intSampFreqs
		dblSamplingFreq = vecSampFreqs(intSampFreqIdx);
		
		fprintf('Running samp freq %d/%d [%s]\n',intSampFreqIdx,intSampFreqs,getTime);
		
		%set indicator properties
		sIndicatorProps = struct;
		sIndicatorProps.dblTimescale = 0;
		sIndicatorProps.dblNoise = 0;
		[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps);
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
				vecUseTrialStart = vecTrialStart + vecJitter -3;
			else
				vecUseTrialStart = vecTrialStart -3;
			end
			%plot if first
			if intNeuron == 1
				intPlot = 0;%2;
			else
				intPlot = 0;
			end
			%run zeta & ts-zeta
			intResampNum = 100;
			[dblTsZetaP,sTsZETA] = zetatstest(vecTimestamps,vecdFoF,vecUseTrialStart,dblTrialDur,intResampNum,intPlot);
			dblTsZeta = sTsZETA.dblZETA;
			matTsZeta(intNeuron,intRunType,intSampFreqIdx) = dblTsZeta;
			
			%save
			if 1%intSampFreqIdx == 1
				[dblZetaP,sZETA] = zetatest(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,intPlot);
				dblZeta = sZETA.dblZETA;
				matZeta(intNeuron,intRunType,intSampFreqIdx) = dblZeta;
			end
		end
	end
end

%% save
if boolSave
	save([strDataTargetPath 'TsZeta' strRec 'Resamp' num2str(intResampNum) '.mat' ],...
		'matZeta','matTsZeta','vecSampFreqs','strRec');
end