%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
intResampNum = 10000;
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = 'SimZeta2';
intNeurons = 100;
%% pre-allocate output variables
matZeta = nan(intNeurons,2);
matZetaTime = nan(intNeurons,2);
		
matZeta2 = nan(intNeurons,2);
matZeta2Time = nan(intNeurons,2);
		
hTicN = tic;
%% generate data
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
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
			intPlot = 0;
		else
			intPlot = 0;
		end
		%run zeta & zeta2
		hTic=tic;
		[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,intPlot);
		dblZeta = sZETA.dblZETA;
		matZeta(intNeuron,intRunType) = dblZeta;
		matZetaTime(intNeuron,intRunType) = toc(hTic);
		
		%save
		hTic2=tic;
		[dblZeta2P,vecLatencies,sZETA2] = getZeta2(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,intPlot);
		dblZeta2 = sZETA2.dblZETA;
		matZeta2(intNeuron,intRunType) = dblZeta2;
		matZeta2Time(intNeuron,intRunType) = toc(hTic2);
	end
end

%% save
if boolSave
	save([strDataTargetPath 'Zeta2' strRec 'Resamp' num2str(intResampNum) '.mat' ],...
		'matZeta','matZeta2','matZetaTime','matZeta2Time','strRec');
end