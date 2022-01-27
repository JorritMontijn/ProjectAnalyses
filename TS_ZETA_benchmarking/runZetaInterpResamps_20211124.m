%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
vecResamps = 100;%[100 200 500 1000 2000];
intResamps= numel(vecResamps);
%vecResamps = [5000 10000];
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = 'SimZeta3';
intNeurons = 1000;
%% pre-allocate output variables
matZeta = nan(intNeurons,2,intResamps);
matZetaTime = nan(intNeurons,2,intResamps);

matZeta2 = nan(intNeurons,2,intResamps);
matZeta2Time = nan(intNeurons,2,intResamps);

%% generate data
hTicN = tic;
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
	dblUseTrialDur = 4.5;
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
	
	%generate offset
	matTrialT2 = cat(2,matTrialT(:,1)+dblStimDur+0.1,matTrialT(:,1)+dblStimDur+dblPostBaseDur-0.1);
	[vecSpikeTimes3,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT2,0,100,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri);
	
	%combine
	vecSpikeTimes = sort(cat(1,vecSpikeTimes1,vecSpikeTimes2,vecSpikeTimes3));
	
	%real+rand
	for intRunType=vecRunTypes
		
		%randomize
		if intRunType ==2
			vecJitter = 2*dblUseTrialDur*((rand(size(vecTrialStart))-1/2)*2);
			vecUseTrialStart = vecTrialStart + vecJitter -3;
		else
			vecUseTrialStart = vecTrialStart -3;
		end
		%plot if first
		if intNeuron == 1
			%f = @() getZeta(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,0);
			%t=timeit(f);
			%f2 = @() getZeta2(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,0);
			%t2=timeit(f2);
			intPlot = 0;
		else
			intPlot = 0;
		end
		
		%run zeta & zeta2
		for intResampIdx = 1:intResamps
			intResampNum = vecResamps(intResampIdx);
			
			hTic1=tic;
			dblZetaP = getZeta(vecSpikeTimes,vecUseTrialStart,dblUseTrialDur,intResampNum,intPlot);
			dblZeta = -norminv(dblZetaP/2);
			matZeta(intNeuron,intRunType,intResampIdx) = dblZeta;
			matZetaTime(intNeuron,intRunType,intResampIdx) = toc(hTic1);
			
			%save
			hTic2=tic;
			dblZeta2P = zetatest(vecSpikeTimes,vecUseTrialStart,dblUseTrialDur,intResampNum,intPlot);
			dblZeta2 = -norminv(dblZeta2P/2);
			matZeta2(intNeuron,intRunType,intResampIdx) = dblZeta2;
			matZeta2Time(intNeuron,intRunType,intResampIdx) = toc(hTic2);
		end
	end
end

%% save
if boolSave
	save([strDataTargetPath 'Zeta2' strRec 'Resamps.mat' ],...
		'matZeta','matZeta2','matZetaTime','matZeta2Time','strRec','vecResamps','sSpikingParams','sTuningParams');
end
