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
strRec = 'SimTwoSampleZeta';
intNeurons = 1000;
intFracDiffSpikes = 0.5;
	
%% pre-allocate output variables
matTtest2 = nan(intNeurons,2,intResamps);
matZeta2 = nan(intNeurons,2,intResamps);

%% generate data
hTicN = tic;
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
	end
	
	%% tuning params
	dblBaseRate = exprnd(5);
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	dblJitter = 5; %in ms'
	
	%% stimulus data
	dblStimDur = 3;
	dblPreBaseDur = 3;
	dblPostBaseDur = 2;
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
	matTrialT = cat(2,vecTrialStart',vecTrialStart'+dblStimDur);
	
	%% generate data
	% generate bursts
	intAddSpikes1 = intTrials/2;
	intDiffSpikes = round(intFracDiffSpikes*intTrials);
	
	% generate peak
	[vecSpikeTimes1,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes1);
	
	%real+rand
	for intRunType=vecRunTypes
		
		%randomize
		if intRunType ==2
			%generate n2, no diff
			intAddSpikes2 = intAddSpikes1;
			[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2);
		else
			%generate n2, diff
			intAddSpikes2 = intAddSpikes1 + intDiffSpikes;
			[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2);
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
			
			[dblZetaP,sZETA] = zetatest2(vecSpikeTimes1,matTrialT,vecSpikeTimes2,matTrialT,dblUseTrialDur,intResampNum,intPlot);
			matTtest2(intNeuron,intRunType,intResampIdx) = -norminv(sZETA.dblMeanP/2);
			matZeta2(intNeuron,intRunType,intResampIdx) = -norminv(dblZetaP/2);
		end
	end
end

%% save
if boolSave
	save([strDataTargetPath 'Zeta2' strRec '.mat' ],...
		'matTtest2','matZeta2','strRec');
end
