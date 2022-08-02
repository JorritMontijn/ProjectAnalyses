%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
vecResamps = 1000;%[100 200 500 1000 2000];
intResamps= numel(vecResamps);
%vecResamps = [5000 10000];
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = 'SimTwoSampleTsZeta';
intNeurons = 100;
dblFracDiffSpikes = 1/2;
dblTau = 2;
dblTau0 = (63/1000);
dblNoise = 0.025;
dblSamplingFreq = 25;
boolQuick = false;

%set indicator properties
sIndicatorProps = struct;
sIndicatorProps.dblTimescale = dblTau;
sIndicatorProps.dblNoise = dblNoise;

%% pre-allocate output variables
matTtest2 = nan(intNeurons,2,intResamps);
matTsZeta2 = nan(intNeurons,2,intResamps);
matAnova2 = nan(intNeurons,2,intResamps);

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
	
	intReps2 = 12;
	intTrials2 = intOris*intReps2;
	vecTrialAngles2 = nan(1,intTrials2);
	for intRep=1:intReps2
		vecTrialAngles2(((intRep-1)*intOris+1):(intRep*intOris)) = vecOris(randperm(intOris));
	end
	vecTrialStart2 = 10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*intTrials2+eps));
	matTrialT2 = cat(2,vecTrialStart2',vecTrialStart2'+dblStimDur);
	matTrialT2 = matTrialT2 + 0.05*rand(size(matTrialT2));
	
	%% generate data
	% generate bursts
	intAddSpikes1 = intTrials/2;
	intDiffSpikes = round(dblFracDiffSpikes*intTrials);
	
	% generate peak
	dblStartDelay = 0.1;
	[vecSpikeTimes1,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT1,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes1,dblStartDelay);
	[vecTimestamps1,vecdFoF1] = getGeneratedFluorescence(vecSpikeTimes1,dblSamplingFreq,sIndicatorProps,boolQuick);
	%add empty end
	dblEndDur = vecTrialStart(end) + dblTrialDur*5;
	vecAddT = vecTimestamps1(end):(1/dblSamplingFreq):dblEndDur;
	vecTimestamps1 = cat(2,vecTimestamps1,vecAddT(2:end));
	vecdFoF1 = cat(2,vecdFoF1,zeros(size(vecAddT(2:end))));
	
	%real+rand
	for intRunType=vecRunTypes
		%randomize
		if intRunType ==2
			%generate n2, no diff
			intAddSpikes2 = intAddSpikes1;
			dblStartDelay2 = dblStartDelay;
			[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay2);
		else
			%generate n2, diff
			intAddSpikes2 = intAddSpikes1 + intDiffSpikes;
			dblStartDelay2 = dblStartDelay;
			[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay2);
		end
		
		%make dF/F
		[vecTimestamps2,vecdFoF2] = getGeneratedFluorescence(vecSpikeTimes2,dblSamplingFreq,sIndicatorProps,boolQuick);
		%add empty end
		dblEndDur = vecTrialStart2(end) + dblTrialDur*5;
		vecAddT = vecTimestamps2(end):(1/dblSamplingFreq):dblEndDur;
		vecTimestamps2 = cat(2,vecTimestamps2,vecAddT(2:end));
		vecdFoF2 = cat(2,vecdFoF2,zeros(size(vecAddT(2:end))));
		intS = min([numel(vecTimestamps1) numel(vecTimestamps2)]);
		
		%run zeta & zeta2
		for intResampIdx = 1:intResamps
			intResampNum = vecResamps(intResampIdx);
			intPlot = 0;
			boolPairwise =0;
			boolDirectQuantile=0;
			%[dblZetaP,sZETA] = zetatstest(vecTimestamps1(1:intS),vecdFoF1(1:intS)-vecdFoF2(1:intS),matTrialT,dblUseTrialDur,intResampNum,0,boolDirectQuantile);
			[dblZeta2P,sZETA] = zetatstest2(vecTimestamps1,vecdFoF1,matTrialT1,vecTimestamps2,vecdFoF2,matTrialT2,vecUseDur,intResampNum,intPlot,boolPairwise,boolDirectQuantile);
			matTtest2(intNeuron,intRunType,intResampIdx) = -norminv(sZETA.dblMeanP/2);
			matTsZeta2(intNeuron,intRunType,intResampIdx) = -norminv(dblZeta2P/2);
			
			%ANOVA
			hTicA = tic;
			[vecRefT1,matTracePerTrial1] = getTraceInTrial(vecTimestamps1,vecdFoF1,matTrialT1(:,1),1/dblSamplingFreq,dblUseTrialDur);
			[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecTimestamps2,vecdFoF2,matTrialT2(:,1),1/dblSamplingFreq,dblUseTrialDur);
			
			vecBin1 = flat(repmat(1:size(matTracePerTrial1,2),[size(matTracePerTrial1,1) 1]));
			vecBin2 = flat(repmat(1:size(matTracePerTrial2,2),[size(matTracePerTrial2,1) 1]));
			vecR1 = matTracePerTrial1(:);
			vecR2 = matTracePerTrial2(:);
			
			%two-sample
			g1 = cat(1,ones(size(vecR1)),2*ones(size(vecR2)));
			g2 = cat(1,vecBin1,vecBin2);
			vecP=anovan(cat(1,vecR1,vecR2),{g1,g2},'model','interaction','display','off');%,'varnames',{'g1','g2'})%
			dblAnova2P = vecP(3);
			
			%one-sample diff
			dblAnovaDur = toc(hTicA);
			matAnova2(intNeuron,intRunType,intResampIdx) = -norminv(dblAnova2P/2);
			
			if 0
				%% plot
				subplot(2,3,1)
				imagesc([],vecRefT1,matTracePerTrial1);
				
				subplot(2,3,2)
				imagesc([],vecRefT2,matTracePerTrial2);
			end
		end
	end
end

%% save
if boolSave
	save([strDataTargetPath strRec 'Q' num2str(boolDirectQuantile) '.mat' ],...
		'matAnova2','matTtest2','matTsZeta2','strRec');
end
