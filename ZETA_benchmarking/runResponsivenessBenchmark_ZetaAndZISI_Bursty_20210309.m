%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;

strDataMasterPath = 'F:\Data\Processed\ePhys\';
strDataTargetPath = 'F:\Data\Processed\ZETA\Inclusion\';
strFigPath = 'F:\Data\Results\ZETA\Examples\';
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 300;%100;%10:10:90;%[10:10:100];
intNeurons = 1000;

for intRandType=vecRandTypes
	%reset vars
	clearvars -except intRandType intNeurons strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRandTypes vecRestrictRange boolSave vecResamples
	strArea = 'BurstyTest';
	
	if intRandType == 1
		strRunType = strArea;
		fprintf('Prepping normal... [%s]\n',getTime);
	elseif intRandType ==2
		strRunType = [strArea '-Rand'];
		fprintf('Prepping random... [%s]\n',getTime);
	end
	
	for intResampleIdx = 1:numel(vecResamples)
		intResampleNum = vecResamples(intResampleIdx);
		%% message
		fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
		hTicMessage=tic;
		
		%% pre-allocate output variables
		vecNumSpikes = nan(1,intNeurons);
		vecZetaP = nan(1,intNeurons);
		vecHzP = nan(1,intNeurons);
		
		vecPoissP = nan(1,intNeurons);
		
		vecISIP_ks = nan(1,intNeurons);
		vecISIP_IntG = nan(1,intNeurons);
		vecISIP_g = nan(1,intNeurons);
		
		vecBISIP_g = nan(1,intNeurons);
		vecBISIP_ks = nan(1,intNeurons);
		
		vecComputTimeZETA = nan(1,intNeurons);
		vecComputTimePoiss = nan(1,intNeurons);
		vecComputTimeISI = nan(1,intNeurons);
		vecComputTimeBISI = nan(1,intNeurons);
		
		
		cellArea = cell(1,intNeurons);
		
		%% analyze
		for intNeuron=[1:intNeurons]%31
			%% generate
			strRecIdx = 'x';
			strMouse = 'Artificial';
			strBlock = '1';
			strDate = getDate();
			intSU = intNeuron;
			intClust = intNeuron;
			
			%% set parameters
			vecUniqueTrialAngles = 0:(360/24):359;
			vecTrialAngles = [];
			for intRep=1:20
				vecTrialAngles = cat(2,vecTrialAngles,vecUniqueTrialAngles(randperm(numel(vecUniqueTrialAngles))));
			end
			vecTrialAngles = deg2rad(vecTrialAngles);
			
			intT = numel(vecTrialAngles);
			dblDur = 1;
			dblITI = 0.5;
			dblPreLead=1;
			vecStarts = dblPreLead:(dblDur+dblITI):((intT-0.5)*(dblDur+dblITI)+dblPreLead);
			vecStops = vecStarts + dblDur;
			matTrialT = cat(2,vecStarts',vecStops');
			
			% spiking params
			sSpikingParams = struct;
			sSpikingParams.dblBaseRate = exprnd(1); %mean baseline single spike rate (Hz) (exponential ISI)
			sSpikingParams.dblBurstEventRate = (abs(normrnd(0.05,0.0125))+0.0125); %mean baseline rate of burst events (Hz) (exponential inter-event times)
			sSpikingParams.dblBurstDuration = normrnd(100,10); %mean duration of burst events (ms) (Gamma, theta=0.5;k=2*dblBurstDuration)
			sSpikingParams.dblBurstISI = (0.5+exprnd(2.4))/1000; %mean ISI during bursts (s) (Gamma, theta=0.5;k=2*dblBurstISI)
			
			% tuning params
			sTuningParams = struct;
			sTuningParams.boolDoublePeaked = false; %orientation or direction tuned
			sTuningParams.dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
			sTuningParams.dblKappa = rand(1)*5+5; %von Mises concentration parameter
			sTuningParams.dblPrefRate = sSpikingParams.dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
			
			%% get visual responsiveness
			%get trial dur
			dblUseMaxDur = round(median(diff(matTrialT(:,1)))*2)/2;
			%set derivative params
			if contains(strRunType,'Rand')
				sSpikingParams.dblBurstEventRate = 2.35*sSpikingParams.dblBurstEventRate; %mean baseline rate of burst events (Hz) (exponential inter-event times)
				sTuningParams.dblPrefBurstEventRate = 0; %mean evoked rate of burst events (Hz) (exponential inter-event times)
				[vecSpikeTimes,dblPrefOri] = getGeneratedBurstingData(vecTrialAngles,matTrialT,sSpikingParams,sTuningParams);
			else
				sTuningParams.dblPrefBurstEventRate = 20*sSpikingParams.dblBurstEventRate; %mean evoked rate of burst events (Hz) (exponential inter-event times)
				[vecSpikeTimes,dblPrefOri] = getGeneratedBurstingData(vecTrialAngles,matTrialT,sSpikingParams,sTuningParams);
			end
			matEventTimes = matTrialT;
			%vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
			
			%pre-alloc
			sZETA = [];
			sZETA.dblMeanP = 1;
			dblZetaP = 1;
			dblPoissP = 1;
			dblP_KS = 1;
			dblP_G = 1;
			dblIntP_G = 1;
			dblP_BISI = 1;
			dblP_BISI_KS = 1;
			dblComputTimePoiss = 0;
			dblComputTimeISI = 0;
			
			%ZETA
			hTic=tic;
			intPlot = 0;
			close
			[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,intPlot,0);
			dblComputTimeZETA = toc(hTic);
			
			%Poisson
			hTic=tic;
			%dblPoissP = getPoissonTest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
			dblComputTimePoiss= toc(hTic);
			
			%ISI varieties
			hTic = tic;
			%[dblP_KS,dblP_G,dblIntP_G] = getISItest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
			dblComputTimeISI = toc(hTic);
			
			%BISI
			hTic = tic;
			[dblP_BISI,dblP_BISI_KS] = getBISI(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
			dblComputTimeBISI = toc(hTic);
			
			intSpikeNum = numel(vecSpikeTimes);
			
			% assign data
			vecNumSpikes(intNeuron) = intSpikeNum;
			vecZetaP(intNeuron) = dblZetaP;
			vecHzP(intNeuron) = sZETA.dblMeanP;
			
			vecPoissP(intNeuron) = dblPoissP;
			
			vecISIP_ks(intNeuron) = dblP_KS;
			vecISIP_IntG(intNeuron) = dblIntP_G;
			vecISIP_g(intNeuron) = dblP_G;
			
			vecBISIP_g(intNeuron) = dblP_BISI;
			vecBISIP_ks(intNeuron) = dblP_BISI_KS;
			
			cellArea{intNeuron} = strArea;
			vecComputTimeZETA(intNeuron) = dblComputTimeZETA;
			vecComputTimePoiss(intNeuron) = dblComputTimePoiss;
			vecComputTimeISI(intNeuron) = dblComputTimeISI;
			vecComputTimeBISI(intNeuron) = dblComputTimeBISI;
			
			%% message
			%if toc(hTicMessage) > 5 && intNeuron > 1
			fprintf('Processed neuron %d/%d, Compute time: ZETA=%.1fs, Poiss=%.1fs, ISI=%.1fs, BISI=%.1fs [%s]\n',intNeuron,intNeurons,...
				vecComputTimeZETA(intNeuron),...
				vecComputTimePoiss(intNeuron),...
				vecComputTimeISI(intNeuron),...
				vecComputTimeBISI(intNeuron),getTime);
			%hTicMessage=tic;
			%end
			%clear vecTrialStarts;
			
		end
		%save
		if boolSave
			save([strDataTargetPath 'ZetaBurst' strRunType 'Resamp' num2str(intResampleNum) '.mat' ],...
				'vecNumSpikes','cellArea',...
				'vecZetaP',...
				'vecHzP',...
				'vecPoissP',...
				'vecISIP_ks',...
				'vecISIP_IntG',...
				'vecISIP_g',...
				'vecBISIP_g',...
				'vecBISIP_ks',...
				'vecComputTimeZETA',...
				'vecComputTimePoiss',...
				'vecComputTimeISI',...
				'vecComputTimeBISI');
		end
	end
end