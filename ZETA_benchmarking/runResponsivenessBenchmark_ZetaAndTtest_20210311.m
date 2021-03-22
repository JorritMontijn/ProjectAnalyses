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
vecResamples = 2.^(1:7);%2.^(3:10);%100;%10:10:90;%[10:10:100];
%difference in spiking
vecSpikingDifference=2.^(-5:3);
intNeurons = 1000;
hTicMessage=tic;
for intRandType=vecRandTypes
	%reset vars
	clearvars -except hTicMessage vecSpikingDifference intRandType intNeurons strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRandTypes vecRestrictRange boolSave vecResamples
	strArea = 'PoissonSpiking';
	
	if intRandType == 1
		strRunType = strArea;
		fprintf('Prepping normal... [%s]\n',getTime);
	elseif intRandType ==2
		strRunType = [strArea '-Rand'];
		fprintf('Prepping random... [%s]\n',getTime);
	end
	
	for intResampleIdx = 1:numel(vecResamples)
		intResampleNum = vecResamples(intResampleIdx);
		
		%% pre-allocate output variables
		intSpDiffNum = numel(vecSpikingDifference);
		vecSpDiff = nan(intSpDiffNum,intNeurons);
		vecNumSpikes = nan(intSpDiffNum,intNeurons);
		vecZetaP = nan(intSpDiffNum,intNeurons);
		vecHzP = nan(intSpDiffNum,intNeurons);
		
		vecPoissP = nan(intSpDiffNum,intNeurons);
		
		vecISIP_ks = nan(intSpDiffNum,intNeurons);
		vecISIP_IntG = nan(intSpDiffNum,intNeurons);
		vecISIP_g = nan(intSpDiffNum,intNeurons);
		
		vecBISIP_g = nan(intSpDiffNum,intNeurons);
		vecBISIP_ks = nan(intSpDiffNum,intNeurons);
		
		vecComputTimeZETA = nan(intSpDiffNum,intNeurons);
		vecComputTimePoiss = nan(intSpDiffNum,intNeurons);
		vecComputTimeISI = nan(intSpDiffNum,intNeurons);
		vecComputTimeBISI = nan(intSpDiffNum,intNeurons);
		
		for intSpikingDiffIdx=1:intSpDiffNum
			dblSpikingDifference = vecSpikingDifference(intSpikingDiffIdx);
			%% message
			fprintf('Processing %s, resampling %d (%d/%d), spiking diff %.3f (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),...
				dblSpikingDifference,intSpikingDiffIdx,intSpDiffNum,getTime);
			
			
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
				sSpikingParams.dblBaseRate = 100; %mean baseline single spike rate (Hz) (exponential ISI)
				sSpikingParams.dblBurstEventRate = 0; %mean baseline rate of burst events (Hz) (exponential inter-event times)
				sSpikingParams.dblBurstDuration = 900; %mean duration of burst events (ms) (Gamma, theta=0.5;k=2*dblBurstDuration)
				sSpikingParams.dblBurstISI = 100/1000; %mean ISI during bursts (ms) (Gamma, theta=0.5;k=2*dblBurstISI)
				
				% tuning params
				sTuningParams = struct;
				sTuningParams.boolDoublePeaked = false; %orientation or direction tuned
				sTuningParams.dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
				sTuningParams.dblKappa = 0.1; %von Mises concentration parameter
				sTuningParams.dblPrefRate = 100; %mean single-spike rate during stimulus (exponential ISI)
				sTuningParams.dblPrefBurstEventRate = 1000; %mean evoked rate of burst events (Hz) (exponential inter-event times)
				
				%% get visual responsiveness
				matGenT = matTrialT;
				dblBase=1;
				dblResp=dblBase+dblSpikingDifference;
				matGenT(:,2) =matTrialT(:,1)+0.1;
				[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingData(vecTrialAngles,matTrialT,dblBase,dblResp,0.1,true);
				%[vecSpikeTimes,dblPrefOri] = getGeneratedBurstingData(vecTrialAngles,matGenT,sSpikingParams,sTuningParams);
				%set derivative params
				if contains(strRunType,'Rand')
					vecJitter = 4*dblDur*rand([numel(matTrialT(:,1)) 1])-dblDur;
					matEventTimes = bsxfun(@plus,matTrialT,vecJitter);
				else
					matEventTimes = matTrialT;
				end
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
				dblUseMaxDur = dblDur+dblITI;
				hTic=tic;
				intPlot = 0;%2;
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
				%[dblP_BISI,dblP_BISI_KS] = getBISI(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
				dblComputTimeBISI = toc(hTic);
				
				intSpikeNum = numel(vecSpikeTimes);
				
				% assign data
				vecSpDiff(intSpikingDiffIdx,intNeuron) = dblSpikingDifference;
				vecNumSpikes(intSpikingDiffIdx,intNeuron) = intSpikeNum;
				vecZetaP(intSpikingDiffIdx,intNeuron) = dblZetaP;
				vecHzP(intSpikingDiffIdx,intNeuron) = sZETA.dblMeanP;
				
				vecPoissP(intSpikingDiffIdx,intNeuron) = dblPoissP;
				
				vecISIP_ks(intSpikingDiffIdx,intNeuron) = dblP_KS;
				vecISIP_IntG(intSpikingDiffIdx,intNeuron) = dblIntP_G;
				vecISIP_g(intSpikingDiffIdx,intNeuron) = dblP_G;
				
				vecBISIP_g(intSpikingDiffIdx,intNeuron) = dblP_BISI;
				vecBISIP_ks(intSpikingDiffIdx,intNeuron) = dblP_BISI_KS;
				
				cellArea{intSpikingDiffIdx,intNeuron} = strArea;
				vecComputTimeZETA(intSpikingDiffIdx,intNeuron) = dblComputTimeZETA;
				vecComputTimePoiss(intSpikingDiffIdx,intNeuron) = dblComputTimePoiss;
				vecComputTimeISI(intSpikingDiffIdx,intNeuron) = dblComputTimeISI;
				vecComputTimeBISI(intSpikingDiffIdx,intNeuron) = dblComputTimeBISI;
				
				%% message
				if toc(hTicMessage) > 5 && intNeuron > 1
				fprintf('Processed neuron %d/%d, Compute time: ZETA=%.1fs, Poiss=%.1fs, ISI=%.1fs, BISI=%.1fs [%s]\n',intNeuron,intNeurons,...
					vecComputTimeZETA(intSpikingDiffIdx,intNeuron),...
					vecComputTimePoiss(intSpikingDiffIdx,intNeuron),...
					vecComputTimeISI(intSpikingDiffIdx,intNeuron),...
					vecComputTimeBISI(intSpikingDiffIdx,intNeuron),getTime);
				hTicMessage=tic;
				end
				%clear vecTrialStarts;
				
			end
			
		end
		%save
			if boolSave
				save([strDataTargetPath 'ZetaTtestSpDResN' strRunType 'Resamp' num2str(intResampleNum) '.mat' ],...
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