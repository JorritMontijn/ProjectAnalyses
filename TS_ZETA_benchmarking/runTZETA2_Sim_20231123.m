%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

vecRandTypes = [1 2];
intResampNum = 500;
boolSave = true;%true;
dblUseDur = 8;
intNeurons = 100;
boolDirectQuantile = false;
intSuperResFactor = 100; %1 or 100
warning('off','zetatstest:InsufficientDataLength');

%% set variables
dblFracDiffSpikes = 1/2;
dblTau = 1;
dblTau0 = (63/1000);
dblNoise = 0.025;
dblSamplingFreq = 51.1;
dblSamplingInterval = 1/dblSamplingFreq;
boolQuick = false;

%set indicator properties
sIndicatorProps = struct;
sIndicatorProps.dblTimescale = dblTau;
sIndicatorProps.dblNoise = dblNoise;

%% load data
for intCompType=2%1:2
	if intCompType == 1
		strCompType = 'PeakHeight';
	else
		strCompType = 'PeakTime';
	end
	fprintf('Running %s [%s]\n',strCompType,getTime);
	
	%% pre-allocate output variables
	matTsZetaP = nan(2,intNeurons);
	matTtestP = nan(2,intNeurons);
	matAnovaP = nan(2,intNeurons);
	
	%% analyze
	hTic = tic;
	for intNeuron1=1:intNeurons%[1:intNeurons]%43%1:27, 2:69
		%% message
		if toc(hTic) > 5
			fprintf('Processing neuron %d/%d [%s]\n',intNeuron1,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		
		%% tuning params
		dblBaseRate = exprnd(1);
		boolDoublePeaked = false; %orientation or direction tuned
		dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
		dblKappa = rand(1)*5+5; %von Mises concentration parameter
		dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
		dblJitter = 50; %in ms'
		
		%% stimulus data
		dblStimDur = 3;
		dblPreBaseDur = 3;
		dblPostBaseDur = 2;
		dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
		dblUseMaxDur = dblStimDur;
		intOris = 24;
		dblStepDeg = 360/intOris;
		vecOris = linspace(0,360-dblStepDeg,intOris);
		intReps = 10;
		intTrials = intOris*intReps;
		vecTrialAngles1 = nan(1,intTrials);
		for intRep=1:intReps
			vecTrialAngles1(((intRep-1)*intOris+1):(intRep*intOris)) = vecOris(randperm(intOris));
		end
		vecTrialStart1 = 10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*intTrials+eps));
		matTrialT1 = cat(2,vecTrialStart1',vecTrialStart1'+dblStimDur);
		
		%unbalanced design
		vecTrialStart2 =  10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*(intTrials/2)+eps));
		matTrialT2 = cat(2,vecTrialStart2',vecTrialStart2'+dblStimDur);
		vecTrialAngles2 = vecTrialAngles1(1:size(matTrialT2,1));
		
		%% generate data
		% generate peak
		dblStartDelay = 0.2;
		dblPeakDelay1 = dblStartDelay;
		intAddSpikes1 = round(intTrials*dblFracDiffSpikes);
		[vecSpikeTimes1,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles1,matTrialT1,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes1,dblStartDelay,dblPeakDelay1);
		intTrials2 = size(matTrialT2,1);
		
		% generate dfof
		[vecTimestamps1,vecdFoF1] = getGeneratedFluorescence(vecSpikeTimes1,dblSamplingFreq,sIndicatorProps,boolQuick);
		%add empty end
		dblEndDur = vecTrialStart1(end) + dblUseMaxDur*5;
		vecAddT1 = vecTimestamps1(end):(1/dblSamplingFreq):dblEndDur;
		vecTimestamps1 = cat(2,vecTimestamps1,vecAddT1(2:end));
		vecdFoF1 = cat(2,vecdFoF1,zeros(size(vecAddT1(2:end))));
		
		for intRandType=vecRandTypes
			%% real+rand
			if intRandType == 1
				strRandType = '';
			elseif intRandType ==2
				strRandType = 'Rand';
			end
			%randomize
			if intRandType ==2
				%generate n2, no diff
				intAddSpikes2 = round((intAddSpikes1/intTrials)*intTrials2);
				dblPeakDelay2 = dblPeakDelay1;
				[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
			else
				if strcmp(strCompType,'PeakTime')
					%generate n2, diff in peak time
					dblPeakDelay2 = dblPeakDelay1+0.025;
					intAddSpikes2 = round(dblFracDiffSpikes*intTrials2);
					[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
				elseif strcmp(strCompType,'PeakHeight')
					%generate n2, diff in peak height
					dblPeakDelay2 = dblPeakDelay1;
					intAddSpikes2 = 2*round(dblFracDiffSpikes*intTrials2);
					[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
				else
					error('not recognized');
				end
			end
			
			%% generate time-series data
			% generate dfof2
			[vecTimestamps2,vecdFoF2] = getGeneratedFluorescence(vecSpikeTimes2,dblSamplingFreq,sIndicatorProps,boolQuick);
			%add empty end
			dblEndDur = vecTrialStart2(end) + dblUseMaxDur*5;
			vecAddT2 = vecTimestamps2(end):(1/dblSamplingFreq):dblEndDur;
			vecTimestamps2 = cat(2,vecTimestamps2,vecAddT2(2:end));
			vecdFoF2 = cat(2,vecdFoF2,zeros(size(vecAddT2(2:end))));
			
			%% run tests
			boolDirectQuantile=0;
			intPlot = 0;
			[dblZeta2P,sZETA] =  zetatstest2(vecTimestamps1,vecdFoF1,matTrialT1,vecTimestamps2,vecdFoF2,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
			if 0%sZETA.dblMeanP > 0.1 && dblZeta2P < 0.01
				intPlot = 4;
				[dblZeta2P,sZETA] = zetatstest2(vecTimestamps1,vecdFoF1,matTrialT1,vecTimestamps2,vecdFoF2,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
				subplot(2,3,5)
				bplot([sZETA.vecMu1' sZETA.vecMu2'],'outliers');
				ylabel('dF/F0');
				set(gca,'xtick',[1 2],'xticklabel',{'Pref stim','Orth stim'});
				title(sprintf('%s N%d',strRecIdx,intNeuron1));
				fixfig;
				pause
			end
			
			matTtestP(intRandType,intNeuron1) = sZETA.dblMeanP;
			matTsZetaP(intRandType,intNeuron1) = dblZeta2P;
			
			%ANOVA
			hTicA = tic;
			[vecRefT1,matTracePerTrial1] = getTraceInTrial(vecTimestamps1,vecdFoF1,matTrialT1(:,1),1/dblSamplingFreq,dblUseMaxDur);
			[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecTimestamps2,vecdFoF2,matTrialT2(:,1),1/dblSamplingFreq,dblUseMaxDur);
			
			vecBin1 = flat(repmat(1:size(matTracePerTrial1,2),[size(matTracePerTrial1,1) 1]));
			vecBin2 = flat(repmat(1:size(matTracePerTrial2,2),[size(matTracePerTrial2,1) 1]));
			vecR1 = matTracePerTrial1(:);
			vecR2 = matTracePerTrial2(:);
			
			%two-sample
			g1 = cat(1,ones(size(vecR1)),2*ones(size(vecR2)));
			g2 = cat(1,vecBin1,vecBin2);
			vecP=anovan(cat(1,vecR1,vecR2),{g1,g2},'model','interaction','display','off');%,'varnames',{'g1','g2'})%
			[h crit_p adj_p]=fdr_bh(vecP([1 3]));
			dblAnova2P = min(adj_p);
			
			%one-sample diff
			matAnovaP(intRandType,intNeuron1) = dblAnova2P;
			
		end
	end
	%% save data
	for intRandType=vecRandTypes
		%% real+rand
		vecAnovaP = matAnovaP(intRandType,:);
		vecTsZetaP = matTsZetaP(intRandType,:);
		vecTtestP = matTtestP(intRandType,:);
		if intRandType == 1
			strRandType = '';
		elseif intRandType ==2
			strRandType = 'Rand';
		end
		if boolSave
			save([strDataPath 'TsZeta2_' strCompType '_Q' num2str(boolDirectQuantile) '_' strRandType 'Resamp' num2str(intResampNum) '.mat' ],...
				'vecAnovaP','vecTsZetaP','vecTtestP','strRandType');
		end
	end
end
