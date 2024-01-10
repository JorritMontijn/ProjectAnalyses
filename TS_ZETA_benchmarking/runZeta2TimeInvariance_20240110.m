%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'HeteroPoissonPeak',...Area 1
	'TriPhasic',...Area 2
	'QuadriPhasic',...Area 3
	'iidGaussian',...Area 4
	'TriHomoPhasic',...Area 5
	'',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	};

if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\DataPreProcessed\';
end
strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
boolSave = true;
intResampleNum = 1000;%250;%10:10:90;%[10:10:100];
intNeurons = 1000;%10000
intTrialNum = 100;
vecDurs = [1 0.1 0.8 0.1];
dblUseMaxDur = sum(vecDurs);
vecTimeShifts = (-dblUseMaxDur/2):0.25:(dblUseMaxDur/2);
intTimeShiftNum = numel(vecTimeShifts);

cellNeuron = cell(intTimeShiftNum,intNeurons,2);
matNumSpikes = nan(intTimeShiftNum,intNeurons,2);
matZeta2P = nan(intTimeShiftNum,intNeurons,2);
matTtest2P = nan(intTimeShiftNum,intNeurons,2);
matKs2P = nan(intTimeShiftNum,intNeurons,2);
%set var
hTic=tic;
for intRandType=vecRandTypes
	%reset vars
	
	if intRandType == 1
		strRunType = 'TI2Z';
		fprintf('Prepping normal... [%s]\n',getTime);
	elseif intRandType ==2
		strRunType = 'TI2Z-Rand';
		fprintf('Prepping random... [%s]\n',getTime);
	end
	
	%% analyze
	for intNeuron=1:intNeurons%26
		%% message
		if toc(hTic) > 5
			
			%% message
			fprintf('Processing %s, neuron %d/%d [%s]\n',strRunType,intNeuron,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		
		%% prep data
		strDate = getDate();
		intSU = intNeuron;
		
		%neuron 1
		dblBaseRate = 0.1;
		dblFirstRate1 = exprnd(3);
		dblFirstRate2 = exprnd(6);
		dblSecondRate = exprnd(4);
		vecRates1 = [0 dblFirstRate1+5 1 dblSecondRate];
		vecRates2 = [0 dblFirstRate2+5 1 dblSecondRate];
		intNumT = intTrialNum+2;
		vecRepStarts = 5+(dblUseMaxDur:dblUseMaxDur:(dblUseMaxDur*intNumT));
		dblEndT = vecRepStarts(end)+5;
		vecSpikeTimes1 = getGeneratedMultiPhasicR(dblBaseRate,vecRates1,vecDurs,vecRepStarts,dblEndT);
		if contains(strRunType,'Rand')
			%generate duplicate
			vecSpikeTimes2 = getGeneratedMultiPhasicR(dblBaseRate,vecRates1,vecDurs,vecRepStarts,dblEndT);
		else
			%generate difference
			vecSpikeTimes2 = getGeneratedMultiPhasicR(dblBaseRate,vecRates2,vecDurs,vecRepStarts,dblEndT);
		end
		matEventTimes = vecRepStarts(:);
		matEventTimes(:,2) = matEventTimes(:,1) + dblUseMaxDur;
		matEventTimes([1 end],:) = [];
		
		%% get visual responsiveness
		%get trial dur
		%set derivative params
		intTrials = size(matEventTimes,1);
		intSpikeNum = numel(vecSpikeTimes1)+numel(vecSpikeTimes2);
		[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes1,matEventTimes(:,1),dblUseMaxDur);
		[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes2,matEventTimes(:,1),dblUseMaxDur);
		if intSpikeNum<3 || numel(vecTimePerSpike1) < 1|| numel(vecTimePerSpike2)< 1,continue;end
		
		%% run time shifts
		for intTimeShiftIdx=1:intTimeShiftNum
			dblTimeShift = vecTimeShifts(intTimeShiftIdx);
			intPlot = 4;
			boolDirectQuantile = false;
			[dblZetaP,sZeta2]=zetatest2(vecSpikeTimes1,matEventTimes+dblTimeShift,vecSpikeTimes2,matEventTimes+dblTimeShift,dblUseMaxDur,intResampleNum,intPlot,boolDirectQuantile);
			dblTtestP = sZeta2.dblMeanP;
			
			%save
drawnow;
export_fig(fullpath(strFigPath,sprintf('ZETA2_TimeShiftWithMuS%.2f.jpg',dblTimeShift)));
export_fig(fullpath(strFigPath,sprintf('ZETA2_TimeShiftWithMuS%.2f.pdf',dblTimeShift)));

			%%
			% assign data
			cellNeuron{intTimeShiftIdx,intNeuron,intRandType} = [strDate 'N' num2str(intSU)];
			matNumSpikes(intTimeShiftIdx,intNeuron,intRandType) = intSpikeNum;
			matZeta2P(intTimeShiftIdx,intNeuron,intRandType) = dblZetaP;
			matTtest2P(intTimeShiftIdx,intNeuron,intRandType) = dblTtestP;
		end
		return
	end
end

if boolSave
	save([strDataTargetPath 'Zeta2TimeInvariance.mat' ],...
		'cellNeuron','matNumSpikes','matTtest2P','matZeta2P','vecTimeShifts');
end