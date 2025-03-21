%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;

vecRunTypes = [1 2];
intResampNum = 500;%[100 200 500 1000 2000];
boolSave = true;

if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta\';
end
strFigPath = [strPath 'Figs\'];
strDataTargetPath = [strPath 'Data\'];

%% set variables
intNeurons = 10000;
dblFracDiffSpikes = 1/2;
dblTau = 2;
dblTau0 = (63/1000);
dblNoise = 0.025;
dblSamplingFreq = 25;
dblSamplingInterval = 1/dblSamplingFreq;
boolQuick = false;
strRec = sprintf('TsZetaQuadriPhasicN%dR%d',intNeurons,intResampNum);
boolDirectQuantile = false;
        
%set indicator properties
sIndicatorProps = struct;
sIndicatorProps.dblTimescale = dblTau;
sIndicatorProps.dblNoise = dblNoise;

%% pre-allocate output variables
matTtest = nan(intNeurons,2);
matTsZeta = nan(intNeurons,2);
matAnova = nan(intNeurons,2);
matZeta_sp = nan(intNeurons,2);

%% generate data
hTicN = tic;
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
	end
	
	%% stimulus data
	dblFactor1=3;
	dblFactor2=3;
	dblBaseRate = exprnd(0.1)+0.1;
	vecDurs = [0.1 0.9 0.1];
	vecRates = [dblFactor1*exprnd(1)+dblFactor1/10 exprnd(0.2)+0.1 dblFactor2*exprnd(1)+dblFactor2/10];
	intNumT = 100;
	vecTrialDur=linspace(0.5,10,intNumT);
	vecRepStarts = 5+cumsum(vecTrialDur);
	dblEndT = vecRepStarts(end)+5;
	vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecRepStarts,dblEndT);
	
	dblUseMaxDur = 1;
	vecTrialStarts = vecRepStarts(:);
	matTrialT1 = cat(2,vecTrialStarts,vecTrialStarts+dblUseMaxDur);
	matTrialT1 = matTrialT1 + 0.05*rand(size(matTrialT1));
	
	%% generate data
	% generate dfof
	[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick);
	%add empty end
	dblEndDur = vecTrialStarts(end) + dblUseMaxDur*5;
	vecAddT = vecTimestamps(end):(1/dblSamplingFreq):dblEndDur;
	vecTimestamps = cat(2,vecTimestamps,vecAddT(2:end));
	vecdFoF = cat(2,vecdFoF,zeros(size(vecAddT(2:end))));
	
	%real+rand
	for intRunType=vecRunTypes
		%% get visual responsiveness
		%set derivative params
		if intRunType ==2
			vecJitterPerTrial = 2*linspace(dblUseMaxDur/intNumT,dblUseMaxDur,intNumT)';
			matEventTimes = bsxfun(@plus,matTrialT1,vecJitterPerTrial(randperm(numel(vecJitterPerTrial))));
		else
			matEventTimes = matTrialT1;
		end
		
		%ANOVA
		hTicA = tic;
		[vecRefT2,matTracePerTrial] = getTraceInTrial(vecTimestamps,vecdFoF,matEventTimes(:,1),dblSamplingInterval,dblUseMaxDur);
		dblBinWidth = median(diff(vecRefT2));
		dblAnovaP=anova1(matTracePerTrial,[],'off');
		dblAnovaDur = toc(hTicA);
		
		%TS-ZETA new
		%zetatstest
		intPlot = 0;
        dblJitterSize = [];
        [dblZetaP,sZETA] = zetatstest(vecTimestamps,vecdFoF,matEventTimes,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize);

		%normal zeta
		[dblZetaP_sp,sZETA_sp] = zetatest(vecSpikeTimes, matEventTimes,dblUseMaxDur,intResampNum,intPlot);
		
		%pause
		% assign data
		dblMeanP = sZETA.dblMeanP;
		dblMeanZ = -norminv(dblMeanP/2);
		dblZetaZ = sZETA.dblZETA;
		matTsZeta(intNeuron,intRunType) = dblZetaP;
		matZeta_sp(intNeuron,intRunType) = dblZetaP_sp;
		matTtest(intNeuron,intRunType) = dblMeanP;
		matAnova(intNeuron,intRunType) = dblAnovaP;
	end
end

%% save
if boolSave
	save([strDataTargetPath strRec '.mat' ],...
		'matAnova','matZeta_sp','matTtest','matTsZeta','strRec');
end
