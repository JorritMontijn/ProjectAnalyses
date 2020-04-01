%% load data
clear all;
strRec = 'Exp2019-11-20';
strExp = [strRec '_MP2_AP.mat'];
strDataPath = ['P:\Montijn\DataPreProcessed\' strRec];
strFigOut = ['D:\Data\Results\ClusterFigures\' strRec];
mkdir(strFigOut);
load(fullfile(strDataPath,strExp));
sCluster = sAP.sCluster;
cellStim = sAP.cellStim;
intClusters = numel(sCluster);

%% select stimulus block
intUseStim = 1;
structEP = cellStim{intUseStim}.structEP;

%% retrieve properties
vecDepth = cell2vec({sCluster(:).Depth});
vecSpikesPerCell = flat(cellfun(@numel,{sCluster(:).SpikeTimes}));

%% prepare variables
%get types
sTypes = getStimulusTypes(structEP,{'Orientation'});
cellSelect = getSelectionVectors(structEP,sTypes);
intTrials = numel(structEP.vecStimOnTime);

%get stim dur
vecStimDur = structEP.vecStimOffTime - structEP.vecStimOnTime;
dblRealStimDur = min(vecStimDur);
dblUseStimDur = 1;
matStimOnOff = [structEP.vecStimOnTime;structEP.vecStimOnTime+dblUseStimDur]';

%get base dur
vecBaseDur =  structEP.vecStimOnTime(2:end) - structEP.vecStimOffTime(1:(end-1));
dblRealBaseDur = min(vecBaseDur);
dblUseBaseDur = 0.3;
matBaseOnOff = [structEP.vecStimOnTime-dblUseBaseDur;structEP.vecStimOnTime]';
	
%pre-allocate spike matrices
matStimHz = nan(intClusters,intTrials);
matBaseHz = nan(intClusters,intTrials);

% pre-allocate latencies
intPlot = 0;
intResamp = 50;
vecOnsetLatencies = nan(intClusters,1); %onset peak
vecPeakWidths = nan(intClusters,1); %onset width

%pre-allocate responsiveness
vecZetaP = nan(intClusters,1);
vecMeanP = nan(intClusters,1);

%% loop through cells
hTic = tic;
for intCluster=1:intClusters
	%msg
	if toc(hTic) > 5
	fprintf('Processing cluster %d/%d [%s]\n',intCluster,intClusters,getTime);
	hTic = tic;
	end
	%get spikes
	vecCellSpikes = sCluster(intCluster).SpikeTimes;
	
	%assign stim
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecCellSpikes,structEP.vecStimOnTime,dblUseStimDur);
	vecSpikeCounts = accumarray(vecTrialPerSpike,1,[intTrials 1]);
	matStimHz(intCluster,:) = vecSpikeCounts ./ dblUseStimDur;
	
	%assign base
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecCellSpikes,structEP.vecStimOnTime-dblUseBaseDur,dblUseBaseDur);
	vecSpikeCounts = accumarray(vecTrialPerSpike,1,[intTrials 1]);
	matBaseHz(intCluster,:) = vecSpikeCounts ./ dblUseBaseDur;
	
	%get zeta
	[dblZETA,vecLatencies,sZETA,sMSD] = getZeta(vecCellSpikes,matStimOnOff,1.5,intResamp,intPlot,3);
	if isempty(sZETA),continue;end
	vecZetaP(intCluster) = sZETA.dblP;
	vecMeanP(intCluster) = sZETA.dblMeanP;
	%get onset peak
	[dblPeakValue,dblPeakTime,dblPeakWidth] = getPeak(sMSD.vecMSD,sMSD.vecT,[0 dblUseStimDur/2]);
	vecOnsetLatencies(intCluster,:) = dblPeakTime;
	vecPeakWidths(intCluster,:) = dblPeakWidth;
	if intPlot > 0
	%show non-stat, contam
	vecPtrs=get(gcf,'Children');
	title(vecPtrs(1),sprintf('%s',strExp),'Interpreter','none');
	title(vecPtrs(6),sprintf('%s, Clust %d',structEP.strFile,intCluster),'Interpreter','none');
	title(vecPtrs(5),sprintf('Depth=%d, Non-st=%.3f, Cont=%.3f',round(sCluster(intCluster).Depth),sCluster(intCluster).NonStationarity,sCluster(intCluster).Contamination),'Interpreter','none');
	title(vecPtrs(2),sprintf('Peak=%dms, Width=%dms',round(dblPeakTime*1000),round(dblPeakWidth*1000)),'Interpreter','none');
	axes(vecPtrs(2));
	hold on;
	scatter(dblPeakTime,dblPeakValue,'rx');
	hold off;
	
	%save figure
	export_fig(fullfile(strFigOut,sprintf('%sS%dClust%dResp.tif',strRec,intUseStim,intCluster)));
	export_fig(fullfile(strFigOut,sprintf('%sS%dClust%dResp.pdf',strRec,intUseStim,intCluster)));
	close;
	end
end

%get tuning curve properties
vecStimOriDegrees = structEP.Orientation;
[sOut] = getTuningCurves(matStimHz,vecStimOriDegrees);

matFittedParams = sOut.matFittedParams;

%transform to [N x S x R]
[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matStimHz,mod(vecStimOriDegrees,180));
[matBaseNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matBaseHz,mod(vecStimOriDegrees,180));
[intN,intS,intR] = size(matRespNSR);

matR = squeeze(matRespNSR(intNeuron,:,:))';
matB = squeeze(matBaseNSR(intNeuron,:,:))';
[pval,T2] = hotellingt2(matR)
matX = repmat(vecUniqueDegs,[intR 1]);
%scatter(flat(matX+5*rand(size(matX))),flat(matR))
hold on
errorbar(vecUniqueDegs,mean(matR),std(matR)/sqrt(40))
hold off

%	- params(1): prefDir
	%	- params(2): kappa
	%	- params(3): direction index
	%	- params(4): baseline
	%	- params(5): gain
	scatter(vecDepth,rad2deg(matFittedParams(:,1)))
	
	
	 [vecR2,vecP] = doCrossValidatedTuningCurves(matStimHz,mod(vecStimOriDegrees,180),2)
return
	%% get tuning curves & parameters
		dblRho = getTuningRho(vecSpikesPerTrial',vecOrientationsRad');
		dblDeltaPrime = getDeltaPrime(vecSpikesPerTrial',vecOrientationsRad');
		dblOPI = getOPI(vecSpikesPerTrial',vecOrientationsRad');
		sTuning = getTuningCurves(vecSpikesPerTrial',vecOrientationsDeg');
		%[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(vecSpikesPerTrial',vecOrientationsRad');
		
	
		%% save data
		cellArea{intNeuron} = strArea;
		vecParams = sTuning.matFittedParams;
		matFitParams(:,intNeuron) = vecParams;
		vecPrefDir(intNeuron) = vecParams(1);
		matFitResp(:,intNeuron) = feval(sTuning.funcFit,vecParams,vecOriTypesRad);
		%vecParams(1) = pi/2;
		%zeta
		vecZeta(intNeuron) = dblZeta;
		vecHzP(intNeuron) = sOptionalOutputs.dblHzP;
		cellSpikeT{intNeuron} = sOptionalOutputs.vecSpikeT;
		cellZeta{intNeuron} = sOptionalOutputs.vecZ;
		
		%rho
		vecRho(intNeuron) = dblRho;
		%d'
		vecDeltaPrime(intNeuron) = dblDeltaPrime;
		vecOPI(intNeuron) = dblOPI;
		%raw ori, mean per ori, and sd per ori
		vecOriTtest(intNeuron) = sTuning.vecOriTtest;
		cellRawOri{intNeuron} = sTuning.vecUniqueRads;
		cellRawMean{intNeuron} = sTuning.matMeanResp;
		cellRawSD{intNeuron} = sTuning.matSDResp;
		%peak CV + BW
		matVariance(:,intNeuron) = sTuning.matVariance;
		matBandwidth(:,intNeuron) = real(sTuning.matBandwidth);
		
%% select cells
vecNonStat = flat(abs(cell2vec({sCluster(:).NonStationarity})));
vecContam = flat(cell2vec({sCluster(:).Contamination}));

indNonStat = vecNonStat < 0.3;
indZ = vecZetaP < 0.025;
indContam = vecContam < 0.2; %cut off 0.1
indKeepCells = sum(indNonStat & indZ & indContam);

%remove cells with few spikes 
%(e.g., 129)

%remove cells with reduction of firing rate