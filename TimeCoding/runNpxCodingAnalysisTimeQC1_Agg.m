%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: can single-neuron variability be explained as population-level gain multiplication?
> estimate tuning curve from real data, then apply trial-by-trial gain to all neurons
>model predicting firing rate as combination of mean-rate multiplied by tuning curve => what do
residuals look like?
A: gain axis and mean-rate axis do not align, but ori tuning is distributed as conic manifold around
gain axis. Using stim-specific gain coupling calculated cross-validated for each neuron, allows pop
response to be predicted with R^2 of 0.72

%see articles:
https://elifesciences.org/articles/8998
https://www.nature.com/articles/nn.3711
https://www.jneurosci.org/content/39/37/7344.abstract
etc

==========
20220803:

As a function of mean-rate, adjacent stimulus discriminability could theoretically be anything, such as:
1) Fixed variance, tuning curves that scale with mean rate (= linear increasing curve)
=> perhaps Poisson cells with a wider dynamic range already show this property
2) Smoothly, slowly saturating Poisson neurons (constant for a range of mean rates, then reduction to d?=0)
=> perhaps this can also be made to show a decreasing curve if saturation starts very early
3) Fixed tuning curves, scaling variance (=linear decreasing curve)

To do: show this in a figure

%}
%% define qualifying areas
strRunType = 'ABI'; %ABI or Npx?
runHeaderPopTimeCoding;
if strcmp(strRunType,'ABI')
	runLoadABI;
else
	runLoadNpx;
end

%what to run?
vecRandomize = 1:10; %5=fixed variance, scaling tuning; 6=smoothly saturating poisson
boolSaveData = true;
boolMakeFigs = true;
boolSaveFigs = true;
boolDoDecodingAnalysis = true;
strRunStim = 'DG';%DG or NM

%% pre-allocate matrices
intAreas = numel(cellUseAreas{1});
dblStartT = 0;
intProjType = 2; %1=train+project per quantile; 2=train overall + project per quantile
intQuantiles = 5;

%% go through recordings
tic
for intRec=1:intRecNum
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		runRecPrepABI;
	elseif strcmp(strRunType,'Npx')
		runRecPrepNpx;
	end
	if intNeuronsInArea == 0 || intNeuronNum < 25
		fprintf('Number of neurons is %d for %s: skipping... [%s]\n',intNeuronNum,strRecOrig,getTime);
		continue;
	end
	
	%% run analysis
	close all;
	for intRandomize=vecRandomize
		%simple "rate code"
		matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
		matMeanRate = matSpikeCounts./dblUseMaxDur;
		
		dblLambda = 1;%1
		intTypeCV = 2;
		%change name
		if intRandomize == 1
			strType = 'Real';
			%real
		elseif intRandomize == 2
			strType = 'Shuff';
			%trial-shuffled
		elseif intRandomize == 3
			strType = 'Poiss';
			%poisson process neurons
		elseif intRandomize == 4
			strType = 'UniStretch';
			%shuffle activity like before (repetitions per stimulus randomly and independently permuted for
			%each neuron), but then rescale each trial to the pop mean (or by gain?) in that trial of the
			%original data set; this will recapture the original distribution of population firing rates, while
			%keeping the noise correlation uniform in all directions except the gain
		elseif intRandomize == 5
			strType = 'SdFixed';
			%uniform sd, mean as normal
		elseif intRandomize == 6
			strType = 'Saturating';
			%saturating poisson
		elseif intRandomize == 7
			strType = 'SdScaling';
			%sd scales as mean
		elseif intRandomize == 8
			strType = 'SdLinear';
		elseif intRandomize == 9
			strType = 'SdQuad';
		elseif intRandomize == 10
			strType = 'SdCube';
		end
		strName = [strRunStim '_' strRec '_' strType '_' strAreaAbbr];
		
		%get population gain
		vecGainAx = mean(matMeanRate,2);
		vecPopGainPerTrial=getProjOnLine(matMeanRate,vecGainAx);
		vecPopGainFactor = vecPopGainPerTrial ./ mean(vecPopGainPerTrial);
		if any(vecPopGainFactor==0),continue;end
		fprintf('Running data type %s for %s [%s]\n',strType,strRecOrig,getTime);
		
		%get population mean
		vecOldMean = mean(matMeanRate,1);
		vecPopMeanFactor = vecOldMean ./ mean(vecOldMean);
		vecOldSd = std(matMeanRate,[],1);
		vecPopSdFactor = vecOldSd ./ mean(vecOldSd);
		
		%randomize per orientation
		if intRandomize > 1
			for intN=1:intRespN
				for intStim=1:intStimNum
					vecUseT = find(vecStimIdx==intStim);
					if intRandomize == 2
						%shuffle spikes
						matMeanRate(intN,vecUseT) = matMeanRate(intN,vecUseT(randperm(numel(vecUseT))));
					elseif intRandomize == 3
						%generate spikes
						dblMean = mean(matMeanRate(intN,vecUseT));
						vecRates = poissrnd(dblMean,size(vecUseT));
						matMeanRate(intN,vecUseT) = vecRates;
					elseif intRandomize == 4
						%shuffle spikes & compensate for pop-rate change later
						matMeanRate(intN,vecUseT) = matMeanRate(intN,vecUseT(randperm(numel(vecUseT))));
					else
						%error
					end
				end
				
				%change for each neuron
				if intRandomize == 6
					%smoothly saturating poisson
					vecR = matMeanRate(intN,:);
					
					%logistic slope is 0.5; both k and L increase slope
					vecOldHz = vecR;
					dblSatStart = mean(vecR)/2;
					indSat = vecR > dblSatStart;
					
					fLogistic = @(x,x0,L) x0 + L*2* (-0.5+1./(1+exp(-(2/L)*(x-x0))));
					x = vecR(indSat);
					x0 = dblSatStart;
					L = dblSatStart + 2*sqrt(dblSatStart);
					vecNewSat = fLogistic(x,x0,L);
					
					vecR(indSat) = vecNewSat;
					matMeanRate(intN,:) = vecR;
					
					%scatter(vecOldHz,vecR);
				end
			end
			if intRandomize == 4
				%unistretch
				vecNewMean = mean(matMeanRate,1);
				vecCompensateBy = vecOldMean./vecNewMean;
				matMeanRate = bsxfun(@times,matMeanRate,vecCompensateBy);
			elseif intRandomize == 5
				%fixed sd, scaling tuning
				
				%remove mean
				matMeanRate = bsxfun(@minus,matMeanRate,vecOldMean);
				
				%make sd uniform
				matMeanRate = bsxfun(@rdivide,matMeanRate,vecOldSd);
				
				%add mean back in
				matMeanRate = bsxfun(@plus,matMeanRate,vecOldMean);
			elseif intRandomize == 7
				%scaling sd as pop mean
				
				%remove mean
				matMeanRate = bsxfun(@minus,matMeanRate,vecOldMean);
				
				%make sd uniform
				matMeanRate = bsxfun(@rdivide,matMeanRate,vecOldSd);
				
				%multiply sd
				matMeanRate = bsxfun(@times,matMeanRate,vecOldMean); %basically recapitulates real data
				
				%add mean back in
				matMeanRate = bsxfun(@plus,matMeanRate,vecOldMean);
			elseif intRandomize == 8
				%linear scaling sd
				
				%remove mean
				matMeanRate = bsxfun(@minus,matMeanRate,vecOldMean);
				
				%make sd uniform
				matMeanRate = bsxfun(@rdivide,matMeanRate,vecOldSd);
				
				%multiply sd
				matMeanRate = bsxfun(@times,matMeanRate,vecPopMeanFactor);
				
				%add mean back in
				matMeanRate = bsxfun(@plus,matMeanRate,vecOldMean);
			elseif intRandomize == 9
				%linear scaling variance (sd quadratic)
				
				%remove mean
				matMeanRate = bsxfun(@minus,matMeanRate,vecOldMean);
				
				%make sd uniform
				matMeanRate = bsxfun(@rdivide,matMeanRate,vecOldSd);
				
				%multiply sd
				%matMeanRate = bsxfun(@times,matMeanRate,vecOldMean); %basically recapitulates real data
				matMeanRate = bsxfun(@times,matMeanRate,vecPopMeanFactor.^2);
				
				%add mean back in
				matMeanRate = bsxfun(@plus,matMeanRate,vecOldMean);
			elseif intRandomize == 10
				%sd cubic
				
				%remove mean
				matMeanRate = bsxfun(@minus,matMeanRate,vecOldMean);
				
				%make sd uniform
				matMeanRate = bsxfun(@rdivide,matMeanRate,vecOldSd);
				
				%multiply sd
				%matMeanRate = bsxfun(@times,matMeanRate,vecOldMean); %basically recapitulates real data
				matMeanRate = bsxfun(@times,matMeanRate,vecPopMeanFactor.^3);
				
				%add mean back in
				matMeanRate = bsxfun(@plus,matMeanRate,vecOldMean);
				
			end
		end
		%ensure positive rates
		matMeanRate(matMeanRate<0.1)=0.1;
		
		%% plot population mean + sd over neurons, compare with neuron mean+sd over trials
		intTrials = size(matMeanRate,2);
		
		%real pop mean+sd
		vecPopMean = mean(matMeanRate,1);
		vecPopSd = std(matMeanRate,[],1);
		vecPopCv = vecPopSd./vecPopMean;
		vecPopFano = (vecPopSd.^2)./vecPopMean;
		
		vecFitX = linspace(min(vecPopMean),max(vecPopMean),100);
		[fitobject,gof] = fit(vecPopMean',vecPopSd','poly1');
		matCoefInt = confint(fitobject);
		dblSlope = fitobject.p1;
		vecSlopeCI = matCoefInt(:,1);
		dblIntercept = fitobject.p2;
		vecIntCI = matCoefInt(:,2);
		dblAdjR2 = gof.adjrsquare;
		[vecPredCI,vecPredY] = predint(fitobject,vecFitX);
		
		%real neuron mean+sd
		vecNeuronMean = mean(matMeanRate,2);
		vecNeuronSd = std(matMeanRate,[],2);
		vecNeuronCv = vecNeuronSd./vecNeuronMean;
		vecNeuronFano = (vecNeuronSd.^2)./vecNeuronMean;
		
		vecNeuronFitX = linspace(min(vecNeuronMean),max(vecNeuronMean),100);
		[fitobjectNeuron,gofNeuron] = fit(vecNeuronMean,vecNeuronSd,'poly1');
		matNeuronCoefInt = confint(fitobjectNeuron);
		dblNeuronSlope = fitobjectNeuron.p1;
		vecNeuronSlopeCI = matNeuronCoefInt(:,1);
		dblNeuronIntercept = fitobjectNeuron.p2;
		vecNeuronIntCI = matNeuronCoefInt(:,2);
		dblNeuronAdjR2 = gofNeuron.adjrsquare;
		[vecNeuronPredCI,vecNeuronPredY] = predint(fitobjectNeuron,vecNeuronFitX);
		
		%pop mean+sd, generated according to poisson cells
		matGenRate = poissrnd(repmat(vecNeuronMean,[1 intTrialNum]));
		vecGenPopMean = mean(matGenRate,1);
		vecGenPopSd = std(matGenRate,[],1);
		vecGenPopCv = vecGenPopSd./vecGenPopMean;
		vecGenPopFano = (vecGenPopSd.^2)./vecGenPopMean;
		
		vecGenFitX = linspace(min(vecGenPopMean),max(vecGenPopMean),100);
		[fitobjectGen,gofGen] = fit(vecGenPopMean',vecGenPopSd','poly1');
		matGenCoefInt = confint(fitobjectGen);
		dblGenSlope = fitobjectGen.p1;
		vecGenSlopeCI = matGenCoefInt(:,1);
		dblGenIntercept = fitobjectGen.p2;
		vecGenIntCI = matGenCoefInt(:,2);
		dblGenAdjR2 = gofGen.adjrsquare;
		[vecGenPredCI,vecGenPredY] = predint(fitobjectGen,vecGenFitX);
		
		%neuron mean+sd, generated according to poisson cells
		vecGenNeuronMean = mean(matGenRate,2);
		vecGenNeuronSd = std(matGenRate,[],2);
		vecGenNeuronCv = vecGenNeuronSd./vecGenNeuronMean;
		vecGenNeuronFano = (vecGenNeuronSd.^2)./vecGenNeuronMean;
		
		vecGenNeuronFitX = linspace(min(vecGenNeuronMean),max(vecGenNeuronMean),100);
		[fitobjectGenNeuron,gofGenNeuron] = fit(vecGenNeuronMean,vecGenNeuronSd,'poly1');
		matGenNeuronCoefInt = confint(fitobjectGenNeuron);
		dblGenNeuronSlope = fitobjectGenNeuron.p1;
		vecGenNeuronSlopeCI = matGenNeuronCoefInt(:,1);
		dblGenNeuronIntercept = fitobjectGenNeuron.p2;
		vecGenNeuronIntCI = matGenNeuronCoefInt(:,2);
		dblGenNeuronAdjR2 = gofGenNeuron.adjrsquare;
		[vecGenNeuronPredCI,vecGenNeuronPredY] = predint(fitobjectGenNeuron,vecGenNeuronFitX);
		
		if boolMakeFigs
			figure;maxfig
			subplot(3,4,1)
			hold on
			scatter(vecPopMean,vecPopSd,[],lines(1),'.');
			plot(vecFitX,vecPredY,'k--');
			%plot(vecFitX,vecPredCI(:,1),'--','color',[0.5 0.5 0.5]);
			%plot(vecFitX,vecPredCI(:,2),'--','color',[0.5 0.5 0.5]);;
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblAdjR2,dblSlope,mean(abs(vecSlopeCI-dblSlope))));
			xlabel('Population mean rate (Hz)');
			ylabel('Pop sd (Hz)');fixfig;
			
			subplot(3,4,5)
			scatter(vecPopMean,vecPopCv,'.')
			xlabel('Population mean rate (Hz)');
			ylabel('Pop cv');
			title(strType);
			fixfig;
			
			subplot(3,4,9)
			scatter(vecPopMean,vecPopFano,'.')
			xlabel('Population mean rate (Hz)');
			ylabel('Pop Fano factor');fixfig;
			
			subplot(3,4,2)
			hold on
			scatter(vecNeuronMean,vecNeuronSd,[],lines(1),'.');
			plot(vecNeuronFitX,vecNeuronPredY,'k--');
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblNeuronAdjR2,dblNeuronSlope,mean(abs(vecNeuronSlopeCI-dblNeuronSlope))));
			xlabel('Neuron mean rate (Hz)');
			ylabel('Neuron sd (Hz)');fixfig;
			
			subplot(3,4,6)
			scatter(vecNeuronMean,vecNeuronCv,'.')
			xlabel('Neuron mean rate (Hz)');
			ylabel('Neuron cv');fixfig;
			
			subplot(3,4,10)
			scatter(vecNeuronMean,vecNeuronFano,'.')
			xlabel('Neuron mean rate (Hz)');
			ylabel('Neuron Fano factor');fixfig;
			
			subplot(3,4,3)
			hold on
			scatter(vecGenPopMean,vecGenPopSd,'.')
			plot(vecGenFitX,vecGenPredY,'k--');
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblGenAdjR2,dblGenSlope,mean(abs(vecGenSlopeCI-dblGenSlope))));
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop sd (Hz)');fixfig;
			
			subplot(3,4,7)
			scatter(vecGenPopMean,vecGenPopCv,'.')
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop cv');fixfig;
			
			subplot(3,4,11)
			scatter(vecGenPopMean,vecGenPopFano,'.')
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop Fano factor');fixfig;
			
			subplot(3,4,4)
			hold on
			scatter(vecGenNeuronMean,vecGenNeuronSd,[],lines(1),'.');
			plot(vecGenNeuronFitX,vecGenNeuronPredY,'k--');
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblGenNeuronAdjR2,dblGenNeuronSlope,mean(abs(vecGenNeuronSlopeCI-dblGenNeuronSlope))));
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron sd (Hz)');fixfig;
			
			subplot(3,4,8)
			scatter(vecGenNeuronMean,vecGenNeuronCv,'.')
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron cv');fixfig;
			
			subplot(3,4,12)
			scatter(vecGenNeuronMean,vecGenNeuronFano,'.')
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron Fano factor');fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC1a%s_PopSpikeStatistics.tif',strName)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC1a%s_PopSpikeStatistics.pdf',strName)));
			end
		end
		
		if boolDoDecodingAnalysis
			%% run decoding
			%run samples with balanced trials
			intIters = 100;
			matSplitPerf = nan(intQuantiles,intIters);
			for intIter=1:intIters
				%select balanced trials
				vecPopRate = sum(matMeanRate,1);
				vecUseBalancedTrials = nan(1,intSplitTrialsPerOri*intStimNum*intQuantiles);
				vecBalancedQ = nan(1,intSplitTrialsPerOri*intStimNum*intQuantiles);
				intCounterBT = 1;
				vecPriorDistributionSplit = ones(size(vecPriorDistribution))*intSplitTrialsPerOri;
				vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
				vecStartTrials(end)=[];
				matUseTrials = nan(intStimNum,intSplitTrialsPerOri,intQuantiles);
				for intQ=1:intQuantiles
					vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
					for intOri=1:intStimNum
						vecThisOri = find(cellSelect{intOri});
						[vecSorted,vecReorder]=sort(vecPopRate(vecThisOri));
						vecQualifyingTrials = vecThisOri(vecReorder(vecUseTrialsTemp));
						matUseTrials(intOri,:,intQ) = vecQualifyingTrials;
						
						vecUseBalancedTrials(intCounterBT:(intCounterBT+intSplitTrialsPerOri-1)) = vecQualifyingTrials(randperm(intSplitTrialsPerOri,intSplitTrialsPerOri));
						vecBalancedQ(intCounterBT:(intCounterBT+intSplitTrialsPerOri-1)) = intQ;
						intCounterBT = intCounterBT + intSplitTrialsPerOri;
					end
				end
				
				%perform overall decoding
				[vecBalancedTrialTypeIdx,vecBalancedUnique,vecBalancedPriorDistribution,cellBalancedSelect,vecBalancedRepetition] = val2idx(vecStimIdx(vecUseBalancedTrials));
				[dblPerf,vecDecodedIndexCV,matPosteriorProbability] = ...
					doCrossValidatedDecodingLR(matMeanRate(:,vecUseBalancedTrials),vecStimIdx(vecUseBalancedTrials),intTypeCV,[],dblLambda);
				vecConfidence = nan(size(vecDecodedIndexCV));
				for intTrial=1:numel(vecBalancedTrialTypeIdx)
					vecConfidence(intTrial) = matPosteriorProbability(vecBalancedTrialTypeIdx(intTrial),intTrial);
				end
				
				%assign per quantile
				matCorrPerQ = nan(intTrialsPerQ,intQuantiles);
				vecCorr = vecDecodedIndexCV(:) == vecBalancedTrialTypeIdx(:);
				for intQ=1:intQuantiles
					matCorrPerQ(:,intQ) = vecCorr(vecBalancedQ==intQ);
				end
				
				%save data
				matSplitPerf(:,intIter) = sum(matCorrPerQ,1)/intTrialsPerQ;
			end
			
			%calc mean and sem
			vecSplitPerfMu = mean(matSplitPerf,2)';
			vecSplitPerfSem = std(matSplitPerf,[],2)';
			
			%perform decoding per quantile
			vecQuantilePerf = nan(1,intQuantiles);
			matQuantileConfusion = nan(intStimNum,intStimNum,intQuantiles);
			%split trials into quantiles
			figure;maxfig;
			vecPlotQ = [];%[1 ceil(intQuantiles/2) intQuantiles];
			vecTrialQuantile = zeros(1,intTrials);
			for intQ=1:intQuantiles
				%divide into quantiles
				vecUseTrials = sort(flat(matUseTrials(:,:,intQ)));
				vecTrialQuantile(vecUseTrials) = intQ;
				
				matUseRate = matMeanRate(:,vecUseTrials);
				vecOriUse = vecStimIdx(vecUseTrials);
				
				%[dblPerfQ,vecDecodedIndexRateCV,matPosteriorProbabilityRate,dblMeanErrorDegsRate,matConfusion,matWeightsRate] = ...
				%	doCrossValidatedDecodingLR(matUseRate,vecOriUse,intTypeCV,vecPriorDistributionSplit,dblLambda);
				[dblPerfQ,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation] = ...
					doCrossValidatedDecodingLR(matUseRate,vecOriUse,intTypeCV,[],dblLambda);
				vecQuantilePerf(intQ) = dblPerfQ;
				matQuantileConfusion(:,:,intQ) = matConfusion;
				intPlot = find(vecPlotQ==intQ);
				if ~isempty(intPlot)
					subplot(2,3,intPlot)
					imagesc(matConfusion)
					xlabel('Real Orientation')
					ylabel('Decoded Orientation')
					title(sprintf('Quantile %d',intQ))
					fixfig;grid off;
				end
			end
			
			% overall train, output split
			clf;
			subplot(2,3,4)
			hold on
			plot([1 numel(vecQuantilePerf)],(1/intStimNum)*[1 1],'--','color',[0.5 0.5 0.5]);
			errorbar(1:intQuantiles,vecSplitPerfMu,vecSplitPerfSem,'color',lines(1));
			xlabel('Pop. activity quantile');
			ylabel('CV decoding accuracy');
			title('Train once on all, test per quantile');
			fixfig;
			
			%split train, output split
			dblAlphaEquivOfSd = normcdf(1)-normcdf(-1);
			[phat2,pci2] = binofit(vecQuantilePerf*intTrialsPerQ,intTrialsPerQ,dblAlphaEquivOfSd);
			subplot(2,3,5)
			hold on
			plot([1 numel(vecQuantilePerf)],(1/intStimNum)*[1 1],'--','color',[0.5 0.5 0.5]);
			errorbar(1:intQuantiles,phat2,phat2'-pci2(:,1),phat2'-pci2(:,2),'color',lines(1));
			xlabel('Pop. activity quantile');
			ylabel('CV decoding accuracy');
			title('Train+test per quantile');
			fixfig;
			normaxes;
			
			%difference
			%[phat3,pci3] = binofit((phat-vecQuantilePerf)*intTrialsPerQ,intTrialsPerQ,dblAlphaEquivOfSd);
			
			subplot(2,3,6)
			hold on
			%plot([1 numel(vecQuantilePerf)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
			plot(1:intQuantiles,(vecSplitPerfMu-vecQuantilePerf));
			xlabel('Pop. activity quantile');
			ylabel('Generalization penalty (\Deltaaccuracy)');
			title('Accuracy difference');
			fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC1b%s_QuantileDecoding.tif',strName)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC1b%s_QuantileDecoding.pdf',strName)));
			end
		else
			%calc mean and sem
			vecSplitPerfMu = [];
			vecSplitPerfSem = [];
			
			%perform decoding per quantile
			matSplitPerf = [];
			vecQuantilePerf = [];
			matQuantileConfusion = [];
		end
		
		%% prep analysis
		%select balanced trials
		dblLambda = 1;
		vecPopRate = sum(matMeanRate,1);
		vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
		vecStartTrials(end)=[];
		matUseTrials = nan(intStimNum,intSplitTrialsPerOri,intQuantiles);
		vecTrialQuantile = zeros(1,intTrialNum);
		for intQ=1:intQuantiles
			vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
			for intStim=1:intStimNum
				vecThisStim = find(cellSelect{intStim});
				[vecSorted,vecReorder]=sort(vecPopRate(vecThisStim));
				vecQualifyingTrials = vecThisStim(vecReorder(vecUseTrialsTemp));
				matUseTrials(intStim,:,intQ) = vecQualifyingTrials;
			end
			%divide into quantiles
			vecUseTrials = sort(flat(matUseTrials(:,:,intQ)));
			vecTrialQuantile(vecUseTrials) = intQ;
		end
		
		%% example
		%for intUseStim=1:intStimNum
		intUseStim = 1;
		vecStim1 = find(vecStimIdx==intUseStim);
		vecStim2 = find(vecStimIdx==modx(round(intUseStim+intStimNum/2),intStimNum));
		vecUseTrials = cat(2,vecStim1,vecStim2);
		vecTQR = vecTrialQuantile(vecUseTrials);
		vecUseTrials(vecTQR==0)=[];
		vecTQR = vecTrialQuantile(vecUseTrials);
		vecStimR = vecStimIdx(vecUseTrials);
		vecPMR = vecPopMean(vecUseTrials);
		matStimPair = matMeanRate(:,vecUseTrials)';
		indStim1=vecStimR==min(vecStimR);
		indStim2=vecStimR==max(vecStimR);
		
		%% plot
		if boolMakeFigs && intUseStim == 1
			
			matColMap = redbluepurple(intQuantiles);
			matColor = matColMap(vecTQR,:);
			
			figure;maxfig;
			%make plots of distributions of activation per quantile
			subplot(2,3,1)
			title(sprintf('%s, %s %s %s', strRec,strArea,strRunStim,strType),'interpreter','none');
			
			subplot(2,3,3)
			hold on;
		end
		dblStep = 1;
		vecBinE = (-10:dblStep:10)/10;
		vecBinC = vecBinE(2:end)-dblStep/2;
		vecAllAct = nan(1,numel(vecTQR));
		vecAbsW = nan(1,intQuantiles);
		vecBinaryPerf = nan(1,intQuantiles);
		beta0 = [0;0];
		
		if intProjType == 2
			matUseResp = matStimPair';
			vecUseOri = val2idx(vecStimR);
			[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation,matAggWeights,vecAggRep] = ...
				doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],1);%dblLambda);
			%calculate normalization factors
			vecBinIdx = val2idx(vecUseOri);
			vecRepWeights = squeeze(sum(abs(matAggWeights(:,1,:)),1));
			vecNormFactors = vecRepWeights(vecAggRep);
		end
		
		for intQ=1:intQuantiles
			vecUseTrialsQ = vecTQR==intQ;
			vecUseOriQ = val2idx(vecStimR(vecUseTrialsQ));
			if intProjType == 1
				%% proj 1
				matUseResp = matStimPair(vecUseTrialsQ,:)';
				[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
					doCrossValidatedDecodingLR(matUseResp,vecUseOriQ,intTypeCV,[],1);%dblLambda);
				
				vecBinaryPerf(intQ) = dblPerfP;
				vecAbsW(intQ) = sum(abs(matWeightsBin(:,1)));
				vecAllAct(vecUseTrialsQ) = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
				
				%split by group & plot
				vecAct = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
			else
				%% proj 2
				%split by group & plot
				
				vecBinaryPerf(intQ) = sum(vecDecodedIndexRateCV(vecUseTrialsQ) == vecBinIdx(vecUseTrialsQ)) / sum(vecUseTrialsQ);
				vecAbsW(intQ) = mean(vecNormFactors(vecUseTrialsQ));
				vecAllAct(vecUseTrialsQ) = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
				
				%split by group & plot
				vecAct = matActivation(1,vecUseTrialsQ)'./vecNormFactors(vecUseTrialsQ);
			end
			%% plot
			vecCounts1 = histcounts(vecAct(vecUseOriQ==1),vecBinE);
			vecCounts2 = histcounts(vecAct(vecUseOriQ==2),vecBinE);
			if boolMakeFigs && intUseStim == 1
				plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
				plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
			end
		end
		
		if boolMakeFigs && intUseStim == 1
			hold off;
			set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
			ylabel('Quantile; y=trials per quantile');
			xlabel('LR activation');
			title('0 vs 90 degrees');
			fixfig;grid off
			
			%plot mean & decision axis
			subplot(2,3,2)
			hold on;
			scatter(vecPMR(indStim1),vecAllAct(indStim1),[],matColor(indStim1,:),'d');
			scatter(vecPMR(indStim2),vecAllAct(indStim2),[],matColor(indStim2,:),'*');
			hold off
			xlabel('Mean pop activity');
			ylabel('LR activation');
			title('point=trial, color=pop activity quantile');
			fixfig;
		end
		
		%% average over all orthogonal (or adjacent?) stimuli
		matBinaryPerf = nan(intQuantiles,intStimNum);
		cellLRActPerQ = cell(intQuantiles,intStimNum,2);
		cellPopMuPerQ = cell(intQuantiles,intStimNum,2);
		
		[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimIdx);
		
		%% get projection onto activation axis, type 1 or 2
		if intProjType == 1
			%% get axes and projection per quantile
			for intQ=1:intQuantiles
				for intOriIdx1 = 1:intStimNum
					intOriIdx2 = intOriIdx1+1;
					%intOriIdx2 = intOriIdx1+floor(intStimNum/2)-1;
					intOriIdx2 = modx(intOriIdx2,intStimNum);
					
					vecUseTrials = vecTrialQuantile==intQ & (vecStimIdx == intOriIdx1 | vecStimIdx == intOriIdx2);
					
					matUseResp = matMeanRate(:,vecUseTrials)';
					vecUseOri = vecStimIdx(vecUseTrials);
					[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
						doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],dblLambda);
					
					
					%split by group & plot
					matBinaryPerf(intQ,intOriIdx1) = dblPerfP;
					vecAct = matActivation(1,:)/sum(abs(matWeightsBin(:,1)));
					vecAct1 = vecAct(vecUseOri==intOriIdx1);
					vecAct2 = vecAct(vecUseOri==intOriIdx2);
					if mean(vecAct1) > mean(vecAct2)
						[vecAct2,vecAct1]= swap(vecAct1,vecAct2);
					end
					cellLRActPerQ{intQ,intOriIdx1,1} = vecAct1;
					cellLRActPerQ{intQ,intOriIdx1,2} = vecAct2;
					cellPopMuPerQ{intQ,intOriIdx1,1} = mean(matUseResp(vecUseOri==intOriIdx1,:),2);
					cellPopMuPerQ{intQ,intOriIdx2,2} = mean(matUseResp(vecUseOri==intOriIdx2,:),2);
				end
			end
		else
			%% get axes from all data, then project per quantile
			for intOriIdx1 = 1:intStimNum
				intOriIdx2 = intOriIdx1+1;
				%intOriIdx2 = intOriIdx1+floor(intStimNum/2)-1;
				intOriIdx2 = modx(intOriIdx2,intStimNum);
				
				indUseTrainTrials = (vecStimIdx == intOriIdx1 | vecStimIdx == intOriIdx2);
				matUseResp = matMeanRate(:,indUseTrainTrials)';
				vecUseOri = vecStimIdx(indUseTrainTrials);
				[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation,matAggWeights,vecAggRep] = ...
					doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],dblLambda);
				%calculate normalization factors
				vecBinIdx = val2idx(vecStimIdx(indUseTrainTrials));
				vecQuantInOri = vecTrialQuantile(indUseTrainTrials);
				vecRepWeights = squeeze(sum(abs(matAggWeights(:,1,:)),1));
				vecNormFactors = vecRepWeights(vecAggRep);
				
				for intQ=1:intQuantiles
					indUseTrials = vecQuantInOri==intQ;
					%split by group & plot
					matBinaryPerf(intQ,intOriIdx1) = sum(vecDecodedIndexCV(indUseTrials) == vecBinIdx(indUseTrials)) / sum(indUseTrials);
					vecAct = matAggActivation(1,indUseTrials)'./vecNormFactors(indUseTrials);
					vecAct1 = vecAct(vecUseOri(indUseTrials)==intOriIdx1);
					vecAct2 = vecAct(vecUseOri(indUseTrials)==intOriIdx2);
					if mean(vecAct1) > mean(vecAct2)
						[vecAct2,vecAct1]= swap(vecAct1,vecAct2);
					end
					cellLRActPerQ{intQ,intOriIdx1,1} = vecAct1;
					cellLRActPerQ{intQ,intOriIdx1,2} = vecAct2;
					cellPopMuPerQ{intQ,intOriIdx1,1} = mean(matUseResp(vecUseOri==intOriIdx1 & indUseTrials,:),2);
					cellPopMuPerQ{intQ,intOriIdx1,2} = mean(matUseResp(vecUseOri==intOriIdx2 & indUseTrials,:),2);
				end
			end
		end
		
		if boolMakeFigs
			dblStep = 1;
			vecBinE = -10:dblStep:10;
			vecBinC = vecBinE(2:end)-dblStep/2;
			
			subplot(2,3,4);
			hold on
			for intQ=1:intQuantiles
				%plot distros
				vecAct1 = cell2vec(cellLRActPerQ(intQ,:,1));
				vecAct2 = cell2vec(cellLRActPerQ(intQ,:,2));
				vecCounts1 = histcounts(vecAct1,vecBinE);
				vecCounts2 = histcounts(vecAct2,vecBinE);
				if boolMakeFigs
					plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
					plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
				end
			end
			%finish plot
			hold off;
			set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
			ylabel('Quantile; y=trials per quantile');
			xlabel('LR activation');
			title('Mean over ~orth stim pairs');
			fixfig;grid off
			
			%plot d', variance and distance in mean
			matDprime = nan(intQuantiles,intStimNum);
			matPooledSd = nan(intQuantiles,intStimNum);
			matMeanD = nan(intQuantiles,intStimNum);
			matQ = nan(intQuantiles,intStimNum);
			for intQ=1:intQuantiles
				for intOriIdx = 1:intStimNum
					matDprime(intQ,intOriIdx) = abs(getdprime2(cellLRActPerQ{intQ,intOriIdx,1},cellLRActPerQ{intQ,intOriIdx,2}));
					matPooledSd(intQ,intOriIdx) = (std(cellLRActPerQ{intQ,intOriIdx,1}) + std(cellLRActPerQ{intQ,intOriIdx,2}))/2;
					matMeanD(intQ,intOriIdx)  = abs(mean(cellLRActPerQ{intQ,intOriIdx,1}) - mean(cellLRActPerQ{intQ,intOriIdx,2}));
					matQ(intQ,intOriIdx) = intQ;
				end
			end
			
			matColMap = redbluepurple(intQuantiles);
			matColor2 = matColMap(matQ(:),:);
			%{
				subplot(2,3,4)
				scatter(matDprime(:),matBinaryPerf(:),[],matColor2)
				xlabel('d''');
				ylabel('Ori-pair decoding accuracy');
				title('Point = stim+quantile combination');
				fixfig;
			%}
			%{
				subplot(2,3,5)
				scatter(matPooledSd(:),matBinaryPerf(:),[],matColor2)
				xlabel('Sd over trials');
				ylabel('Ori-pair decoding accuracy');
				title('Point = stim+quantile combination');
				fixfig;
			%}
			
			h=subplot(2,3,5);
			colormap(h,matColMap);
			%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
			%calc mean+sem per q
			vecMeanDprime = mean(matDprime,2);
			vecSemDprime = std(matDprime,[],2)./sqrt(intStimNum);
			vecMeanSd = mean(matPooledSd,2);
			vecSemSd = std(matPooledSd,[],2)./sqrt(intStimNum);
			vecMeanMu = mean(matMeanD,2);
			vecSemMu = std(matMeanD,[],2)./sqrt(intStimNum);
			hold on
			cline(h,vecMeanSd,vecMeanMu,[],1:5);
			for intQ=1:intQuantiles
				errorbar(vecMeanSd(intQ),vecMeanMu(intQ),vecSemSd(intQ)/2,vecSemSd(intQ)/2,vecSemMu(intQ)/2,vecSemMu(intQ)/2,'x','color',matColMap(intQ,:));
			end
			hold off
			xlabel('Sd over trials');
			ylabel('Mean over trials');
			title('Point = stim+quantile mu+/-sem');
			fixfig;
			
			h=subplot(2,3,6);
			colormap(h,matColMap);
			hold on
			cline(h,vecMeanSd,vecMeanDprime,[],1:5);
			for intQ=1:intQuantiles
				errorbar(vecMeanSd(intQ),vecMeanDprime(intQ),vecSemSd(intQ)/2,vecSemSd(intQ)/2,vecSemDprime(intQ)/2,vecSemDprime(intQ)/2,'x','color',matColMap(intQ,:));
			end
			hold off
			xlabel('Sd over trials');
			ylabel('d''');
			title('Point = stim+quantile mu+/-sem');
			fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC1c%s_DynamicCoding.tif',strName)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC1c%s_DynamicCoding.pdf',strName)));
			end
		end
		
		%% save data
		if boolSaveData
			save([strTargetDataPath 'QC1Data' strName '.mat'],...
				'strRec','strArea','strRunStim','strType',...
				'vecStimIdx','vecTrialQuantile',...
				'matMeanRate','matBinaryPerf','vecQuantilePerf','matSplitPerf','cellPopMuPerQ','cellLRActPerQ',...
				'dblBC','dblMaxDevFrac');
		end
	end
end
toc
