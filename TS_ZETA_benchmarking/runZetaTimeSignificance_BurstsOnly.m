%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
vecTestResamps = 1:100;
intResampNum = 250;
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = sprintf('SimBurstsOnlyResamp%d',intResampNum);
intNeurons = 100;

%% pre-allocate output variables
matTtest = nan(intNeurons,2);
matZeta = nan(intNeurons,2);
matAnova = nan(intNeurons,2);
matZetaTime = nan(intNeurons,2,numel(vecTestResamps));
matCompTimeZeta = nan(intNeurons,2);
matCompTimeZetaTime = nan(intNeurons,2);
matCompTimeAnova = nan(intNeurons,2);

%% generate data
hTicN = tic;
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
	end
	
	%% spiking params
	dblBaseRate = exprnd(1)+0.1;
	dblBurstEventRate = exprnd(1);
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
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	dblPrefBurstEventRate = 1+dblBurstEventRate; %mean evoked rate of burst events (Hz) (exponential inter-event times)
	
	sTuningParams.boolDoublePeaked = boolDoublePeaked;
	sTuningParams.dblPrefOri = dblPrefOri;
	sTuningParams.dblKappa = dblKappa;
	sTuningParams.dblPrefRate = dblPrefRate;
	sTuningParams.dblPrefBurstEventRate = dblPrefBurstEventRate;
	
	%% stimulus data
	dblStimDur = 1;
	dblPreBaseDur = 0.25;
	dblPostBaseDur = 0.25;
	dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
	intOris = 12;
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
	vecSpikeTimes = vecSpikeTimes1;
	
	%% test
	%real+rand
	for intRunType=vecRunTypes
		
		%randomize
		if intRunType ==2
			vecJitter = 2*dblTrialDur*((rand(size(vecTrialStart))-1/2)*2);
			vecUseTrialStart = vecTrialStart + vecJitter;
			strRand = 'Rand';
		else
			vecUseTrialStart = vecTrialStart;
			strRand = 'Real';
		end
		matUseTrialT = cat(2,vecUseTrialStart',vecUseTrialStart'+dblStimDur);
		
		intPlot = 0;
		hTicZ = tic;
		dblUseTrialDur = dblTrialDur;
		[dblZetaP,sZETA] = zetatest(vecSpikeTimes,matUseTrialT,dblUseTrialDur,intResampNum,0);
		dblTimeZ = toc(hTicZ);
		dblMeanP = sZETA.dblMeanP;
		
		%% anova
		hTicA=tic;
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,matUseTrialT(:,1),dblUseTrialDur);
		if numel(vecTimePerSpike) < 3,continue;end
		[optN, dblC, allN, allC] = opthist(vecTimePerSpike);
		if optN==1,optN=2;end %at least 2 bins
		dblBinWidth = dblUseTrialDur/optN;
		vecBins = 0:dblBinWidth:dblUseTrialDur;
		matPSTH = nan(intTrials,numel(vecBins)-1);
		for intTrial=1:intTrials
			matPSTH(intTrial,:) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
		end
		dblAnovaP=anova1(matPSTH,[],'off');
		dblAnovaDur = toc(hTicA);
		
		%% calculate IFR of random jitters to set threshold for significant epochs
		%get average of multi-scale derivatives, and rescaled to instantaneous spiking rate
		hTicTZ = tic;
		intRandIters = numel(sZETA.cellRandT);
		matRandRate = nan(numel(sZETA.vecSpikeT),intRandIters);
		vecT = sZETA.vecSpikeT;
		%figure
		%hold on
		for intRandIter=1:intRandIters
			vecRandT = sZETA.cellRandT{intRandIter};
			vecRandD = sZETA.cellRandDiff{intRandIter};
			
			
			dblMeanRandRate = (numel(vecRandT)/(dblUseTrialDur*size(matUseTrialT,1)));
			[vecRandRate,sRate] = getMultiScaleDeriv(vecRandT,vecRandD,[],[],[],0,dblMeanRandRate,dblUseTrialDur);
			vecInterpR = interp1(vecRandT,vecRandRate,vecT);
			matRandRate(:,intRandIter) = vecInterpR;
			%plot(vecT,vecInterpR,'color',[0.5 0.5 0.5]);
		end
		%get real rate
		dblMeanRealRate = (numel(vecT)/(dblUseTrialDur*size(matUseTrialT,1)));
		[vecRealRate,sRate] = getMultiScaleDeriv(vecT,sZETA.vecD,[],[],[],0,dblMeanRealRate,dblUseTrialDur);
		dblTimeTZ = toc(hTicTZ);
		
		%% test multiple resamp nums
		for intResampIdx=1:numel(vecTestResamps)
			intUseResamps = vecTestResamps(intResampIdx);
			vecRandRate = flat(matRandRate(:,randperm(intResampNum,intUseResamps)));
			
			[vecQ,dblZETA] = getZetaP(vecRealRate,vecRandRate,false);
			vecP = 1-abs(1-(vecQ)*2);
			vecP_corr = (1 - (1 - vecP).^(numel(vecQ))); %sidak
			
			%new
			vecP_corr3 = nan(size(vecQ));
			indP_down = vecQ<0.5;
			vecP_corr3(indP_down) = (1 - (1 - vecQ(indP_down)).^(numel(vecQ)));
			%vecP_corr3(indP_down) = (1 - (1 - 2*vecQ(indP_down)).^(numel(vecRandRate)));
			indP_up = vecQ>=0.5;
			vecP_corr3(indP_up) = 1 - ((1 - (1 - vecQ(indP_up))).^numel(vecQ));
			%vecP_corr3(indP_up) = 1 - ((1 - (1 - vecQ(indP_up))*2).^numel(vecRandRate));
			
			%extra
			intType = 2;
			if intType==2
				vecZ = -norminv(vecP_corr); %double the statistic for two-sided test
				%vecZ(vecZ<0)=0;
				vecP_corr2 = (1-normcdf(vecZ))*2; %transform back to p-value
			end
			vecZ = -norminv(vecP_corr/2); %double the statistic for two-sided test
			vecZ2 = -norminv(vecP_corr2/2); %double the statistic for two-sided test
			vecZ3 = -norminv(vecP_corr3/2); %double the statistic for two-sided test
				
			%save
			dblZetaTimeP = min(vecP_corr);
			dblZetaTimeP2 = min(vecP_corr2);
			dblZetaTimeP3 = min(vecP_corr3);
			
			%save
			matZetaTime(intNeuron,intRunType,intResampIdx) = dblZetaTimeP3;
		end
		vecZ_corr = -norminv(vecP_corr/2);
		%real spiking
		if 0
			%% plot
			figure
			hold on
			for intRandIter=1:2;%intRandIters
				plot(vecT,matRandRate,'color',[0.5 0.5 0.5]);
			end
			hold off
			colormap(gca,parula);
			set(gca,'clim',[min(vecZ_corr),max(vecZ_corr)]);
			hold on
			vecH = cline(gca,sZETA.vecSpikeT,vecRealRate,[],vecZ_corr,true);
			indSign = vecP_corr<0.05;
			scatter(sZETA.vecSpikeT(indSign),vecRealRate(indSign),'k.');
			hold off
			h=colorbar;
			clabel(h,'Z-score');
			xlabel('Time (s)')
			ylabel('Spiking rate (Hz)')
			title(sprintf('ZETA+, p=%.3f',dblZetaTimeP3))
			fixfig;
			return
		end
		
		%% save data
		matTtest(intNeuron,intRunType) = dblMeanP;
		matZeta(intNeuron,intRunType) = dblZetaP;
		matAnova(intNeuron,intRunType) = dblAnovaP;
		matCompTimeZeta(intNeuron,intRunType) = dblTimeZ;
		matCompTimeZetaTime(intNeuron,intRunType) = dblTimeTZ;
		matCompTimeAnova(intNeuron,intRunType) = dblAnovaDur;
		
		
		%{
	%% check if match
	dblRandMu = mean(matRandRate(:));
	dblRandVar = var(matRandRate(:));
	
	%derive beta parameter from variance
	dblBeta = (sqrt(6).*sqrt(dblRandVar))./(pi);
	
	%derive mode from mean, beta and E-M constant
	dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
	dblMode = dblRandMu - dblBeta.*dblEulerMascheroni;
	
	%define Gumbel cdf
	fGumbelCDF = @(x) exp(-exp(-((x(:)-dblMode)./dblBeta)));
	fGumbelPDF = @(x) (1./dblBeta)*exp(-(x-dblMode./dblBeta+exp(-(x-dblMode./dblBeta))));
	
	
	dblStep=0.01;
	vecE=0:dblStep:3;
	vecX = vecE(2:end)-dblStep/2;
	[vecN,e]=histcounts(matRandRate(:),vecE);
	vecN = cumsum(vecN);
	vecN = vecN./sum(vecN(:));
	plot(vecX(2:end),diff(vecN))
	hold on
	vecY = fGumbelCDF(vecX);
	vecY = vecY./sum(vecY(:));
	plot(vecX(2:end),diff(vecY));
	hold off
		%}
	end
end

%% save
if boolSave
	save([strDataTargetPath 'ZetaTime' strRec '.mat' ],...
		'matAnova','matCompTimeAnova','vecTestResamps','matTtest','matZeta','matZetaTime','matCompTimeZeta','matCompTimeZetaTime','intResampNum','strRec');
end
