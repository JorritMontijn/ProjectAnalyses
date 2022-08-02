%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
intResampNum = 500;
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% set variables
strRec = 'SimNoSuperResolution';
dblBaseFreq = 25;%10
intNeurons = 100;
if contains(strRec,'NoSuper')
	dblSamplingFreq = dblBaseFreq;
else
	dblSamplingFreq = dblBaseFreq+(1+rand(1))/10;
end
vecTau = 1;%[2.^(-3:8)];
dblTau0 = (63/1000);
dblNoise = 0.025;
intTauNum = numel(vecTau);
boolQuick = false;

%% pre-allocate output variables
matZeta = nan(intNeurons,2);
matTsZeta = nan(intNeurons,2,intTauNum);
matSrTsZeta = nan(intNeurons,2,intTauNum);
matSrInTsZeta = nan(intNeurons,2,intTauNum);
matAnova = nan(intNeurons,2,intTauNum);
matTtest = nan(intNeurons,2,intTauNum);
	
hTic = tic;
%% generate data
for intNeuron=1:intNeurons
	%% message
	if toc(hTic) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTic=tic;
	end
	
	%% stimulus data
	dblStimDur = 1;
	dblPreBaseDur = 0.25;
	dblPostBaseDur = 0.25;
	dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
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
	dblBaseRate = exprnd(1);
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	
	% generate peak
	dblJitter = 5; %in ms'
	intAddSpikes = round(intTrials*0.2);
	[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes,100/1000);
	
	%% go through sampling freqs
	for intTauIdx=1:intTauNum
		dblTau = vecTau(intTauIdx);
		
		fprintf('Running tau %d/%d [%s]\n',intTauIdx,intTauNum,getTime);
		
		%set indicator properties
		sIndicatorProps = struct;
		sIndicatorProps.dblTimescale = dblTau;
		sIndicatorProps.dblNoise = dblNoise;
		[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick);
		
		%add empty end
		dblEndDur = vecTrialStart(end) + dblTrialDur*5;
		vecAddT = vecTimestamps(end):(1/dblSamplingFreq):dblEndDur;
		vecTimestamps = cat(2,vecTimestamps,vecAddT(2:end));
		vecdFoF = cat(2,vecdFoF,zeros(size(vecAddT(2:end))));
		
		%real+rand
		for intRunType=vecRunTypes
			
			%randomize
			if intRunType ==2
				vecJitter = 2*dblTrialDur*((rand(size(vecTrialStart))-1/2)*2);
				vecUseTrialStart = vecTrialStart + vecJitter -3;
				strRand = 'Rand';
			else
				vecUseTrialStart = vecTrialStart -3;
				strRand = 'Real';
			end
			matUseTrialT = [vecUseTrialStart(:) vecUseTrialStart(:)+dblStimDur];
			
			%plot if first
			if intNeuron == 1 && intRunType == 1
				intPlot = 0;
			else
				intPlot = 0;
			end
			%anova
			dblSampTime = median(diff(vecTimestamps));
			vecWindow = 0:dblSampTime:(dblTrialDur+dblSampTime);
			[matTE,vecWindowBinCenters] = getRespMat(vecTimestamps,vecdFoF,vecUseTrialStart,vecWindow);
			matTrial = repmat((1:(numel(vecWindow)-1)),[numel(vecUseTrialStart) 1]);
			[dblP_Anova,tbl]=anova1(matTE(:),matTrial(:),'off');
			
			%run zeta & ts-zeta
			[dblTsZetaP,sTsZETA] = zetatstest(vecTimestamps,vecdFoF,matUseTrialT,dblTrialDur,intResampNum,intPlot,[],[],0);
			[dblSrTsZetaP,sSrTsZETA] = zetatstest(vecTimestamps,vecdFoF,vecUseTrialStart,dblTrialDur,intResampNum,intPlot,[],[],2);
			[dblSrInTsZetaP,sSrInTsZETA] = zetatstest(vecTimestamps,vecdFoF,vecUseTrialStart,dblTrialDur,intResampNum,intPlot,[],[],1);
			
			
			matAnova(intNeuron,intRunType,intTauIdx) = dblP_Anova;
			dblMeanP = sTsZETA.dblMeanP;
			matTtest(intNeuron,intRunType,intTauIdx) = dblMeanP;
			
			dblTsZeta = sTsZETA.dblZETA;
			matTsZeta(intNeuron,intRunType,intTauIdx) = dblTsZeta;
			dblSrTsZeta = sSrTsZETA.dblZETA;
			matSrTsZeta(intNeuron,intRunType,intTauIdx) = dblSrTsZeta;
			dblSrInTsZeta = sSrInTsZETA.dblZETA;
			matSrInTsZeta(intNeuron,intRunType,intTauIdx) = dblSrInTsZeta;
			%save
			if intTauIdx == 1
				[dblZetaP,sZETA] = zetatest(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,0);
				dblZeta = sZETA.dblZETA;
				matZeta(intNeuron,intRunType) = dblZeta;
			end
			
			if intPlot > 0
				subplot(2,3,5)
				title(sprintf('Samp freq=%.1fHz, Tau=%.3fs, Noise=%.3f, Binned Spikes=%d',...
					dblSamplingFreq,dblTau*dblTau0,dblNoise,boolQuick));
				%subplot(2,3,1)
				%title(sprintf('%s; TS-ZETA=%.3f; ZETA=%.3f',...
				%	strRand,matTsZeta(intNeuron,intRunType,intTauIdx),matZeta(intNeuron,intRunType)));

				%save
				drawnow;
				export_fig(fullpath(strFigPath,[strRec,strRand,sprintf('%.3f',dblTau*dblTau0),'Example.tif']));
				export_fig(fullpath(strFigPath,[strRec,strRand,sprintf('%.3f',dblTau*dblTau0),'Example.pdf']));
				
				if intTauIdx == 1
					[dblZetaP,sZETA] = zetatest(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,3);
					dblZeta = sZETA.dblZETA;
					matZeta(intNeuron,intRunType) = dblZeta;
					hFig=gcf;
					title(hFig.Children(end),strRand);
					
					%save
					drawnow;
					export_fig(fullpath(strFigPath,[strRec,strRand,'ZETA_','Example.tif']));
					export_fig(fullpath(strFigPath,[strRec,strRand,'ZETA_','Example.pdf']));
					
				end
				
			end
		end
	end
end

%% save
if boolSave
	save([strDataTargetPath 'TsZeta' strRec 'Resamp' num2str(intResampNum) 'SampFreq' num2str(round(dblSamplingFreq)) '.mat' ],...
		'matAnova','matTtest','matZeta','matTsZeta','matSrTsZeta','matSrInTsZeta','vecTau','strRec','sIndicatorProps');
end
