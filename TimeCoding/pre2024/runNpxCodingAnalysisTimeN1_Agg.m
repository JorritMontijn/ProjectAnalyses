%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?
%}

%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
boolHome = false;
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx('','natural',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep nm data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	[vecFrameIdx,vecUniqueFrames,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimIdx);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecStimIdx));
	intRepNum = intTrialNum/intStimNum;
	
	%original movie time
	vecOrigStimOnTime = sThisRec.cellBlock{1}.vecStimOnTime;
	vecOrigStimOffTime = sThisRec.cellBlock{1}.vecStimOffTime;
	dblStimDur = min(diff(vecOrigStimOnTime));
	intOrigTrialNum = numel(vecOrigStimOnTime);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%% calculate spike cross-correlation
		%as in https://www.nature.com/articles/ncomms13805
		%C_spike_cross(omega) = C_PSTH_auto(omega)
		% omega is frequency, or equivalently 1/tau, for C(tau)=1/t * SUM_t (r_s(t) * r_s(t+tau)
		%where t is time bin, and r_s(t) is the average over trials
		
		%get PSTHs
		intNumN = numel(cellSpikeTimes);
		
		%% get spike cross-correlations
		intTauSteps = 20;%floor(numel(vecEventBins)/4);
		vecTau = -intTauSteps:intTauSteps;
		vecTimescales = (logspace(log10(1),log10(512),10))./1000;
		matCrossCorr = nan(numel(vecTau),numel(vecTimescales),intRepNum,intNumN);
		matPeakFit = nan(numel(vecTimescales),intNumN,2); %max, width
		matPeakFitR2 = nan(numel(vecTimescales),intNumN,2); %R2, p
		for intTimescale=1:numel(vecTimescales)
			dblTimescale = vecTimescales(intTimescale)
			vecEventBins = 0:dblTimescale:dblStimDur;
			[dummy,dummy,vecEventBinsC,matPET] = doPEP(cellSpikeTimes,vecEventBins,vecOrigStimOnTime);
			
			for intCell=1:intNumN
				intCell
				for intRep=1:intRepNum
					indOtherTrials = true(1,intRepNum);
					indOtherTrials(intRep) = false;
					vecR = matPET(intRep,:,intCell) - mean(matPET(intRep,:,intCell));
					vecOtherR = mean(matPET(indOtherTrials,:,intCell),1) - mean(mean(matPET(indOtherTrials,:,intCell),1));
					dblNormFactor = std(vecR).*std(vecOtherR);
					for intTauIdx=1:numel(vecTau)
						intTau=vecTau(intTauIdx);
						
						matCrossCorr(intTauIdx,intTimescale,intRep,intCell) = (nanmean(vecR .* circshift(vecOtherR,intTau)))./(dblNormFactor);
					end
				end
				
				%fit peak with gaussian
				vecR =nanmean(matCrossCorr(:,intTimescale,:,intCell),3);
				x0 = [max(vecR) numel(vecR)/4];
				fGauss = @(x) (...
					(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2))...
					./max(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2)))...
					)*(x(1)));
				
				fMin = @(x) sum((vecR(:)'-fGauss(x)).^2);
				[x,fval,exitflag,output] = fminsearch(fMin,x0);
				matPeakFit(intTimescale,intCell,:) = x;
				
				[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecR(:)',fGauss(x),2);
				matPeakFitR2(intTimescale,intCell,:) = [dblR2 dblP];
			end
			matAvgCrossCorr = squeeze(nanmean(matCrossCorr(:,intTimescale,:,:),3));
			matNorm = matAvgCrossCorr-nanmean(matAvgCrossCorr,1);
			
			
			%plot
			figure;maxfig;
			intPlotCell = 1;
			vecT = vecTau*dblTimescale;
			
			subplot(2,3,1)
			imagesc(vecT,1:intNumN,matAvgCrossCorr')
			colorbar;
			title(sprintf('timescale %.1fms',dblTimescale*1000));
			xlabel('Time (s)');
			ylabel('Neuron #');
			
			subplot(2,3,2)
			imagesc(vecT,1:intNumN,matNorm')
			colorbar;
			xlabel('Time (s)');
			ylabel('Neuron #');
			
			subplot(2,3,3)
			plot(vecT,matAvgCrossCorr(:,intPlotCell))
			hold on
			vecR=matAvgCrossCorr(:,intPlotCell);
			fGauss = @(x) (...
				(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2))...
				./max(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2)))...
				)*(x(1)));
			
			plot(vecT,fGauss(matPeakFit(intTimescale,intPlotCell,:)));
			xlabel('Time (s)');
			title(sprintf('PSTH CCG; R.^2=%.3f,p=%.3f',matPeakFitR2(intTimescale,intPlotCell,:)));
			ylabel('Cross-correlation (rho)');
			
			vecSds = matPeakFit(intTimescale,:,2)*dblTimescale;
			vecSds(matPeakFitR2(intTimescale,:,2)>0.05)=nan;
		end
		
		%% plot fits
		figure;maxfig;
		matRem = matPeakFitR2(:,:,2)>0.05;
		matR2 = matPeakFitR2(:,:,1);
		matR2(matRem)=nan;
		subplot(2,3,1)
		plot(matR2,'x-')
		
		%take maximum R2 for each cell
		[vecR2,vecIdxT]=max(matR2);
		vecSd = nan(1,intNumN);
		for intCell=1:intNumN
			intIdxT = vecIdxT(intCell);
			vecSd(intCell)=vecTimescales(intIdxT)*matPeakFit(intIdxT,intCell,2);
		end
		subplot(2,3,2)
		vecBins = 2.^(4:0.5:10);
		vecBinsC = vecBins(2:end)-diff(vecBins);
		[vecCounts]=histcounts(vecSd,vecBins/1000);
		plot(vecBinsC,vecCounts);
		set(gca,'xscale','log')
		xlabel('Temporal precision (ms)');
		ylabel('# of neurons');
		
		%% run spike-time based decoding
		%create high-resolution binned matrix and 
		dblBinDur = 1/1000;
		intBinsPerRep = round(dblStimDur/dblBinDur);
		vecBinEdges = linspace(0,dblBinDur*intBinsPerRep,intBinsPerRep+1);
		
		%generate fake stimulus vectors
		vecUniqueStims = 1:intBinsPerRep;
		vecFrameIdx = repmat(vecUniqueStims,[1 numel(vecOrigStimOnTime)]);
		vecStimOnTime = flat(vecBinEdges(1:(end-1))' + vecOrigStimOnTime)';
		vecStimOffTime = vecStimOnTime + dblBinDur;
		[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecFrameIdx);
		
		%get matrix
		matSpikeCounts = getSpikeCounts(cellSpikeTimes,vecStimOnTime,vecStimOffTime);
		matSpikeCounts(matSpikeCounts>0)=1;
		
		%blur single spikes with time constants from above
		for intN=1:intNumN
			intN
			dblSd = vecSd(intN)/dblBinDur;
			vecRange = round(-2*dblSd):round(2*dblSd);
			gVecS = gpuArray(matSpikeCounts(intN,:));
			gVecFilt = gpuArray(normpdf(vecRange,0,dblSd));
			
			gVecS = padarray(gVecS,floor(size(gVecFilt)/2),'circular');
	
			gVecSmoothed = conv(gVecS,gVecFilt,'valid');
			
			matSpikeCounts(intN,:) = gather(gVecSmoothed);
		end
		%run standard decoder
		error out of memory
		intTypeCV = 2;
		dblLambda = 1;
		matData = single(matSpikeCounts);
		vecTrialTypes = single(vecStimIdx);
		[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion] = ...
			doCrossValidatedDecodingLR(matData(1,:),vecTrialTypes,intTypeCV,[],dblLambda);
		
		%% repeat for whole population
		matCrossCorr_Pop = nan(numel(vecTau),numel(vecTimescales),intRepNum);
		matPeakFit_Pop = nan(numel(vecTimescales),2); %max, width
		matPeakFitR2_Pop = nan(numel(vecTimescales),2); %R2, p
		vecPopSpikes = sort(cell2vec(cellSpikeTimes));
		vecPopSpikes(vecPopSpikes < vecOrigStimOnTime(1) | vecPopSpikes > vecOrigStimOnTime(end)+dblStimDur) = [];
		for intTimescale=1:numel(vecTimescales)
			dblTimescale = vecTimescales(intTimescale)
			vecEventBins = 0:dblTimescale:dblStimDur;
			[dummy,dummy,vecEventBinsC,matPET] = doPEP(vecPopSpikes,vecEventBins,vecOrigStimOnTime,-1);
			
			
			for intRep=1:intRepNum
				indOtherTrials = true(1,intRepNum);
				indOtherTrials(intRep) = false;
				vecR = matPET(intRep,:) - mean(matPET(intRep,:));
				vecOtherR = mean(matPET(indOtherTrials,:),1) - mean(mean(matPET(indOtherTrials,:),1));
				dblNormFactor = std(vecR).*std(vecOtherR);
				for intTauIdx=1:numel(vecTau)
					intTau=vecTau(intTauIdx);
					
					matCrossCorr_Pop(intTauIdx,intTimescale,intRep) = (nanmean(vecR .* circshift(vecOtherR,intTau)))./(dblNormFactor);
				end
			end
			
			%fit peak with gaussian
			vecR =nanmean(matCrossCorr_Pop(:,intTimescale,:),3);
			x0 = [max(vecR) numel(vecR)/4];
			fGauss = @(x) (...
				(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2))...
				./max(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2)))...
				)*(x(1)));
			
			fMin = @(x) sum((vecR(:)'-fGauss(x)).^2);
			[x,fval,exitflag,output] = fminsearch(fMin,x0);
			matPeakFit_Pop(intTimescale,:) = x;
			
			[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecR(:)',fGauss(x),2);
			matPeakFitR2_Pop(intTimescale,:) = [dblR2 dblP];
			
			matAvgCrossCorr_Pop = squeeze(nanmean(matCrossCorr_Pop(:,intTimescale,:,:),3));
			matNorm_Pop = matAvgCrossCorr_Pop-nanmean(matAvgCrossCorr_Pop,1);
			
			
			%plot
			figure;maxfig;
			intPlotCell = 1;
			vecT = vecTau*dblTimescale;
			
			subplot(2,3,3)
			plot(vecT,matAvgCrossCorr_Pop(:,intPlotCell))
			hold on
			vecR=matAvgCrossCorr_Pop(:,intPlotCell);
			fGauss = @(x) (...
				(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2))...
				./max(normpdf(1:numel(vecR),ceil(numel(vecR)/2),x(2)))...
				)*(x(1)));
			
			plot(vecT,fGauss(matPeakFit_Pop(intTimescale,:)));
			xlabel('Time (s)');
			title(sprintf('PSTH CCG; R.^2=%.3f,p=%.3f',matPeakFitR2_Pop(intTimescale,:)));
			ylabel('Cross-correlation (rho)');
		end
		
		%% plot fits
		figure;maxfig;
		matRem_Pop = matPeakFitR2_Pop(:,2)>0.05;
		matR2_Pop = matPeakFitR2_Pop(:,1);
		matR2_Pop(matRem_Pop)=nan;
		subplot(2,3,1)
		plot(matR2_Pop,'x-')
		
		%take maximum R2 for each cell
		[dblR2_Pop,intIdxT_Pop]=max(matR2_Pop);
		dblPopSd=vecTimescales(intIdxT_Pop)*matPeakFit_Pop(intIdxT_Pop,2);
		subplot(2,3,2)
		vecBins = 2.^(4:0.5:10);
		vecBinsC = vecBins(2:end)-diff(vecBins);
		[vecCounts]=histcounts(dblPopSd,vecBins/1000);
		plot(vecBinsC,vecCounts);
		set(gca,'xscale','log')
		xlabel('Temporal precision (ms)');
		ylabel('# of neurons');
		
	end
end
toc
