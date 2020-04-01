%load data
vecFrac=[];
for intMouse=5
	close all
	drawnow;
	clearvars -except vecFrac intMouse
	
	%use neuropil subtraction?
	boolUseNeuropilSubtraction = false;
	boolOnlyTuned = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == 4
			strSes = '20140507';
		elseif intMouse == 5
			strSes = '20140530';
		elseif intMouse == 6
			strSes = '20140604';
		elseif intMouse == 7
			strSes = '20140711';
		elseif intMouse == 8
			strSes = '20140715';
		end
		if boolUseNeuropilSubtraction
			strSes = ['NPS' strSes];
		end
		if boolOnlyTuned
			strSes = [strSes 'OT'];
		end
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['E:\UvA_Backup\Data\Results\spikeAnalysis\dataPreProAssemblies' strSes '.mat']);
	end
	
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['E:\UvA_Backup\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = false;
	sParams.boolSaveData = false;
	%strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	%#ok<*AGROW>
	strSes = ['AA' strSes]; 
	
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
		%{
		vecFrac(end+1) = (numel(sAssemblies{intPopulation}.vecAssemblies)/numel(cellMultiSes{intPopulation}.neuron(1).dFoF)/cellMultiSes{intPopulation}.samplingFreq)*1000
		if intPopulation==2
			vecFrac(end-1) = (vecFrac(end-1) + vecFrac(end))/2;
			vecFrac(end) = []
		end
		%}
		
		%define stimulus duration parameter
		dblStimSecs = 2;
		
		% remove trials with reaction time <150ms
		dblRemSecs = 0.15;
		indTooFast = (cellMultiSes{intPopulation}.structStim.vecTrialRespSecs-cellMultiSes{intPopulation}.structStim.SecsOn)<dblRemSecs & cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		cellMultiSes{intPopulation}.structStim = remel(cellMultiSes{intPopulation}.structStim,~indTooFast);
		
		%get stimulus variables
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Contrast'});
		cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		intContrasts = length(cellSelectC);
		
		sTypesO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation'});
		cellSelectO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesO);
		intOrientations = length(cellSelectO);
		
		sTypesCO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation','Contrast'});
		cellSelectCO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesCO);
		
		% get stimulus data
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		intOris = length(vecOrientations);
		[vecStimOris,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Contrast,unique(cellMultiSes{intPopulation}.structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn<3;
		
		
		%get behavioral variables
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs -  cellMultiSes{intPopulation}.structStim.SecsOn;
		vecRTs(indMiss) = nan;
		indFast = false(size(indMiss));
		indSlow = false(size(indMiss));
		for intContrast=unique(vecStimContrasts)
			indTheseContrastTrials = vecStimContrasts==intContrast;
			indFast = (vecRTs <= nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indFast;
			indSlow = (vecRTs > nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indSlow;
		end
		
		%get assembly variables
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		%matRespAssemblies = sAssemblies{intPopulation}.matRespAssemblies; 
		matReorderedAssemblyCorrs=sAssemblies{intPopulation}.matAssemblyCorrelations(sAssemblies{intPopulation}.vecReorder,sAssemblies{intPopulation}.vecReorder);
		matM = getBlockMeans(matReorderedAssemblyCorrs,sAssemblies{intPopulation}.vecAssemblyIdentity);
		
		%get assembly activity during trials
		dblSamplingFreq = (cellMultiSes{intPopulation}.structStim.FrameOff(end)-cellMultiSes{intPopulation}.structStim.FrameOff(1))/(cellMultiSes{intPopulation}.structStim.SecsOff(end)-cellMultiSes{intPopulation}.structStim.SecsOff(1));
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		sStimTemp.FrameOff = sStimTemp.FrameOn + round(dblSamplingFreq*dblStimSecs);
		%sStimTemp.FrameOff = sStimTemp.FrameOn + round(((sStimTemp.FrameOff- sStimTemp.FrameOn)/2)*dblStimSecs);
		matRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp);%[assemblies x trials]
		
		%get response matrix
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},sStimTemp);
		
		%get pre-resp assemblies
		structParamsAss.intStopOffset = 0;
		structParamsAss.intStartOffset = round(-1*cellMultiSes{intPopulation}.samplingFreq);
		matPreRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp,[],[],structParamsAss);%[assemblies x trials]

		%detect spikes
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		dblTau = 0.5;
		intTotDur = length(cellMultiSes{intPopulation}.neuron(1).dFoF);
		dblTotDurSecs = intTotDur/dblSamplingFreq;
		if isfield(sAssemblies{intPopulation},'matSpikeCountsHD')
			matSpikeCountsHD = sAssemblies{intPopulation}.matSpikeCountsHD;
		else
			matSpikeCountsHD = nan(intNeurons,intTotDur);
			for intNeuron=1:intNeurons
				%msg
				fprintf('Performing detection for neuron %d/%d...\n',intNeuron,intNeurons)
				ptrTime = tic;
				
				%AE detection
				[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(cellMultiSes{intPopulation}.neuron(intNeuron).dFoF,dblSamplingFreq,dblTau);
				intTotSpikes = sum(apSpikes);
				
				%assign to neuron
				cellMultiSes{intPopulation}.neuron(intNeuron).apFrames = apFrames;
				cellMultiSes{intPopulation}.neuron(intNeuron).apSpikes = apSpikes;
				cellMultiSes{intPopulation}.neuron(intNeuron).vecSpikes = vecSpikes;
				cellMultiSes{intPopulation}.neuron(intNeuron).expFit = expFit;
				matSpikeCountsHD(intNeuron,:) = vecSpikes;
				
				%msg
				fprintf('\b	Done! Took %.1f seconds; %d transients; mean rate is %.2f Hz\n',toc(ptrTime),intTotSpikes,intTotSpikes/dblTotDurSecs)
			end
		end
		matSpikeCounts = getNeuronResponseRespMat(matSpikeCountsHD,sStimTemp);
		matSpikeCounts = round(matSpikeCounts*dblSamplingFreq);
		
		%% create miss/slow/fast selection vectors
		intFrames = size(matSpikeCountsHD,2);
		%miss
		vecMissOn = sStimTemp.FrameOn(indMiss);
		vecMissOff = vecMissOn + round(dblSamplingFreq*dblStimSecs);
		indMissHD = false(1,intFrames);
		for intTrial=1:length(vecMissOn)
			indMissHD(vecMissOn(intTrial):vecMissOff(intTrial)) = true;
		end
		vecMissHD = find(indMissHD);
		%slow
		vecSlowOn = sStimTemp.FrameOn(indSlow);
		vecSlowOff = vecSlowOn + round(dblSamplingFreq*dblStimSecs);
		indSlowHD = false(1,intFrames);
		for intTrial=1:length(vecSlowOn)
			indSlowHD(vecSlowOn(intTrial):vecSlowOff(intTrial)) = true;
		end
		vecSlowHD = find(indSlowHD);
		%fast
		vecFastOn = sStimTemp.FrameOn(indFast);
		vecFastOff = vecFastOn + round(dblSamplingFreq*dblStimSecs);
		indFastHD = false(1,intFrames);
		for intTrial=1:length(vecFastOn)
			indFastHD(vecFastOn(intTrial):vecFastOff(intTrial)) = true;
		end
		vecFastHD = find(indFastHD);
		
		%% get members of assemblies
		intNeurons = size(sAssemblies{intPopulation}.matAssemblies,1);
		matAssemblyMembers = sAssemblies{intPopulation}.matAssemblies > 0;
		intAssemblies = max(sAssemblies{intPopulation}.vecAssemblies);
		matAssemblyCoreMembers = false(intAssemblies,intNeurons);
		for intAssembly=1:intAssemblies
			%get random memberships
			intRandIters = 1000;
			indAssembly = sAssemblies{intPopulation}.vecAssemblies == intAssembly;
			matRandMembership = nan(intNeurons,intRandIters);
			for intIter=1:intRandIters
				indRandomAssembly = indAssembly(randperm(length(indAssembly)));
				matRandMembership(:,intIter) = sum(matAssemblyMembers(:,indRandomAssembly),2)/sum(indAssembly);
			end
			dblAlpha = 0.01;
			dblCorrectedSDs = qnorm(1-(dblAlpha/2)/intNeurons);
			vecUpperBound = mean(matRandMembership,2)+std(matRandMembership,[],2)*dblCorrectedSDs;
			vecData = sum(matAssemblyMembers(:,indAssembly),2)/sum(indAssembly);
			indHigher = vecData>vecUpperBound;
			
			%save output
			matAssemblyCoreMembers(intAssembly,:) = indHigher;
		end
		%matAssemblyCoreMembers = true(size(matAssemblyCoreMembers));
		%get assembly types
		vecCoreMemberSize = sum(matAssemblyCoreMembers,2);
		vecAssemblyAutoCorr = diag(matM);
		indRecurringClusters = zscore(vecAssemblyAutoCorr)>-2;
		[dummy,intNonRecurringCluster]=min(zscore(vecAssemblyAutoCorr));
		indRecurringClusters(intNonRecurringCluster) = false;
		indMultiNeuronClusters = vecCoreMemberSize>1;
		indRealAssemblies = indRecurringClusters & indMultiNeuronClusters;
		%indRealAssemblies = true(size(indRealAssemblies));
		
		%% get assembly data
		vecAssemblyStarts = sAssemblies{intPopulation}.vecAssemblyStarts;
		vecAssemblyStops = sAssemblies{intPopulation}.vecAssemblyStops;
		vecAssemblies = sAssemblies{intPopulation}.vecAssemblies';
		vecRealAssemblies = ismember(vecAssemblies,find(indRealAssemblies));
		vecNonRecurringAssemblies = vecAssemblies == intNonRecurringCluster;
		vecAssemblyAEs = zeros(size(vecAssemblyStarts));
		dblSecsCCG = (max(vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq);
		vecPlotCCG = (-ceil(dblSamplingFreq*dblSecsCCG):ceil(dblSamplingFreq*dblSecsCCG));
		intMiddle = ceil(dblSamplingFreq)+1;
		
		%% assembly occurrence probability drop-off (PSTH relative to occurrences, or to stimulus onset); self vs. other
		close all
		dblSamplingFreq = 25.4;
		dblNumSecsPlot = 60;
		vecPSTH_Frames1 = -round(dblNumSecsPlot*dblSamplingFreq):round(dblNumSecsPlot*dblSamplingFreq);
		cellPSTH_self = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_others = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_Stim = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_Offset = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		vecLaggingRatio = nan(1,intAssemblies);
		vecStimDrivenRatio = nan(1,intAssemblies);
		for intAssembly=find(indRealAssemblies)'
			vecOccurrences = matAssemblyActivity(intAssembly,:);
			indOtherAssemblies = true(1,intAssemblies);
			indOtherAssemblies(intAssembly) = false;
			vecOtherOccurrences = sum(matAssemblyActivity(indOtherAssemblies,:),1);
			for intOccFrame=find(logical(vecOccurrences))
				vecGetFrames = vecPSTH_Frames1+intOccFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_self{intAssembly}(vecAssgnFrames) = cellPSTH_self{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
				cellPSTH_others{intAssembly}(vecAssgnFrames) = cellPSTH_others{intAssembly}(vecAssgnFrames) + vecOtherOccurrences(vecGetFrames);
			end
			
			vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn;
			for intStimFrame=vecStimOn
				vecGetFrames = vecPSTH_Frames1+intStimFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_Stim{intAssembly}(vecAssgnFrames) = cellPSTH_Stim{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
			end
			
			vecStimOff = cellMultiSes{intPopulation}.structStim.FrameOff;
			for intStimFrame=vecStimOff(indHit)
				vecGetFrames = vecPSTH_Frames1+intStimFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_Offset{intAssembly}(vecAssgnFrames) = cellPSTH_Offset{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
			end
			
			%get data
			vecSelf = cellPSTH_self{intAssembly}/sum(vecOccurrences);
			vecOthers = cellPSTH_others{intAssembly}/sum(vecOccurrences);
			dblCorrFactorOthers = sum(vecOtherOccurrences)/length(vecOtherOccurrences);
			dblCorrFactorSelf = sum(vecOccurrences)/length(vecOccurrences);
			[vecPowerSelf,vecF] = pwelch(vecSelf/dblCorrFactorSelf,[],[],length(vecSelf),dblSamplingFreq,'onesided');
			[vecPowerOthers,vecFOthers] = pwelch(vecOthers/dblCorrFactorOthers,[],[],length(vecOthers),dblSamplingFreq,'onesided');
			
			intMid = ceil(length(vecOthers)/2);
			vecPre = (intMid-round(dblSamplingFreq)):(intMid-1);
			vecPost = (intMid+1):(intMid+round(dblSamplingFreq));
			vecData = vecOthers/dblCorrFactorOthers;
			vecPreData = vecData(vecPre);
			vecPostData = vecData(vecPost);
			
			dblLaggingRatio = mean(vecPreData)/mean(vecPostData);
			vecLaggingRatio(intAssembly) = dblLaggingRatio;
			
			vecSelfStim = cellPSTH_Stim{intAssembly};
			vecSelfOffset = cellPSTH_Offset{intAssembly};
			intStartStim = ceil(length(vecSelfStim)/2);
			intStimDur = round(dblSamplingFreq*dblStimSecs);
			vecStimDur = intStartStim:(intStartStim+intStimDur-1);
			dblBaseAct = mean(vecSelfStim(1:(intStartStim-1)));
			dblBaseActOffset = mean(vecSelfOffset(1:(intStartStim-1)));
			dblMeanStim = mean(vecSelfStim(vecStimDur));
			vecStimDrivenRatio(intAssembly) = dblMeanStim/dblBaseAct;
			
			%save figure
			if sParams.boolSavePlots
				figure
				subplot(2,2,1)
				plot([-dblNumSecsPlot dblNumSecsPlot],[1 1],'k--');
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecOthers/dblCorrFactorOthers,'b')
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelf/dblCorrFactorSelf,'r')
				
				hold off;
				title(sprintf('Assembly %d; Blue=others, red=self',intAssembly));
				set(gca,'yscale','log');
				xlim([-dblNumSecsPlot dblNumSecsPlot])
				ylim([0.5 10]);
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				subplot(2,2,2)
				pwelch(vecSelf/dblCorrFactorSelf,[],[],length(vecSelf),dblSamplingFreq,'onesided');
				title('Power spectrum self-occurrence (Welch)')
				
				%{
				subplot(2,2,3)
				errorbar(1:2,[mean(vecPreData) mean(vecPostData)],[std(vecPreData)/sqrt(dblSamplingFreq) std(vecPostData)/sqrt(dblSamplingFreq)],'x');
				title(sprintf('Lagging ratio, pre/post: %.3f',dblLaggingRatio));
				ylabel('Occurrence probability');
				set(gca,'xtick',[1 2],'xticklabel',{'Pre-occurrence','Post-occurrence'});
				%}
				subplot(2,2,3)
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelfOffset/dblBaseAct,'r')
				hold off;
				title(sprintf('centered on stim offset'));
				%set(gca,'yscale','log');
				xlim([-8 8])
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				
				
				subplot(2,2,4)
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelfStim/dblBaseAct,'r')
				hold off;
				title(sprintf('centered on stim onset'));
				%set(gca,'yscale','log');
				xlim([-5 10])
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				export_fig(sprintf('%spop%d_ass%doccurrence_probability.tif',strSes,intPopulation,intAssembly));
				export_fig(sprintf('%spop%d_ass%doccurrence_probability.pdf',strSes,intPopulation,intAssembly));
			end
			
			%assign data
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,1} = vecF;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,2} = vecPowerSelf;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,3} = vecPowerOthers;
			
			%occurrence probability centered on event
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,4} = vecPSTH_Frames1/dblSamplingFreq;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,5} = vecSelf/dblCorrFactorSelf;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,6} = vecOthers/dblCorrFactorOthers;
			
			%occurrence probability centered on stim presentation
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,7} = vecSelfStim;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,8} = dblBaseAct;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,9} = vecSelfOffset;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,10} = dblBaseActOffset;
		end
		continue
		
		%% perform decoding; are stimuli more accurately represented during population events?
		indStimTrials = ~cellSelectContrasts{1};
		vecOriStimTrials = vecStimOris(indStimTrials);
		vecContrastStimTrials = vecStimContrasts(indStimTrials);
		matRespRealAssemblies = matRespAssemblies(indRealAssemblies,:);
		matPreRespRealAssemblies  = matPreRespAssemblies(indRealAssemblies,:);
		matRespNonRecurring =  matRespAssemblies(~indRealAssemblies,:);
		matPreRespNRAssemblies = matPreRespAssemblies(~indRealAssemblies,:);
		intStimTrials = length(vecOriStimTrials);
		
		cellSelectSubO = cell(1,4);
		for i=1:length(cellSelectO)
			cellSelectSubO{i} = cellSelectO{i}(indStimTrials);
		end
		
		
		%decode orientation with cross-validated template matching algorithm
		structStimDecode = remel(sStimTemp,indStimTrials);
		sMetaOut_dFoF = doCrossValidationRespMat(matTrialResponse(:,indStimTrials),structStimDecode);
		dblPerformanceML_dFoF = sum(sMetaOut_dFoF.vecDecodedStimType==sMetaOut_dFoF.vecStimType)/sum(indStimTrials);
		sMetaOut_AEs = doCrossValidationRespMat(matSpikeCounts(:,indStimTrials),structStimDecode);
		dblPerformanceML_AEs = sum(sMetaOut_AEs.vecDecodedStimType==sMetaOut_AEs.vecStimType)/sum(indStimTrials);
		[phat,dblCI_MLdFoF] = binofit(dblPerformanceML_dFoF*intStimTrials,intStimTrials);
		[phat,dblCI_MLAEs] = binofit(dblPerformanceML_AEs*intStimTrials,intStimTrials);
		

		% decode orientation with template matching
		[dblPerformance_dFoF,vecDecodedIndexCV,matDists] = doCrossValidatedTemplateMatching(matTrialResponse(:,indStimTrials),vecOriStimTrials);
		[dblPerformance_AEs,vecDecodedIndexCV_AEs_TM,matDists] = doCrossValidatedTemplateMatching(matSpikeCounts(:,indStimTrials),vecOriStimTrials);
		[phat,dblCI_dFoF] = binofit(dblPerformance_dFoF*intStimTrials,intStimTrials);
		[phat,dblCI_AEs] = binofit(dblPerformance_AEs*intStimTrials,intStimTrials);
		
		%get neuron response properties
		structStimCorrsNAE = calcStimCorrsRespMat(matSpikeCounts(:,indStimTrials),cellSelectSubO);
		sTuningNAE = calcTuningRespMat(matSpikeCounts(:,indStimTrials),cellSelectSubO,vecOrientations);
		
		%get assembly response properties
		structStimCorrsAss = calcStimCorrsRespMat(matRespRealAssemblies(:,indStimTrials),cellSelectSubO);
		sTuningAss = calcTuningRespMat(matRespRealAssemblies(:,indStimTrials),cellSelectSubO,vecOrientations);
		
		%plot decoding
		if sParams.boolSavePlots
			figure
			%TM
			subplot(2,2,1)
			vecDecodingPerformance = [dblPerformance_dFoF dblPerformance_AEs];
			errorbar([0.3 0.7],vecDecodingPerformance,[dblCI_dFoF(1) dblCI_AEs(1)]-vecDecodingPerformance,[dblCI_dFoF(2) dblCI_AEs(2)]-vecDecodingPerformance,'xb')
			xlim([0 1])
			ylim([0 1])
			hold on
			plot(get(gca,'xlim'),[1 1]/length(cellSelectSubO),'k--')
			set(gca,'XTick',[0.3 0.7],'xticklabel',{'Neuron dF/F0','Neuron AE'})
			ylabel('Orientation decoding accuracy (TM)')
			
			%ML
			subplot(2,2,2)
			vecDecodingPerformance = [dblPerformanceML_dFoF dblPerformanceML_AEs];
			errorbar([0.3 0.7],vecDecodingPerformance,[dblCI_MLdFoF(1) dblCI_MLAEs(1)]-vecDecodingPerformance,[dblCI_MLdFoF(2) dblCI_MLAEs(2)]-vecDecodingPerformance,'xb')
			xlim([0 1])
			ylim([0 1])
			hold on
			plot(get(gca,'xlim'),[1 1]/length(cellSelectSubO),'k--')
			set(gca,'XTick',[0.3 0.7],'xticklabel',{'Neuron dF/F0','Neuron AE'})
			ylabel('Orientation decoding accuracy (ML)')
		end
		
		%% get dependence on nr of PEs
		%get population events at each trial
		vecRecurringTrialPEs = nan(1,intTrials);
		vecNonRecurringTrialPEs = nan(1,intTrials);
		vecAllTrialPEs = nan(1,intTrials);
		
		for intTrial = 1:intTrials
			intStart = sStimTemp.FrameOn(intTrial);
			intStop = sStimTemp.FrameOff(intTrial);
			
			vecPEs = (vecAssemblyStarts > intStart & vecAssemblyStarts < intStop) | (vecAssemblyStops > intStart & vecAssemblyStops < intStop);
			vecRecurringTrialPEs(intTrial) = sum(vecRealAssemblies(vecPEs));
			vecNonRecurringTrialPEs(intTrial) = sum(vecNonRecurringAssemblies(vecPEs));
			vecAllTrialPEs(intTrial) = sum(vecPEs);
		end
		
		%transform to correct
		vecDecodedCorrect = vecDecodedIndexCV_AEs_TM == vecOriStimTrials;
		
		%group by how many PEs
		intMaxPEs = 4;
		vecRecurAccMean = nan(1,intMaxPEs+2);
		matRecurAccCI = nan(2,intMaxPEs+2);
		vecNonRecurAccMean = nan(1,intMaxPEs+2);
		matNonRecurAccCI = nan(2,intMaxPEs+2);
		vecAllAccMean = nan(1,intMaxPEs+2);
		matAllAccCI = nan(2,intMaxPEs+2);
		for intPEs = 0:(intMaxPEs+1)
			%recurring
			if intPEs == (intMaxPEs+1)
				vecTheseTrials = vecRecurringTrialPEs(indStimTrials) > intMaxPEs;
			else
				
				vecTheseTrials = vecRecurringTrialPEs(indStimTrials) == intPEs;
			end
			vecRecurAccMean(intPEs+1) = mean(vecDecodedCorrect(vecTheseTrials));
			[phat,matRecurAccCI(:,intPEs+1)] = binofit(vecRecurAccMean(intPEs+1)*sum(vecTheseTrials),sum(vecTheseTrials),0.32);
			
			%non-recurring
			if intPEs == (intMaxPEs+1)
				vecTheseTrials = vecNonRecurringTrialPEs(indStimTrials) > intMaxPEs;
			else
				
				vecTheseTrials = vecNonRecurringTrialPEs(indStimTrials) == intPEs;
			end
			vecNonRecurAccMean(intPEs+1) = mean(vecDecodedCorrect(vecTheseTrials));
			[phat,matNonRecurAccCI(:,intPEs+1)] = binofit(vecNonRecurAccMean(intPEs+1)*sum(vecTheseTrials),sum(vecTheseTrials),0.32);
			
			%all
			if intPEs == (intMaxPEs+1)
				vecTheseTrials = vecAllTrialPEs(indStimTrials) > intMaxPEs;
			else
				
				vecTheseTrials = vecAllTrialPEs(indStimTrials) == intPEs;
			end
			vecAllAccMean(intPEs+1) = mean(vecDecodedCorrect(vecTheseTrials));
			[phat,matAllAccCI(:,intPEs+1)] = binofit(vecAllAccMean(intPEs+1)*sum(vecTheseTrials),sum(vecTheseTrials),0.32);
		end
		%save data
		cellSaveDASPE{intPopulation,1} = [vecRecurAccMean;vecNonRecurAccMean;vecAllAccMean];
		
		%save figure
		if sParams.boolSavePlots
			%plot
			vecX = 0:(intMaxPEs+1);
			vecLimX = [-0.5 intMaxPEs+1.5];
			cellTicks = num2cell(num2str(vecX'));
			cellTicks{end} = ['>' cellTicks{end-1}];
			cellTicks = cellTicks';
			
			figure
			subplot(2,2,1)
			plot(vecLimX,[1 1]/4,'k--');
			hold on
			errorbar(vecX,vecRecurAccMean,matRecurAccCI(1,:)-vecRecurAccMean,matRecurAccCI(2,:)-vecRecurAccMean)
			hold off
			ylim([0 1]);
			xlim(vecLimX);
			set(gca,'xtick',0:7,'xticklabel',cellTicks);
			xlabel('Number of PEs during stim. pres.');
			ylabel('Decoding accuracy');
			title('Recurring PEs');
			
			subplot(2,2,2)
			plot(vecLimX,[1 1]/4,'k--');
			hold on
			errorbar(vecX,vecNonRecurAccMean,matNonRecurAccCI(1,:)-vecNonRecurAccMean,matNonRecurAccCI(2,:)-vecNonRecurAccMean)
			hold off
			ylim([0 1]);
			xlim(vecLimX);
			set(gca,'xtick',vecX,'xticklabel',cellTicks);
			xlabel('Number of PEs during stim. pres.');
			ylabel('Decoding accuracy');
			title('Non-recurring PEs');
			
			subplot(2,2,3)
			plot(vecLimX,[1 1]/4,'k--');
			hold on
			errorbar(vecX,vecAllAccMean,matAllAccCI(1,:)-vecAllAccMean,matAllAccCI(2,:)-vecAllAccMean)
			hold off
			ylim([0 1]);
			xlim(vecLimX);
			set(gca,'xtick',vecX,'xticklabel',cellTicks);
			xlabel('Number of PEs during stim. pres.');
			ylabel('Decoding accuracy');
			title('All PEs');
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dDecAccStimPEs.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dDecAccStimPEs.pdf',strSes,intPopulation));
		end
		
		%% get latency matrix [neuron x PE]
		intAssOccs = numel(vecAssemblyStarts);
		matLatencies = nan(intNeurons,intAssOccs);
		vecLatencySDs = nan(1,intAssOccs);
		vecAssDuringTrial = nan(1,intAssOccs);
		vecAssPreTrial = nan(1,intAssOccs);
		intPreFrames = round(2*dblSamplingFreq);
		for intOccurrence=1:intAssOccs
			%get spike times for this occurrence
			if ~indMultiNeuronClusters(vecAssemblies(intOccurrence)),continue;end
			intStart = vecAssemblyStarts(intOccurrence);
			intStop = vecAssemblyStops(intOccurrence)-1;
			matAEs = matSpikeCountsHD(:,intStart:intStop);
			
			%get mean and std
			matVals = repmat(1:size(matAEs,2),[size(matAEs,1) 1]);
			vecVals = getValByIdx(matVals,matAEs);
			dblMean = mean(vecVals);
			dblSD = std(vecVals);
			
			%get overall sd
			matLatencies(:,intOccurrence) = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
			vecLatencySDs(intOccurrence) = dblSD;
			
			%get trial
			vecDuringTrials = find((sStimTemp.FrameOn > intStart & sStimTemp.FrameOn < intStop) |...
				(sStimTemp.FrameOff > intStart & sStimTemp.FrameOff < intStop));
			if ~isempty(vecDuringTrials)
				vecAssDuringTrial(intOccurrence) = vecDuringTrials;
			end
			
			%get pre-trial
			vecPreTrials = find((sStimTemp.FrameOn-intPreFrames > intStart & sStimTemp.FrameOn-intPreFrames < intStop) |...
				(sStimTemp.FrameOff-intPreFrames > intStart & sStimTemp.FrameOff-intPreFrames < intStop));
			if ~isempty(vecPreTrials)
				vecAssPreTrial(intOccurrence) = vecPreTrials;
			end
		end
		%matLatencies = (matLatencies / cellMultiSes{intPopulation}.samplingFreq)*1000;
		
		%get correlation with overall
		[h,vecSignificantOffset]=ttest(matLatencies');
		vecMeanOffset = nanmean(matLatencies,2);
		vecCorr = nan(1,intAssOccs);
		for intOccurrence=1:intAssOccs
			vecCorr(intOccurrence) = corr(vecMeanOffset,matLatencies(:,intOccurrence),'rows','pairwise','type','Pearson');
		end
		
		%%
		
		close all
		[vecMeanSort,vecSort] = sort(vecMeanOffset);
		matSortedLatencies = matLatencies(vecSort,:);
		vecMeanSort = vecMeanSort-mean(vecMeanSort);
		
		nancolorbar(matSortedLatencies,[],redbluepurple);
		subplot(2,2,1)
		scatter(vecMeanSort,1:numel(vecMeanSort),[],vecMeanSort,'filled')
		colorbar
		caxis([-0.2 0.2])
		xlim([-0.7 0.7])
		ylim([0 intNeurons])
		axis ij
		
		
		
		%%
		close all
		figure
		colormap(redbluepurple)
		intFigC = 0;
		vecDurTrials = find(~isnan(vecAssDuringTrial));
		intOcc = 0;
		cellRespType = {'Miss','Slow','Fast'};
		%vecDurTrials = [9 33 237 604];
		vecDurTrials = [255 2291 1034 2071 1557 1440 2116 3307];
		
		%hit: 255 2291 1034 2071
		%miss: 1557 1440 2116 3307
		while intOcc<length(vecDurTrials)
			intOcc = intOcc + 1;
			intOccTrial = vecDurTrials(intOcc);
			intStartFrame = vecAssemblyStarts(intOccTrial);
			intTrial = find(sStimTemp.FrameOn<intStartFrame,1,'last');
			intResp =  indSlow(intTrial)+indFast(intTrial)*2;
			intFigC = intFigC + 1;
			subplot(2,4,intFigC)
			
			scatter(matSortedLatencies(:,intOccTrial),1:numel(vecMeanSort),[],matSortedLatencies(:,intOccTrial),'filled')
			axis ij
			caxis([-1.5 1.5])
			xlim([-2 2])
			ylim([0 intNeurons])
			%nancolorbar(matSortedLatencies(:,intOccTrial),[],redbluepurple);
			
			title(sprintf('Trial %d, %s; occ %d, corr: %.3f',intTrial,cellRespType{intResp+1},intOccTrial,vecCorr(intOccTrial)))
			if intFigC == 8
				intFigC = 0;
				pause;
			end
		end
		
		%}
		
		return
		%% get variability per trial
		vecConsistencyPerTrial = nan(1,intTrials);
		vecPreConsistencyPerTrial = nan(1,intTrials);
		for intTrial=1:intTrials
			%get consistency
			vecOccurrences = find(vecAssDuringTrial==intTrial);
			if ~isempty(vecOccurrences)
				vecConsistencyPerTrial(intTrial) = nanmean(vecCorr(vecOccurrences));
			end
			
			%get pre-consistency
			vecPreOccurrences = find(vecAssPreTrial==intTrial);
			if ~isempty(vecPreOccurrences)
				vecPreConsistencyPerTrial(intTrial) = nanmean(vecCorr(vecPreOccurrences));
			end
		end
		%get consistency per response type
		vecHitConsistency = vecConsistencyPerTrial(indHit);
		vecMissConsistency = vecConsistencyPerTrial(indMiss);
		vecSlowConsistency = vecConsistencyPerTrial(indSlow);
		vecFastConsistency = vecConsistencyPerTrial(indFast);
		vecPreConsistency = vecPreConsistencyPerTrial;
		
		%save data resp corr
		cellSaveConsistencyPerRespType{intPopulation,1} = [nanmean(vecHitConsistency) nanmean(vecMissConsistency) nanmean(vecSlowConsistency) nanmean(vecFastConsistency) nanmean(vecPreConsistency)];
		
		%get consistencies per stimulus trial, split by highest/lowest 50%
		vecConsistencyPerStimTrial = vecConsistencyPerTrial(indStimTrials);
		vecLowConsistencyTrials = vecConsistencyPerStimTrial < nanmedian(vecConsistencyPerStimTrial);
		vecHighConsistencyTrials = vecConsistencyPerStimTrial >= nanmedian(vecConsistencyPerStimTrial);
		
		%get decoding accuracy for lowest and highest 50% consistencytrials
		vecDecAccLowConsist = vecDecodedCorrect(vecLowConsistencyTrials);
		vecDecAccHighConsist = vecDecodedCorrect(vecHighConsistencyTrials);
		
		%get response types during stimulus trials
		indHitStim = indHit(indStimTrials);
		indMissStim = indMiss(indStimTrials);
		
		%save data
		cellSaveConsistencyPerTrial{intPopulation,1} = vecDecAccLowConsist;
		cellSaveConsistencyPerTrial{intPopulation,2} = vecDecAccHighConsist;
		cellSaveConsistencyPerTrial{intPopulation,3} = vecDecodedCorrect(indMissStim);
		cellSaveConsistencyPerTrial{intPopulation,4} = vecDecodedCorrect(indHitStim);
		
		%% check whether neurons close in time (in sequence) are more similarly tuned
		%get anatomical locations
		dblPix2Micron = cellMultiSes{intPopulation}.xml.dblActualImageSizeX / cellMultiSes{intPopulation}.xml.intImageSizeX;
		vecLocX = nan(1,intNeurons);
		vecLocY = nan(1,intNeurons);
		for intNeuron=1:intNeurons
			vecLocX(intNeuron) = cellMultiSes{intPopulation}.neuron(intNeuron).x;
			vecLocY(intNeuron) = cellMultiSes{intPopulation}.neuron(intNeuron).y;
		end
		matLocDists = dblPix2Micron*sqrt(bsxfun(@minus,vecLocX,vecLocX').^2 + bsxfun(@minus,vecLocY,vecLocY').^2);
		
		%sort
		[vecSorted,vecSort] = sort(vecMeanOffset);
		
		%get OSIs of neurons
		matSpikeCountsStim = matSpikeCounts(:,indStimTrials);
		sTuning = calcTuningRespMat(matSpikeCountsStim,cellSelectSubO,vecOrientations);
		vecNeuronPA = 2*ang2rad(sTuning.vecPrefAngle);
		mat_dPA = abs(circ_dist2(vecNeuronPA,vecNeuronPA));
		
		%get difference in temporal offset
		mat_dTO = abs(bsxfun(@minus,vecMeanOffset,vecMeanOffset'));
		
		%get lower half
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vec_dPA = mat_dPA(matSelect);
		vec_dTO = mat_dTO(matSelect);
		vec_dLoc = matLocDists(matSelect);
		
		%save figure
		if sParams.boolSavePlots
			figure
			subplot(2,2,1)
			imagesc(mat_dTO(vecSort,vecSort),[-0.4 0.4]);colormap(redblue);colorbar;freezeColors;cbfreeze;
			title('Difference in preferred temporal position (sd)');
			xlabel('Neuron');
			ylabel('Neuron');
			
			subplot(2,2,2)
			imagesc(matLocDists(vecSort,vecSort));colormap(hot);colorbar;freezeColors;cbfreeze;
			title('Difference in preferred orientation');
			xlabel('Neuron');
			ylabel('Neuron');
			
			
			subplot(2,2,3)
			scatter(vec_dLoc,vec_dTO,'.')
			%set(gca,'yscale','log')
			[r,p,low,hi]=corrcoef(vec_dLoc,vec_dTO);
			[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vec_dLoc,vec_dTO,0:25:250);
			errorbar(12.5:25:250,meanVec,stdVec./sqrt(nVec))
			ylabel('Pairwise latency difference');
			xlabel('Pairwise distance between somata (micron)');
			ylim([0 max(get(gca,'ylim'))])
			title(sprintf('Anatomical distance vs temporal position, r=%.3f,p=%.3f',r(2,1),p(2,1)));
			
			
			subplot(2,2,4)
			vecDiffRad = roundi(vec_dPA/pi,1);
			vecOriD0 = vecDiffRad == 0;
			vecOriD45 = vecDiffRad == 0.5;
			vecOriD90 = vecDiffRad == 1;
			
			boxplot(vec_dTO,vecDiffRad*90);
			ylabel('Pairwise latency difference');
			xlabel('Pairwise difference in preferred orientation');
			
			[h,p0_45]=ttest2(vec_dTO(vecOriD0),vec_dTO(vecOriD45));
			[h,p0_90]=ttest2(vec_dTO(vecOriD0),vec_dTO(vecOriD90));
			[h,p45_90]=ttest2(vec_dTO(vecOriD45),vec_dTO(vecOriD90));
			[dummy,dummy2,vecP_corr]=fdr_bh([p0_45 p0_90 p45_90]);
			title(sprintf('t-test, 0-45,p=%.3f; 0-90,p=%.3f; 90-45,p=%.3f',vecP_corr))
			drawnow;
			
			
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dLatencyDifference.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dLatencyDifference.pdf',strSes,intPopulation));
		end
		
		%% collect ordering of all neurons in each population event during miss/slow/fast trials, then check if ordering is more consistent during fast than slow trial assembly occurrences
		cellNeuronSequencesAssemblies = cell(intAssemblies,3);
		cellNeuronSequences = cellfill(nan(intNeurons,length(vecAssemblies)),[1 4]);
		intMissCounterOverall = 0;
		intSlowCounterOverall = 0;
		intFastCounterOverall = 0;
		intAllCounterOverall = 0;
		vecRunAssemblies = find(indMultiNeuronClusters)';%find(indRealAssemblies)';%1:intAssemblies
		for intAssembly=vecRunAssemblies
			indOccurrences = vecAssemblies == intAssembly;
			vecStarts = vecAssemblyStarts(indOccurrences);
			vecStops = vecAssemblyStops(indOccurrences);
			intOccurrences = numel(vecStarts);
			
			%pre-allocate these matrices
			cellNeuronSequencesAssemblies{intAssembly,1} = nan(intNeurons,intOccurrences);
			cellNeuronSequencesAssemblies{intAssembly,2} = nan(intNeurons,intOccurrences);
			cellNeuronSequencesAssemblies{intAssembly,3} = nan(intNeurons,intOccurrences);
			intMissCounterAssembly = 0;
			intSlowCounterAssembly = 0;
			intFastCounterAssembly = 0;
			
			%get sequences
			for intOccurrence=1:intOccurrences
				%get spike times for this occurrence
				intStart = vecStarts(intOccurrence);
				intStop = vecStops(intOccurrence)-1;
				matAEs = matSpikeCountsHD(:,intStart:intStop);
				
				%get mean and std
				matVals = repmat(1:size(matAEs,2),[size(matAEs,1) 1]);
				vecVals = getValByIdx(matVals,matAEs);
				dblMean = mean(vecVals);
				dblSD = std(vecVals);
				
				%check if during trial & get relative timing to center of mass of this occurrence
				if any(ismember(intStart:intStop,vecMissHD))
					intMissCounterAssembly = intMissCounterAssembly + 1;
					intMissCounterOverall = intMissCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,1}(:,intMissCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{1}(:,intMissCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecSlowHD))
					intSlowCounterAssembly = intSlowCounterAssembly + 1;
					intSlowCounterOverall = intSlowCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,2}(:,intSlowCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{2}(:,intSlowCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecFastHD))
					intFastCounterAssembly = intFastCounterAssembly + 1;
					intFastCounterOverall = intFastCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,3}(:,intFastCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{3}(:,intFastCounterOverall) = vecThisMuCC;
				end
				
				%get overall sd
				vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
				cellNeuronSequences{4}(:,intOccurrence) = vecThisMuCC;
				intAllCounterOverall = intAllCounterOverall + 1;
			end
			
			%remove pre-allocated nans
			if intMissCounterAssembly == 0,intMissCounterAssembly = 1;end
			if intSlowCounterAssembly == 0,intSlowCounterAssembly = 1;end
			if intFastCounterAssembly == 0,intFastCounterAssembly = 1;end
			cellNeuronSequencesAssemblies{intAssembly,1}(:,(intMissCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,2}(:,(intSlowCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,3}(:,(intFastCounterAssembly+1):end) = [];
		end
		%remove pre-allocated nans
		cellNeuronSequences{1}(:,(intMissCounterOverall+1):end) = [];
		cellNeuronSequences{2}(:,(intSlowCounterOverall+1):end) = [];
		cellNeuronSequences{3}(:,(intFastCounterOverall+1):end) = [];
		cellNeuronSequences{4}(:,(intAllCounterOverall+1):end) = [];
		
		%calculate mean+sd
		vecMeanLatencyMiss = nanmean(cellNeuronSequences{1},2);
		vecMeanLatencySlow = nanmean(cellNeuronSequences{2},2);
		vecMeanLatencyFast = nanmean(cellNeuronSequences{3},2);
		vecMeanLatencyOverall = nanmean(cellNeuronSequences{4},2);
		vecSDLatencyMiss = nanstd(cellNeuronSequences{1},[],2);
		vecSDLatencySlow = nanstd(cellNeuronSequences{2},[],2);
		vecSDLatencyFast = nanstd(cellNeuronSequences{3},[],2);
		vecSDLatencyOverall = nanstd(cellNeuronSequences{4},[],2);
		
		
		%select neurons with sufficient data points
		intMin = 5;
		indSelect1 = (sum(~isnan(cellNeuronSequences{1}),2)>intMin);
		indSelect2 = (sum(~isnan(cellNeuronSequences{2}),2)>intMin);
		indSelect3 = (sum(~isnan(cellNeuronSequences{3}),2)>intMin);
		indSelect4 = (sum(~isnan(cellNeuronSequences{4}),2)>intMin);
		
		vecMeanLatencyMiss(~indSelect1) = [];
		vecMeanLatencySlow(~indSelect2) = [];
		vecMeanLatencyFast(~indSelect3) = [];
		vecMeanLatencyOverall(~indSelect4) = [];
		vecSDLatencyMiss(~indSelect1) = [];
		vecSDLatencySlow(~indSelect2) = [];
		vecSDLatencyFast(~indSelect3) = [];
		vecSDLatencyOverall(~indSelect4) = [];
		
		%save data
		cellSaveTemporalSequenceStability{intPopulation} = [nanmean(vecSDLatencyMiss) nan nanmean(vecSDLatencySlow) nanmean(vecSDLatencyFast) nanmean(vecSDLatencyOverall)]; %correlation temporal order within assemblies with temporal order overall
		
		%{
		%plot
		figure
		subplot(2,2,1)
		scatter(vecMeanLatencyMiss,vecSDLatencyMiss,'r')
		
		subplot(2,2,2)
		scatter(vecMeanLatencySlow,vecSDLatencySlow,'m')
		
		subplot(2,2,3)
		scatter(vecMeanLatencyFast,vecSDLatencyFast,'g')
		
		subplot(2,2,4)
		scatter(vecMeanLatencyFast,vecSDLatencyFast,'k')
		%}
		if sParams.boolSavePlots
			%plot
			hSequenceStability = figure;
			
			matLatenciesMiss = cellNeuronSequences{1};
			vecLatenciesMiss = nanmean(matLatenciesMiss,2);
			matLatenciesHit = [cellNeuronSequences{2} cellNeuronSequences{3}];
			vecLatenciesHit = nanmean(matLatenciesHit,2);
			vecLatencies = mean(cat(2,vecLatenciesMiss,vecLatenciesHit),2);
			[dummy,vecReorder] = sort(vecLatencies);
			
			subplot(2,2,1);
			
			vecLimMiss = 2*(nanmean(matLatenciesMiss(:))+nanstd(matLatenciesMiss(:)))*[-1 1];
			matColormap = colormap('redbluepurple');
			implot(matLatenciesMiss(vecReorder,:),vecLimMiss,matColormap,3);
			%colorbar;
			title('Neuronal response latency relative to pop. CoM [Miss]');
			ylabel('Neuron')
			xlabel('Population event during miss trial')
			
			subplot(2,2,2);
			vecLimHit = 2*(nanmean(matLatenciesHit(:))+nanstd(matLatenciesHit(:)))*[-1 1];
			matColormap = colormap('redbluepurple');
			implot(matLatenciesHit(vecReorder,:),vecLimHit,matColormap,3);
			%colorbar;
			title('Neuronal response latency relative to pop. CoM [Hit]');
			ylabel('Neuron')
			xlabel('Population event during hit trial')
			
			
			subplot(2,2,3);
			errorbar(1,nanmean(vecSDLatencyMiss),nanstd(vecSDLatencyMiss)/sqrt(sum(~isnan(vecSDLatencyMiss))),'xr');
			hold on;
			errorbar(2,nanmean(vecSDLatencySlow),nanstd(vecSDLatencySlow)/sqrt(sum(~isnan(vecSDLatencySlow))),'xm');
			errorbar(3,nanmean(vecSDLatencyFast),nanstd(vecSDLatencyFast)/sqrt(sum(~isnan(vecSDLatencyFast))),'xg');
			hold off;
			%ylim([0 max(get(gca,'ylim'))]);
			set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Slow','Fast'});
			ylabel('Temporal sequence variability (s)');
			title('Mean consistency of activation order in population events across neurons');
			
			subplot(2,2,4);
			errorbar(1,nanmean(vecMeanLatencyMiss),nanstd(vecMeanLatencyMiss)/sqrt(sum(~isnan(vecMeanLatencyMiss))),'xr');
			hold on;
			errorbar(2,nanmean(vecMeanLatencySlow),nanstd(vecMeanLatencySlow)/sqrt(sum(~isnan(vecMeanLatencySlow))),'xm');
			errorbar(3,nanmean(vecMeanLatencyFast),nanstd(vecMeanLatencyFast)/sqrt(sum(~isnan(vecMeanLatencyFast))),'xg');
			hold off;
			%ylim([0 max(get(gca,'ylim'))]);
			set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Slow','Fast'});
			ylabel('Temporal sequence variability (s)');
			title('Mean consistency of activation order in population events across neurons');
			
			%save figure
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dtemporal_sequence_stability.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dtemporal_sequence_stability.pdf',strSes,intPopulation));
		end
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnalFour' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end