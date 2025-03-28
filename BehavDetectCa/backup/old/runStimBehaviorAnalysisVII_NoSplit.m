%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
for intMouse=5
	close all
	clearvars -except intMouse
	boolUseNeuropilSubtraction = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == -1
			%return; %exclude?; bad behavior, weird signals
			strSes = '20140430';
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
		load(['C:\Data\Processed\StimDetectionAgg\dataPreProAggregate' strSes '.mat']);
	end
	
	%load separate ses files
	for intFile=1:numel(cellAggregate)
		%load data
		sLoad = load([cellAggregate{intFile} '.mat']);
		ses = sLoad.ses;
		
		%transform orientation to
		ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
		
		%assign to multi-cell
		cellSes{intFile} = ses;
	end
	
	%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	vecBlockTypes = unique(vecBlock);
	intNumBlocks = length(vecBlockTypes);
	%vecNeuronNum = zeros(1,intNumBlocks);
	%cellKeepList = cell(1,intNumBlocks);
	%#ok<*ASGLU>
	%#ok<*AGROW>
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['C:\Data\ResultsStimDetectionCa' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	strSes = ['NS' strSes];
	
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		%take opposite directions as the same
		cellMultiSes{intPopulation}.structStim.Orientation = mod(cellMultiSes{intPopulation}.structStim.Orientation,180);
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		intContrasts = length(cellSelectContrasts);
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
		cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		%calculate preferred stimulus orientation for all neurons
		vecOriPrefContrasts = 4:6; % 8% - 100%
		matPref = nan(length(vecOriPrefContrasts),intNeurons);
		intCounter = 0;
		for intContrast = vecOriPrefContrasts
			intCounter = intCounter + 1;
			matResp = matTrialNormResponse(:,cellSelectContrasts{intContrast});
			structStimC{intContrast} = cellMultiSes{intPopulation}.structStim;
			cellFields = fieldnames(cellMultiSes{intPopulation}.structStim);
			for intField=1:length(cellFields)
				strField = cellFields{intField};
				structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
			end
			cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
			sTuning{intContrast} = calcTuningRespMat(matResp,cellSelect,vecOrientations);
			matPref(intCounter,:) = sTuning{intContrast}.vecPrefIndex;
		end
		vecNeuronPrefStim = nan(1,intNeurons);
		for intOri=1:length(vecOrientations);
			vecNeuronPrefStim(sum(matPref == intOri,1) > 1) = intOri;
		end
		%vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = randi(length(vecOrientations),[1 sum(isnan(vecNeuronPrefStim))]);
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(isnan(vecNeuronPrefStim)) = [];
		vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		intOris = length(vecOrientations);
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		intContrasts = length(cellSelectContrasts);
		
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Contrast,unique(cellMultiSes{intPopulation}.structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn<3;
		
		
		%{
		%make resp mat
		matRespTime = nan(intNeurons,length(cellMultiSes{intPopulation}.neuron(intNeuron).dFoF));
		vecGauss = normpdf(-2:1:2,0,1);
		for intNeuron=1:intNeurons
			matRespTime(intNeuron,:) = conv(cellMultiSes{intPopulation}.neuron(intNeuron).dFoF,vecGauss,'same');
		end
		
		%%
		
		h=figure;
		for intTrial=1:intTrials
			intStimOn = cellMultiSes{intPopulation}.structStim.FrameOn(intTrial);
			intStimOff = cellMultiSes{intPopulation}.structStim.FrameOff(intTrial);
			intStart = intStimOn-round(3*cellMultiSes{intPopulation}.samplingFreq);
			intStop = intStimOff+round(60*cellMultiSes{intPopulation}.samplingFreq);
			
			clf;
			hold on
			vecStims = (cellMultiSes{intPopulation}.structStim.FrameOn < intStop & cellMultiSes{intPopulation}.structStim.FrameOn > intStart) | (cellMultiSes{intPopulation}.structStim.FrameOff < intStop & cellMultiSes{intPopulation}.structStim.FrameOff > intStart);
			vecLimY = [-1 0];
			for intStim=find(vecStims);
				dblStart = cellMultiSes{intPopulation}.structStim.FrameOn(intStim);
				dblStop = cellMultiSes{intPopulation}.structStim.FrameOff(intStim);
				
				fill([dblStart dblStop dblStop dblStart]/cellMultiSes{intPopulation}.samplingFreq,[0 0 -1 -1],[0.5 0.5 0.5]);
			end
			
			hold off
			
			hold all
			dblMinValCounter=0;
			vecNeurons = [3 9 66 randperm(intNeurons,7)];
			for intNeuron=vecNeurons
				vecTrace = matRespTime(intNeuron,intStart:intStop);
				vecTrace=vecTrace-min(vecTrace);
				dblRange = range(vecTrace);
				dblMinValCounter=dblMinValCounter+dblRange;
				
				plot((intStart:intStop)/cellMultiSes{intPopulation}.samplingFreq,vecTrace-dblMinValCounter);
			end
			hold off
			legend(cellfun(@num2str,num2cell(vecNeurons),'UniformOutput',false))
			
			drawnow;
			pause;
		end
		close(h);
		
		%}
		
		
		
		%% calculate relative increase per neuron from miss => hit trials
		%mean over all test contrasts
		matMeanMissAct = nan(intContrasts,intOris,intNeurons);
		matMeanHitAct = nan(intContrasts,intOris,intNeurons);
		matMeanMissSD = nan(intContrasts,intOris,intNeurons);
		matMeanHitSD = nan(intContrasts,intOris,intNeurons);
		for intContrast=1:intContrasts
			indSelectContrastTrials = cellSelectContrasts{intContrast};
			for intOri=1:intOris
				indSelectOriTrials = cellSelectOri{intOri};
				matMeanMissAct(intContrast,intOri,:) = mean(matTrialResponse(:,indSelectContrastTrials & indSelectOriTrials & ~indStimResp),2);
				matMeanHitAct(intContrast,intOri,:) = mean(matTrialResponse(:,indSelectContrastTrials & indSelectOriTrials & indStimResp),2);
				matMeanMissSD(intContrast,intOri,:) = std(matTrialResponse(:,indSelectContrastTrials & indSelectOriTrials & ~indStimResp),[],2);
				matMeanHitSD(intContrast,intOri,:) = std(matTrialResponse(:,indSelectContrastTrials & indSelectOriTrials & indStimResp),[],2);
			end
		end
		
		
		%get activity per trial relative to mean miss act for that type
		matTrialRelAct = nan(size(matTrialResponse));
		for intTrial=1:intTrials
			vecMissAct = squeeze(matMeanMissAct(vecStimContrasts(intTrial),vecStimOris(intTrial),:));
			vecMissSD = squeeze(matMeanMissSD(vecStimContrasts(intTrial),vecStimOris(intTrial),:));
			vecThisAct = matTrialResponse(:,intTrial);
			matTrialRelAct(:,intTrial) = (vecThisAct-vecMissAct)./vecMissSD;
		end
		
		matTrialRelActHit = matTrialRelAct(:,indStimResp&(~cellSelectContrasts{1} & ~cellSelectContrasts{end}));
		
		%calculate inter-trial correlation of relative hit-increase
		matCorr=corr(matTrialRelAct);
		
		%split per trial type
		cellRelAct = cell(intContrasts,intOris);
		cellRelHitIncCorrs = cell(intContrasts,intOris);
		cellRelMissIncCorrs = cell(intContrasts,intOris);
		for intContrast=1:intContrasts
			indSelectContrastTrials = cellSelectContrasts{intContrast};
			for intOri=1:intOris
				indSelectOriTrials = cellSelectOri{intOri};
				matSelect = tril(true(size(matCorr)),-1);
				matSelect2 = false(size(matSelect));
				matSelect2(indSelectContrastTrials&indSelectOriTrials&indStimResp,indSelectContrastTrials&indSelectOriTrials&indStimResp) = true;
				matSelect3 = false(size(matSelect));
				matSelect3(indSelectContrastTrials&indSelectOriTrials&~indStimResp,indSelectContrastTrials&indSelectOriTrials&~indStimResp) = true;
				cellRelHitIncCorrs{intContrast,intOri} = matCorr(matSelect&matSelect2);
				cellRelMissIncCorrs{intContrast,intOri} = matCorr(matSelect&matSelect3);
			end
		end
		
		%% NEW FIG
		figure
		
		%plot relative hit enhancement per trial
		subplot(2,2,1)
		%imagesc(matTrialRelActHit,[-1 1]*max(abs(matTrialRelActHit(:))));
		vecRange = [-5 5];
		intHitTrials = size(matTrialRelActHit,2);
		imagesc(matTrialRelActHit,vecRange);
		colormap('redblue');
		colorbar;
		title('Relative hit enhancement per trial')
		ylabel('Neuron ID')
		xlabel('Hit trial number')
		
		%plot prediction by neurons
		subplot(2,2,2)
		matNeuronPredictor = repmat(mean(matTrialRelActHit,2),[1 intHitTrials]);
		imagesc(matNeuronPredictor,vecRange);
		colormap('redblue');
		colorbar;
		
		dblVarRes=sum((matTrialRelActHit(:)-matNeuronPredictor(:)).^2);
		dblVarTot=sum((matTrialRelActHit(:)).^2);
		dblPercExplainedNeurons = 1-dblVarRes/dblVarTot;
		dblPercExplainedNeuronsAdj = dblPercExplainedNeurons*(intHitTrials/intNeurons);
		title(sprintf('Predicted by neuron ID (R^2=%.3f; R^2-adj=%.3f)',dblPercExplainedNeurons,dblPercExplainedNeuronsAdj))
		ylabel('Neuron ID')
		xlabel('Hit trial number')
		
		%plot prediction by trials
		subplot(2,2,3)
		matTrialPredictor = repmat(mean(matTrialRelActHit,1),[size(matTrialRelActHit,1) 1]);
		imagesc(matTrialPredictor,vecRange);
		colormap('redblue');
		colorbar;
		
		dblVarRes=sum((matTrialRelActHit(:)-matTrialPredictor(:)).^2);
		dblVarTot=sum((matTrialRelActHit(:)).^2);
		dblPercExplainedTrials = 1-dblVarRes/dblVarTot;
		dblPercExplainedTrialsAdj = dblPercExplainedTrials*(intNeurons/intHitTrials);
		
		%dblPercExplained = corr(matTrialPredictor(:),matTrialRelActHit(:))*100;
		
		title(sprintf('Predicted by trial ID (R^2=%.3f; R^2-adj=%.3f)',dblPercExplainedTrials,dblPercExplainedTrialsAdj))
		ylabel('Neuron ID')
		xlabel('Hit trial number')
		
		%plot dstirubiton
		subplot(2,2,4)
		matAllPredictor = matTrialPredictor+matNeuronPredictor;
		imagesc(matAllPredictor,vecRange);
		colormap('redblue');
		colorbar;
		
		dblVarRes=sum((matTrialRelActHit(:)-matAllPredictor(:)).^2);
		dblVarTot=sum((matTrialRelActHit(:)).^2);
		dblPercExplainedAll = 1-dblVarRes/dblVarTot;
		
		title(sprintf('Predicted by both trial and neuron ID (R^2=%.3f)',dblPercExplainedAll))
		ylabel('Neuron ID')
		xlabel('Hit trial number')
		
		
		%% shuffle analysis
		%neurons
		intIters = 1000;
		vecSigNr = nan(1,intIters);
		vecNeuronR2 = nan(1,intIters);
		for intIter=1:intIters
			matShuffledNeurons = matTrialRelActHit;
			for intTrial=1:intHitTrials
				matShuffledNeurons(:,intTrial) = matTrialRelActHit(randperm(intNeurons),intTrial);
			end
			%get nr of sign neurons
			[h,p]=ttest(matShuffledNeurons,0,0.05,'both',2);
			vecSigNr(intIter) = sum(h);
			
			%get explained variance
			matNeuronPredictor = repmat(mean(matShuffledNeurons,2),[1 intHitTrials]);
			dblVarRes=sum((matShuffledNeurons(:)-matNeuronPredictor(:)).^2);
			dblVarTot=sum((matShuffledNeurons(:)).^2);
			vecNeuronR2(intIter) = 1-dblVarRes/dblVarTot;
		end
		
		%trials
		intIters = 1000;
		vecTrialR2 = nan(1,intIters);
		for intIter=1:intIters
			matShuffledTrials = matTrialRelActHit;
			for intNeuron=1:intNeurons
				matShuffledTrials(intNeuron,:) = matTrialRelActHit(intNeurons,randperm(intHitTrials));
			end

			%get explained variance
			matTrialPredictor = repmat(mean(matShuffledTrials,1),[intNeurons 1]);
			dblVarRes=sum((matShuffledTrials(:)-matTrialPredictor(:)).^2);
			dblVarTot=sum((matShuffledTrials(:)).^2);
			vecTrialR2(intIter) = 1-dblVarRes/dblVarTot;
		end
		
		%% NEW FIG
		figure
		[h,p]=ttest(matTrialRelActHit,0,0.05,'both',2);
		intNrSig = sum(h);
		
		%CI explained variance neurons
		subplot(2,2,1)
		histx(vecNeuronR2);
		xlim([0 roundi(dblPercExplainedNeurons,2,'ceil')]);
		hold on
		plot([dblPercExplainedNeurons dblPercExplainedNeurons],get(gca,'ylim'),'r')
		hold off
		dblNrSD = (dblPercExplainedNeurons-mean(vecNeuronR2))/std(vecNeuronR2);
		title(sprintf('R^2 by neuron ID when shuffled [real=%.3f; %.1fsd from mean shuffled]',dblPercExplainedNeurons,dblNrSD))
		xlabel('Explained variance')
		ylabel('Number of iterations (count)')
		
		%plot dstirubiton
		subplot(2,2,2)
		histx(mean(matTrialRelActHit,2));
		title(sprintf('Distribution neurons hit dF/F0 relative to miss,mean=%.3f',mean(mean(matTrialRelActHit,2))))
		xlabel('Relative enhancement hit trial activity (sd)')
		ylabel('Number of neurons (count)')
		
		%CI nr of significant neurons
		subplot(2,2,3)
		intMaxS = max(vecSigNr);
		vecCounts = 0:intMaxS;
		bar(vecCounts,sum(bsxfun(@eq,vecSigNr,vecCounts'),2));
		xlim([0 max([intMaxS intNrSig+1])]);
		hold on
		plot([intNrSig intNrSig],get(gca,'ylim'),'r')
		hold off
		dblNrSD = (intNrSig-mean(vecSigNr))/std(vecSigNr);
		title(sprintf('Consistently hit-modulated neurons when shuffled [real=%d;%.1fsd from mean]',intNrSig,dblNrSD))
		xlabel('Number of significantly hit-modulated neurons')
		ylabel('Number of iterations (count)')
		set(gca,'XTickMode','auto')
		
		%CI nr of significant neurons
		subplot(2,6,10)
		errorbar(0.5,mean(vecNeuronR2),std(vecNeuronR2)*2,'bx');
		hold on
		scatter(0.5,dblPercExplainedNeurons,'rx')
		hold off
		set(gca,'xtick',0.5,'xticklabel','R^2 neurons')
		ylabel('Explained variance by neuron ID vs shuffled (+/- 2sd)')
		
		subplot(2,6,11)
		[vecSigNeurons,col]=find(bsxfun(@eq,vecSigNr,vecCounts'));
		errorbar(0.5,mean(vecSigNeurons),std(vecSigNeurons)*2,'bx');
		hold on
		scatter(0.5,intNrSig,'rx')
		hold off
		set(gca,'xtick',0.5,'xticklabel','Nr sig. neurons')
		ylabel('Nr of sign. modulated neurons vs shuffled (+/- 2sd)')
		
		subplot(2,6,12)
		errorbar(0.5,mean(vecTrialR2),std(vecTrialR2)*2,'bx');
		hold on
		scatter(0.5,dblPercExplainedTrials,'rx')
		hold off
		set(gca,'xtick',0.5,'xticklabel','R^2 trials')
		ylabel('Explained variance by trial ID vs shuffled (+/- 2sd)')
		
		return
		
		
		%% make sigmoid
		%get updated response matrices
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		
		%make sigmoid plots
		cellHitActPrefPop = cell(1,intContrasts);
		cellMissActPrefPop = cell(1,intContrasts);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			cellHitActPrefPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissActPrefPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		
		
		%loop through trials
		for intTrial = 1:intTrials
			intContrast = vecStimContrasts(intTrial);
			if intContrast == 0
				indNeurons = true(size(vecNeuronPrefStim));
			else
				indNeurons = vecNeuronPrefStim==vecStimOris(intTrial);
			end
			
			vecAct = matTrialResponse(indNeurons,intTrial);
			
			if indStimResp(intTrial)
				cellHitActPrefPop{intContrast}(find(isnan(cellHitActPrefPop{intContrast}),1)) = mean(vecAct);
			else
				cellMissActPrefPop{intContrast}(find(isnan(cellMissActPrefPop{intContrast}),1)) = mean(vecAct);
			end
		end
		
		%remove trailing nans
		for intContrast=1:length(cellHitActPrefPop)
			cellHitActPrefPop{intContrast} = cellHitActPrefPop{intContrast}(~isnan(cellHitActPrefPop{intContrast}));
			cellMissActPrefPop{intContrast} = cellMissActPrefPop{intContrast}(~isnan(cellMissActPrefPop{intContrast}));
		end
		vecMeanHit = cellfun(@mean,cellHitActPrefPop);
		vecSEHit = cellfun(@std,cellHitActPrefPop)./cellfun(@numel,cellHitActPrefPop);
		vecMeanMiss = cellfun(@mean,cellMissActPrefPop);
		vecSEMiss = cellfun(@std,cellMissActPrefPop)./cellfun(@numel,cellMissActPrefPop);
		
		%plot
		h=figure;
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		errorfill(vecPlotC,vecMeanMiss,vecSEMiss,[1 0 0],[1 0.5 0.5]);
		hold on
		errorfill(vecPlotC,vecMeanHit,vecSEHit,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean pref pop dF/F0')
		xlabel('Stimulus contrast (%)')
		
		%%
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_dFoF_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% create hit-correlated activity region plot
		%{
	loop through contrasts, select preferred stimulus trials per neuron,
	calculate df/f activity for pref stim response for hits and for misses,
	calculate hit/miss difference per neuron per contrast per trial and
	plot binned matrices
	
	calculate correlation between 20% most hit-correlated activity enhanced
	neurons for different contrasts, make plot
	
	calculate overall HCAE for all neurons; vecStimDetectActInc
	select neurons with overall 20% highest HCAE as mean over different
	contrasts; indSelectHitCorrelatedNeurons and 20% lowest HCAE;
	indSelectHitAnticorrelatedNeurons
		%}

		
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		cellHitPrefAct = cell(intContrasts,1);
		cellMissPrefAct = cell(intContrasts,1);
		cellHitNonPrefAct = cell(intContrasts,1);
		cellMissNonPrefAct = cell(intContrasts,1);
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
			cellHitPrefAct{intContrastIndex} = nan(1,intNeurons*(intTrials/intContrasts));
			cellMissPrefAct{intContrastIndex} = nan(1,intNeurons*(intTrials/intContrasts));
			cellHitNonPrefAct{intContrastIndex} = nan(1,intNeurons*(intTrials/intContrasts));
			cellMissNonPrefAct{intContrastIndex} = nan(1,intNeurons*(intTrials/intContrasts));
		end
		
		%split trial responses for pref pop / non-pref pop
		intCountHP = 0;
		intCountHNP = 0;
		intCountMP = 0;
		intCountMNP = 0;
		for intTrial=1:intTrials
			indPrefNeurons = vecNeuronPrefStim==vecStimOris(intTrial);
			vecPrefAct = matRespNormPerContrast(indPrefNeurons,intTrial);
			vecNonPrefAct = matRespNormPerContrast(indPrefNeurons,intTrial);
			intContrastIndex = vecStimContrasts(intTrial);
			
			if indStimResp(intTrial)
				cellHitPrefAct{intContrastIndex}((intCountHP+1):intCountHP+length(vecPrefAct)) = vecPrefAct;
				cellHitNonPrefAct{intContrastIndex}((intCountHNP+1):intCountHNP+length(vecNonPrefAct)) = vecNonPrefAct;
				intCountHP = intCountHP + length(vecPrefAct);
				intCountHNP = intCountHNP + length(vecNonPrefAct);
			else
				cellMissPrefAct{intContrastIndex}((intCountMP+1):intCountMP+length(vecPrefAct)) = vecPrefAct;
				cellMissNonPrefAct{intContrastIndex}((intCountMNP+1):intCountMNP+length(vecNonPrefAct)) = vecNonPrefAct;
				intCountMP = intCountMP + length(vecPrefAct);
				intCountMNP = intCountMNP + length(vecNonPrefAct);
			end
		end
		
		
		%remove trailing nans
		for intContrastIndex=1:length(cellSelectContrasts)
			cellHitPrefAct{intContrastIndex} = cellHitPrefAct{intContrastIndex}(~isnan(cellHitPrefAct{intContrastIndex}));
			cellMissPrefAct{intContrastIndex} = cellMissPrefAct{intContrastIndex}(~isnan(cellMissPrefAct{intContrastIndex}));
			cellHitNonPrefAct{intContrastIndex} = cellHitNonPrefAct{intContrastIndex}(~isnan(cellHitNonPrefAct{intContrastIndex}));
			cellMissNonPrefAct{intContrastIndex} = cellMissNonPrefAct{intContrastIndex}(~isnan(cellMissNonPrefAct{intContrastIndex}));
		end
		
		%data
		intSwitchZ = 1;
		if intSwitchZ~=1
			vecBins = -0.05:0.01:0.15;
			vecTickY = [6 16];
			matHit = zeros(length(vecBins),intContrasts)+eps;
			matMiss = zeros(length(vecBins),intContrasts)+eps;
			for intRepRec = 1:intRepRecs
				for intContrast=1:intContrasts
					matHit(:,intContrast) = matHit(:,intContrast) + hist(cellHitPrefResp{intContrast,intRepRec},vecBins)';
					matMiss(:,intContrast) = matMiss(:,intContrast) + hist(cellMissPrefResp{intContrast,intRepRec},vecBins)';
				end
			end
			strLabelY = sprintf('dF/F0 distribution');
		else
			vecBins = -1:0.1:2;
			vecTickY = [1 11 21 31];
			matHit = zeros(length(vecBins),intContrasts);
			matMiss = zeros(length(vecBins),intContrasts);
			for intRepRec = 1:intRepRecs
				for intContrast=1:intContrasts
					matHit(:,intContrast) = matHit(:,intContrast) + hist(zscore(cellHitPrefResp{intContrast,intRepRec}),vecBins)';
					matMiss(:,intContrast) = matMiss(:,intContrast) + hist(zscore(cellMissPrefResp{intContrast,intRepRec}),vecBins)';
				end
			end
			strLabelY = sprintf('Z-scored dF/F0 distribution');
		end
		%if intMouse == 8 %only 2 hit resps at 100%, so remove from analysis
		%	matMiss(:,end) = mean(matHit(:,end));
		%end
		
		
		
		%normalize
		%%{
		matHit = buildConvolutedMatrix(matHit,1,(1/3)*[1;1;1]);
		matMiss = buildConvolutedMatrix(matMiss,1,(1/3)*[1;1;1]);
		matHit = matHit(3:(end-2),:);
		matMiss = matMiss(3:(end-2),:);
		matHit = matHit./repmat(sum(matHit,1),[size(matHit,1) 1]);
		matMiss = matMiss./repmat(sum(matMiss,1),[size(matMiss,1) 1]);
		matDiff = (matHit ./ (matMiss + matHit))-0.5;
		matFilt = normpdf(-1:1,0,1)'*normpdf(-1:1,0,1);
		matFilt=matFilt/sum(matFilt(:));
		matDiff = conv2(matDiff,matFilt,'same');
		cellSaveMatrices{1,intPopulation} = matHit;
		cellSaveMatrices{2,intPopulation} = matMiss;
		cellSaveMatrices{3,intPopulation} = matDiff;
		%%}
		
		%vars
		vecContrasts = [0 0.5 2 8 32 100];
		vecTickLabelY = vecBins(vecTickY);
		vecTickLabelX = vecContrasts;
		
		%plot over contrasts
		hActMat = figure;
		set(hActMat,'Color',[1 1 1]);
		figure(hActMat);
		
		%hits
		subplot(2,2,1)
		imagesc(matHit);
		axis xy;
		title('Pref / Hit; z=normalized count')
		set(gca,'YTick',vecTickY-2,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%misses
		subplot(2,2,2)
		imagesc(matMiss);
		axis xy;
		title('Pref / Miss; z=normalized count')
		set(gca,'YTick',vecTickY-2,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%difference
		subplot(2,2,3)
		imagesc(matDiff,max(abs(matDiff(:)))*[-1 1]);
		axis xy;
		title('Diff')
		set(gca,'YTick',vecTickY-2,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(redblue(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%select which bins to plot
		intTestType = 3;
		intBins = length(matDiff(:,1));
		matH = false(size(matDiff));
		matP = nan(size(matDiff));
		if intTestType < 3
			for intContrastIndex=1:size(matDiff,2)
				if intTestType == 1
					%bonferroni corrected ttest
					alpha=0.05/intBins;
					for intBin=1:length(matDiff(:,intContrastIndex))
						[h,p]=ttest(matDiff(intBin,intContrastIndex)-0.5,alpha);
						matP(intBin,intContrastIndex) = p;
						matH(intBin,intContrastIndex) = h;
					end
				elseif intTestType == 2
					%z-score bin count percentage, then take  >2sd's
					matP(:,intContrastIndex) = zscore(matDiff(:,intContrastIndex));
					matH(:,intContrastIndex) = matP(:,intContrastIndex)>1.5;
				end
			end
		elseif intTestType == 3
			matDiffT = matDiff;
			matDiffT(isnan(matDiffT)) = nanmean(matDiff(:));
			matP = reshape(zscore(matDiffT(:)),size(matDiff));
			matH = matP>2;
		end
		
		%select for calculation all bins that are significant on any contrast
		matPlotP = ones(size(matDiff));
		matPlotP(matH) = matP(matH);
		subplot(2,2,4)
		imagesc(matH)
		cellSaveMatrices{4,intPopulation} = matH;
		axis xy
		title('Significant outlier bins')
		set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(2));colorbar;drawnow;freezeColors;cbfreeze;
		%%
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_heatmaps_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% calculate noise + signal correlations FOR HIT VERSUS MISS
		%get miss/hit/slow
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		
		structStimMiss = remel(cellMultiSes{intPopulation}.structStim,indMiss);
		structStimHit = remel(cellMultiSes{intPopulation}.structStim,indHit);
		structStimFast = remel(cellMultiSes{intPopulation}.structStim,indFast);
		structStimSlow = remel(cellMultiSes{intPopulation}.structStim,indSlow);
		
		sesMiss = cellMultiSes{intPopulation};
		sesMiss.structStim = structStimMiss;
		sesHit = cellMultiSes{intPopulation};
		sesHit.structStim = structStimHit;
		sesFast = cellMultiSes{intPopulation};
		sesFast.structStim = structStimFast;
		sesSlow = cellMultiSes{intPopulation};
		sesSlow.structStim = structStimSlow;
		
		%calc correlations
		sCorrMiss = calcStimCorrsAsym(sesMiss,{'Orientation','Contrast'});
		matSC_Miss = sCorrMiss.matSignalCorrs;
		matNC_Miss = sCorrMiss.matNoiseCorrs;
		sCorrHit = calcStimCorrsAsym(sesHit,{'Orientation','Contrast'});
		matSC_Hit = sCorrHit.matSignalCorrs;
		matNC_Hit = sCorrHit.matNoiseCorrs;
		sCorrFast = calcStimCorrsAsym(sesFast,{'Orientation','Contrast'});
		matSC_Fast = sCorrFast.matSignalCorrs;
		matNC_Fast = sCorrFast.matNoiseCorrs;
		sCorrSlow = calcStimCorrsAsym(sesSlow,{'Orientation','Contrast'});
		matSC_Slow = sCorrSlow.matSignalCorrs;
		matNC_Slow = sCorrSlow.matNoiseCorrs;
		matSelect = tril(true(size(matSC_Miss)),-1);
		
		%signal correlations
		hCorr = figure;
		colormap(redblue(128))
		subplot(4,4,1)
		imagesc(matSC_Miss,[-1 1])
		title('SC Miss')
		subplot(4,4,2)
		imagesc(matSC_Hit,[-1 1])
		title('SC Hit')
		subplot(4,4,5)
		imagesc(matSC_Fast,[-1 1])
		title('SC Fast')
		subplot(4,4,6)
		imagesc(matSC_Slow,[-1 1])
		title('SC Slow')
		subplot(2,2,3)
		vecSC_M = matSC_Miss(matSelect);
		vecSC_H = matSC_Hit(matSelect);
		vecSC_F = matSC_Fast(matSelect);
		vecSC_S = matSC_Slow(matSelect);
		errorbar(1:4,[mean(vecSC_M) mean(vecSC_H) mean(vecSC_F) mean(vecSC_S)],[std(vecSC_M) std(vecSC_H) std(vecSC_F) std(vecSC_S)]/sqrt(sum(matSelect(:))),'Linestyle','none','Marker','x','Color','k');
		set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
		title('Signal correlations')
		
		%noise correlations
		subplot(4,4,3)
		imagesc(matNC_Miss,[-1 1])
		title('NC Miss')
		subplot(4,4,4)
		imagesc(matNC_Hit,[-1 1])
		title('NC Hit')
		subplot(4,4,7)
		imagesc(matNC_Fast,[-1 1])
		title('NC Fast')
		subplot(4,4,8)
		imagesc(matNC_Slow,[-1 1])
		title('NC Slow')
		subplot(2,2,4)
		vecNC_M = matNC_Miss(matSelect);
		vecNC_H = matNC_Hit(matSelect);
		vecNC_F = matNC_Fast(matSelect);
		vecNC_S = matNC_Slow(matSelect);
		errorbar(1:4,[mean(vecNC_M) mean(vecNC_H) mean(vecNC_F) mean(vecNC_S)],[std(vecNC_M) std(vecNC_H) std(vecNC_F) std(vecNC_S)]/sqrt(sum(matSelect(:))),'Linestyle','none','Marker','x','Color','k');
		set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
		title('Noise correlations')
		
		%put in output
		cellSaveStimCorrs{1,intPopulation} = matSC_Miss;
		cellSaveStimCorrs{2,intPopulation} = matNC_Miss;
		cellSaveStimCorrs{3,intPopulation} = matSC_Hit;
		cellSaveStimCorrs{4,intPopulation} = matNC_Hit;
		cellSaveStimCorrs{5,intPopulation} = matSC_Fast;
		cellSaveStimCorrs{6,intPopulation} = matNC_Fast;
		cellSaveStimCorrs{7,intPopulation} = matSC_Slow;
		cellSaveStimCorrs{8,intPopulation} = matNC_Slow;
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_HCAEneurons_correlations_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		%% heterogeneity for hit/miss
		%group per pref stim
		vecReorder = [];
		vecReorderedPrefStim = [];
		for intStim=unique(vecNeuronPrefStim)
			%neurons
			vecSubPop = find(vecNeuronPrefStim==intStim);
			vecReorder = [vecReorder vecSubPop];
			vecReorderedPrefStim = [vecReorderedPrefStim ones(1,length(vecSubPop))*intStim];
		end
		
		%loop through contrasts & stim types for trial reordering
		vecTrialOriVector = [];
		vecTrialContrastVector = [];
		
		for intStim=unique(vecNeuronPrefStim)
			for intContrastIndex=1:length(cellSelectContrasts)
				%trials
				indSelectOri = cellSelectOri{intStim};
				vecTrialOriVector = [vecTrialOriVector find(indSelectOri)];
				
				vecSelectContrastTrials = cellSelectContrasts{intContrastIndex} & indSelectOri;
				vecTrialContrastVector = [vecTrialContrastVector find(vecSelectContrastTrials)];
			end
		end
		
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		%plot
		figure
		colormap(hot(256))
		%imagesc(matRespNormPerContrast)
		imagesc(matRespNormPerContrast(vecReorder,vecTrialContrastVector))
		set(gca,'XTick',1:sum(vecSelectContrastTrials):length(vecTrialContrastVector))
		colorbar
		title('Z-scored activation level per neuron per trial')
		ylabel('Neuron number')
		xlabel('Trial number, sorted by stimulus ori&contrast')
		%save plot
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%high/low DCAE
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intValuesTotal = (intNeurons*intNeurons-intNeurons)/2;
		matHeterogeneity = nan(intTrials,intValuesTotal);
		matRespDistNoRem = nan(intTrials,intValuesTotal);
		matRespDistLoZRem = nan(intTrials,intValuesTotal);
		matRespDistHiZRem = nan(intTrials,intValuesTotal);
		matRespDistLoARem = nan(intTrials,intValuesTotal);
		matRespDistHiARem = nan(intTrials,intValuesTotal);
		
		matTrialIndex = nan(intTrials,intValuesTotal);
		%calculate inverse correlation matrix per trial
		for intTrial=1:intTrials
			%general
			matTrialIndex(intTrial,:) = intTrial;
			
			%% STANDARD
			%perform calculation for all neurons
			vecActivity = matRespNormPerContrast(:,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeurons]);
			matZ2 = repmat(vecActivity',[intNeurons 1]);
			matDistAll = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectAll = tril(true(size(matDistAll)),-1);
			vecRespDistAll = matDistAll(matSelectAll);
			
			matRespDistNoRem(intTrial,:) = vecRespDistAll;
			matHeterogeneity(intTrial,:) = vecRespDistAll;
			
			%% REMOVE MOST ACTIVE NEURONS (z-score)
			%do calculation
			vecActivity = matRespNormPerContrast(:,intTrial);
			vecActivity = findmin(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% most active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDist = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelect = tril(true(size(matDist)),-1);
			vecRespDist = matDist(matSelect);
			matRespDistHiZRem(intTrial,:) = [vecRespDist;nan(intValuesTotal-length(vecRespDist),1)];
			
			%% REMOVE LEAST ACTIVE NEURONS (z-score)
			%do calculation
			vecActivity = matRespNormPerContrast(:,intTrial);
			vecActivity = findmax(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% least active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDist = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelect = tril(true(size(matDist)),-1);
			vecRespDist = matDist(matSelect);
			matRespDistLoZRem(intTrial,:) = [vecRespDist;nan(intValuesTotal-length(vecRespDist),1)];
			
			%% REMOVE MOST ACTIVE NEURONS (df/f)
			%do calculation
			vecRawActivity = matTrialResponse(:,intTrial);
			vecActivity = matRespNormPerContrast(:,intTrial);
			[dummy,vecIndex] = findmin(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDist = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelect = tril(true(size(matDist)),-1);
			vecRespDist = matDist(matSelect);
			matRespDistHiARem(intTrial,:) = [vecRespDist;nan(intValuesTotal-length(vecRespDist),1)];
			
			%% REMOVE LEAST ACTIVE NEURONS (df/f)
			%do calculation
			vecRawActivity = matTrialResponse(:,intTrial);
			vecActivity = matRespNormPerContrast(:,intTrial);
			[dummy,vecIndex] = findmax(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% least active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDist = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelect = tril(true(size(matDist)),-1);
			vecRespDist = matDist(matSelect);
			matRespDistLoARem(intTrial,:) = [vecRespDist;nan(intValuesTotal-length(vecRespDist),1)];
		end
		
		for intRemType=1:5
			if intRemType == 1
				%no rem
				matRespDist = matRespDistNoRem;
				strRemType = 'No';
			elseif intRemType == 2
				%rem hi
				indKeep = ~isnan(matRespDistHiZRem(1,:));
				matRespDistHiZRem = matRespDistHiZRem(:,indKeep);
				matRespDist = matRespDistHiZRem;
				strRemType = 'High Z';
			elseif intRemType == 3
				%rem lo
				indKeep = ~isnan(matRespDistLoZRem(1,:));
				matRespDistLoZRem = matRespDistLoZRem(:,indKeep);
				matRespDist = matRespDistLoZRem;
				strRemType = 'Low Z';
			elseif intRemType == 4
				%rem hi
				indKeep = ~isnan(matRespDistHiARem(1,:));
				matRespDistHiARem = matRespDistHiARem(:,indKeep);
				matRespDist = matRespDistHiARem;
				strRemType = 'High A';
			elseif intRemType == 5
				%rem lo
				indKeep = ~isnan(matRespDistLoARem(1,:));
				matRespDistLoARem = matRespDistLoARem(:,indKeep);
				matRespDist = matRespDistLoARem;
				strRemType = 'Low A';
			end
			
			%{
		for intTrial=1:intTrials
			%do calculation for high response correlated neurons
			vecActivity = matRespNormPerContrast(:,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeurons]);
			matZ2 = repmat(vecActivity',[intNeurons 1]);
			matDist = abs(matZ1 - matZ2);

			figure
			imagesc(matDist)
			colormap(flipud(jet))
			colorbar
			xlabel('Neuron ID')
			ylabel('Neuron ID')
			title(sprintf('Neuronal pairwise distance in z-scored activation level\nExample trial (# %d)',intTrial))
			pause
			close
		end
			%}
			
			%plot activity dissimilarity over trials
			figure
			subplot(2,2,1)
			vecX = 1:intTrials;
			vecInvX = intTrials:-1:1;
			vecMeanH = nanmean(matRespDist,2)';
			vecErrH = nanstd(matRespDist,[],2)'/sqrt(intValuesTotal);
			vecMinTrace = vecMeanH-vecErrH;
			vecMaxTrace = vecMeanH+vecErrH;
			
			fill([vecX vecInvX],[vecMaxTrace vecMinTrace(vecInvX)],[0.7 0.7 0.7],'EdgeColor','none');
			plot(vecX,vecMeanH,'-','LineWidth',2,'Color',[0 0 0]);
			title(sprintf('%s rem;',strRemType))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%split data set
			indSelectContrasts = cellMultiSes{intPopulation}.structStim.Contrast == 0.005 | cellMultiSes{intPopulation}.structStim.Contrast == 0.02  | cellMultiSes{intPopulation}.structStim.Contrast == 0.08  | cellMultiSes{intPopulation}.structStim.Contrast == 0.32;
			indSelectResp = cellMultiSes{intPopulation}.structStim.vecTrialResponse;
			
			%hit trials only
			subplot(2,2,3)
			indSelectRespTrials = indSelectContrasts & indSelectResp;
			intRespTrials = sum(indSelectRespTrials);
			matRespDistSub = matRespDist(indSelectRespTrials,:);
			vecX = 1:intRespTrials;
			vecInvX = intRespTrials:-1:1;
			vecMeanResp = nanmean(matRespDistSub,2)';
			vecErrResp = nanstd(matRespDistSub,[],2)'/sqrt(intValuesTotal);
			vecMinTraceResp = vecMeanResp-vecErrResp;
			vecMaxTraceResp = vecMeanResp+vecErrResp;
			
			fill([vecX vecInvX],[vecMaxTraceResp vecMinTraceResp(vecInvX)],[0.7 0.7 0.7],'EdgeColor','none');
			hold on
			plot(vecX,vecMeanResp,'-','LineWidth',2,'Color',[0 0 0]);
			hold off
			dblMeanR = nanmean(matRespDistSub(:));
			dblSdR = nanstd(matRespDistSub(:));
			intValsR = length(matRespDistSub(:));
			title(sprintf('%s rem; Resp;',strRemType))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%miss trials only
			indSelectNoRespTrials = indSelectContrasts & ~indSelectResp;
			subplot(2,2,4)
			intNoRespTrials = sum(indSelectNoRespTrials);
			matNoResp = matRespDist(indSelectNoRespTrials,:);
			vecX = 1:intNoRespTrials;
			vecInvX = intNoRespTrials:-1:1;
			vecMeanNoResp = nanmean(matNoResp,2)';
			vecErrNoResp = nanstd(matNoResp,[],2)'/sqrt(intValuesTotal);
			vecMinTraceNoResp = vecMeanNoResp-vecErrNoResp;
			vecMaxTraceNoResp = vecMeanNoResp+vecErrNoResp;
			
			fill([vecX vecInvX],[vecMaxTraceNoResp vecMinTraceNoResp(vecInvX)],[0.7 0.7 0.7],'EdgeColor','none');
			hold on
			plot(vecX,vecMeanNoResp,'-','LineWidth',2,'Color',[0 0 0]);
			hold off
			dblMeanNR = nanmean(matNoResp(:));
			dblSdNR = nanstd(matNoResp(:));
			intValsNR = length(matNoResp(:));
			
			title(sprintf('%s rem; No Resp',strRemType))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%comparison of hit/miss and hit corr/uncorr
			subplot(2,2,2)
			dblErrNR = dblSdNR/sqrt(intValsNR);
			dblErrR = dblSdR/sqrt(intValsR);
			errorbar(1:2,[dblMeanNR dblMeanR],[dblErrNR dblErrR],'Linestyle','none','Marker','x');
			set(gca,'XTick',1:2,'XTickLabel',{'No Resp','Resp'})
			xlim([0.5 2.5])
			ylabel('Mean within-group z-scored activation dissimilarity')
			title(sprintf('%s rem; Block %d',strRemType,intPopulation));
			
			%put in output
			cellSaveNormActDissim{1,1,intPopulation,intRemType} = matRespDistSub(:); %within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals) with every cell = vector of dissimilarity values
			cellSaveNormActDissim{1,2,intPopulation,intRemType} = matNoResp(:);
			
			if sParams.boolSavePlots
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_%srem_raw',strSes,intPopulation,strRemType);
				export_fig([strFig '.tif']);
				export_fig([strFig '.pdf']);
			end
		end
		%no rem
		matRespDist = matRespDistNoRem;
		strRemType = 'No';
		
		%% calculate the face of god
		%loop through contrasts & stim types for trial reordering
		intStimType = 0;
		vecStimType = nan(1,length(cellSelectOri{1}));
		for intStim=unique(vecNeuronPrefStim)
			for intContrastIndex=2:length(cellSelectContrasts)
				intStimType = intStimType + 1;
				%trials
				indSelectTrials = cellSelectOri{intStim} & cellSelectContrasts{intContrastIndex};
				vecStimType(indSelectTrials) = intStimType;
			end
		end
		vecOriType = mod(vecStimType,4);
		vecOriType(vecOriType==0) = 4;
		
		%create figs
		hPopCorr = figure;
		hPopAct = figure;
		hPopITD = figure;
		
		%source matrices
		%matTrialResponse(intNeuron,intTrial);
		%matRespNormPerContrastintNeuron,intTrial);
		
		
		%for whole pop or only preferred pop
		for intPopType=1:3
			if intPopType == 1
				strPopType = 'AllNeurons';
			elseif intPopType == 2
				strPopType = 'PrefNeurons';
			elseif intPopType == 3
				strPopType = 'NonPrefNeurons';
			end
			if intPopType == 1 %all neurons
				%calculate trial-by-trial correlations of population activity
				matCorr=corr(matRespNormPerContrast(:,:));
				[vecHet,vecAct] = calcMatRespHeteroGen(matTrialResponse);
			else %only pref-pop
				for intPrefType=1:4
					indSelect = vecNeuronPrefStim' == intPrefType;
					if intPopType == 3
						indSelect = ~indSelect;
					end
					%strPopType
					%sum(indSelect)
					
					%calculate trial-by-trial correlations of population activity
					if sum(indSelect) > 0
						cellCorr{intPrefType}=corr(matRespNormPerContrast(indSelect,:));
						[cellHet{intPrefType},cellAct{intPrefType}] = calcMatRespHeteroGen(matRespNormPerContrast(indSelect,:));
					else
						cellCorr{intPrefType} = nan(intTrials,intTrials);
						cellHet{intPrefType} = nan(intTrials,1);
						cellAct{intPrefType} = nan(intTrials,1);
					end
				end
				matCorr = nan(size(cellCorr{1}));
				vecHet = nan(size(cellHet{1}));
				vecAct = nan(size(cellAct{1}));
				for intOriType=vecOriType
					if isnan(intOriType),continue;end
					%get trials
					indSelectTrials = vecOriType == intOriType;
					%assign
					matCorr(indSelectTrials,indSelectTrials) = cellCorr{intPrefType}(indSelectTrials,indSelectTrials);
					vecHet(indSelectTrials) = cellHet{intPrefType}(indSelectTrials);
					vecAct(indSelectTrials) = cellAct{intPrefType}(indSelectTrials);
				end
			end
			
			
			%compute absolute trial distance
			matTrialDist = abs(repmat(1:intTrials,[intTrials 1])-repmat((1:intTrials)',[1 intTrials]));
			%matSelect = tril(true(size(matTrialDist)),-1);
			indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
			indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
			indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			
			%separate across-trial correlation correlations for different trial types (1-4)
			%also split for hit/miss trials
			intHitTrials = sum(indHit);
			intMissTrials = sum(indMiss);
			intFastTrials = sum(indFast);
			intSlowTrials = sum(indSlow);
			vecDistHit = [];
			vecDistMiss = [];
			vecDistFast = [];
			vecDistSlow = [];
			vecCorrHit = [];
			vecCorrMiss = [];
			vecCorrFast = [];
			vecCorrSlow = [];
			
			vecHetHit = [];
			vecHetMiss = [];
			vecHetFast = [];
			vecHetSlow = [];
						
			vecActHit = [];
			vecActMiss = [];
			vecActFast = [];
			vecActSlow = [];
			for intOriType=vecOriType
				%% hit; select trials
				if isnan(intOriType),continue;end
				indSelectHitTrials = vecOriType == intOriType & indHit;
				if sum(indSelectHitTrials) > 2
					%get trials
					matDistSubHit = matTrialDist(indSelectHitTrials,indSelectHitTrials);
					matCorrSubHit = matCorr(indSelectHitTrials,indSelectHitTrials);
					
					%subselection matrix hit
					matSubSelectHit = tril(true(size(matDistSubHit)),-1);
					
					%subselect
					vecDistHit = [vecDistHit matDistSubHit(matSubSelectHit)'];
					vecCorrHit = [vecCorrHit matCorrSubHit(matSubSelectHit)'];
					vecHetHit = [vecHetHit vecHet(indSelectHitTrials)'];
					vecActHit = [vecActHit vecAct(indSelectHitTrials)'];
				end
				
				%% miss; select trials
				indSelectMissTrials = vecOriType == intOriType & indMiss;
				if sum(indSelectMissTrials) > 2
					%get trials
					matDistSubMiss = matTrialDist(indSelectMissTrials,indSelectMissTrials);
					matCorrSubMiss = matCorr(indSelectMissTrials,indSelectMissTrials);
					
					%subselection matrix hit
					matSubSelectMiss = tril(true(size(matDistSubMiss)),-1);
					
					%subselect
					vecDistMiss = [vecDistMiss matDistSubMiss(matSubSelectMiss)'];
					vecCorrMiss = [vecCorrMiss matCorrSubMiss(matSubSelectMiss)'];
					vecHetMiss = [vecHetMiss vecHet(indSelectMissTrials)'];
					vecActMiss = [vecActMiss vecAct(indSelectMissTrials)'];
				end
				
				%% fast; select trials
				indSelectFastTrials = vecOriType == intOriType & indFast;
				if sum(indSelectFastTrials) > 2
					%get trials
					matDistSubFast = matTrialDist(indSelectFastTrials,indSelectFastTrials);
					matCorrSubFast = matCorr(indSelectFastTrials,indSelectFastTrials);
					
					%subselection matrix hit
					matSubSelectFast = tril(true(size(matDistSubFast)),-1);
					
					%subselect
					vecDistFast = [vecDistFast matDistSubFast(matSubSelectFast)'];
					vecCorrFast = [vecCorrFast matCorrSubFast(matSubSelectFast)'];
					vecHetFast = [vecHetFast vecHet(indSelectFastTrials)'];
					vecActFast = [vecActFast vecAct(indSelectFastTrials)'];
				end
				
				%% slow; select trials
				indSelectSlowTrials = vecOriType == intOriType & indSlow;
				if sum(indSelectSlowTrials) > 2
					%get trials
					matDistSubSlow = matTrialDist(indSelectSlowTrials,indSelectSlowTrials);
					matCorrSubSlow = matCorr(indSelectSlowTrials,indSelectSlowTrials);
					
					%subselection matrix hit
					matSubSelectSlow = tril(true(size(matDistSubSlow)),-1);
					
					%subselect
					vecDistSlow = [vecDistSlow matDistSubSlow(matSubSelectSlow)'];
					vecCorrSlow = [vecCorrSlow matCorrSubSlow(matSubSelectSlow)'];
					vecHetSlow = [vecHetSlow vecHet(indSelectSlowTrials)'];
					vecActSlow = [vecActSlow vecAct(indSelectSlowTrials)'];
				end
			end
			
			%plot
			figure(hPopCorr);
			subplot(3,2,1+(intPopType-1)*2)
			matPlotCorr = matCorr(~isnan(vecOriType),~isnan(vecOriType));
			matPlotCorr(diag(true(1,length(matPlotCorr)))) = nan;
			vecLim = [-nanmax(abs(matPlotCorr(:))) nanmax(abs(matPlotCorr(:)))];
			imagesc(matPlotCorr,vecLim);
			matC = redbluepurple(128);
			nancolorbar(matPlotCorr,vecLim,matC);
			title(sprintf('Inter-trial correlation of %s population activity (block %d)',strPopType,intPopulation))
			xlabel('Trial')
			ylabel('Trial')
			
			subplot(3,2,2+(intPopType-1)*2)
			dblMeanCorrHit = mean(vecCorrHit);
			dblMeanCorrMiss = mean(vecCorrMiss);
			dblMeanCorrFast = mean(vecCorrFast);
			dblMeanCorrSlow = mean(vecCorrSlow);
			dblErrCorrHit = std(vecCorrHit)/sqrt(length(vecCorrHit));
			dblErrCorrMiss = mean(vecCorrMiss)/sqrt(length(vecCorrMiss));
			dblErrCorrFast = mean(vecCorrFast)/sqrt(length(vecCorrFast));
			dblErrCorrSlow = mean(vecCorrSlow)/sqrt(length(vecCorrSlow));
			
			
			[h,pHM] = ttest2(vecCorrHit,vecCorrMiss);
			[h,pMF] = ttest2(vecCorrMiss,vecCorrFast);
			[h,pMS] = ttest2(vecCorrMiss,vecCorrSlow);
			[h,pFS] = ttest2(vecCorrFast,vecCorrSlow);
			
			errorbar(1:4,[dblMeanCorrMiss dblMeanCorrHit dblMeanCorrFast dblMeanCorrSlow],[dblErrCorrMiss dblErrCorrHit dblErrCorrFast dblErrCorrSlow],'Linestyle','none','Marker','x','Color','k');
			set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
			title(sprintf('Mean inter-trial correlation of %s population activity\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
			ylabel('Pearson correlation of population activity')
			
			%plot heterogeneity + activity
			figure(hPopAct);

			subplot(3,2,1+(intPopType-1)*2)
			
			[h,pHM] = ttest2(vecHetHit,vecHetMiss);
			[h,pMF] = ttest2(vecHetMiss,vecHetFast);
			[h,pMS] = ttest2(vecHetMiss,vecHetSlow);
			[h,pFS] = ttest2(vecHetFast,vecHetSlow);
			
			errorbar(1:4,[mean(vecHetMiss) mean(vecHetHit) mean(vecHetFast) mean(vecHetSlow)],...
				[std(vecHetMiss)/sqrt(length(vecHetMiss)) std(vecHetHit)/sqrt(length(vecHetHit)) std(vecHetFast)/sqrt(length(vecHetFast)) std(vecHetSlow)/sqrt(length(vecHetSlow))],...
				'Linestyle','none','Marker','x','Color','k');
			set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
			title(sprintf('Mean heterogeneity of %s population\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
			ylabel('Heterogeneity')
			
			subplot(3,2,2+(intPopType-1)*2)
			
			[h,pHM] = ttest2(vecActHit,vecActMiss);
			[h,pMF] = ttest2(vecActMiss,vecActFast);
			[h,pMS] = ttest2(vecActMiss,vecActSlow);
			[h,pFS] = ttest2(vecActFast,vecActSlow);
			
			errorbar(1:4,[mean(vecActMiss) mean(vecActHit) mean(vecActFast) mean(vecActSlow)],...
				[std(vecActMiss)/sqrt(length(vecActMiss)) std(vecActHit)/sqrt(length(vecActHit)) std(vecActFast)/sqrt(length(vecActFast)) std(vecActSlow)/sqrt(length(vecActSlow))],...
				'Linestyle','none','Marker','x','Color','k');
			set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
			title(sprintf('Mean dF/F0 activity of %s population\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
			ylabel('Mean dF/F0 activity')
			
			
			%% plot inter-trial correlation vs distance; + relation dF/F0 vs heterogeneity
			figure(hPopITD)
			
			%plot ITD
			subplot(3,2,1+(intPopType-1)*2)
			
			intStep = 25;
			intMax = ceil(intTrials/intStep)*intStep;
			vecBinX = 0:intStep:intMax;
			vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
			[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHit,vecBinX);
			errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'r')
			sStats=regstats(vecCorrHit,vecDistHit,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStats.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
			ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
			xlabel('Inter-trial distance')
			
			
			%plot relation dF/F0 vs heterogeneity
			subplot(3,2,2+(intPopType-1)*2)
			scatter(vecActMiss,vecHetMiss,'xr')
			hold on
			scatter(vecActFast,vecHetFast,'xg')
			scatter(vecActSlow,vecHetSlow,'xy')
			hold off
			legend(gca,'Miss','Hit-Fast','Hit-Slow')
			xlim([-0.05 0.15])
			xlabel('Mean population dF/F0')
			ylabel('Heterogeneity')
			sStatsM=regstats(vecActMiss,vecHetMiss,'linear',{'beta','rsquare','tstat'});
			sStatsF=regstats(vecActFast,vecHetFast,'linear',{'beta','rsquare','tstat'});
			sStatsS=regstats(vecActSlow,vecHetSlow,'linear',{'beta','rsquare','tstat'});
			title(sprintf('Lin reg: Miss; slope=%.3f; Fast; slope=%.3f; Slow; slope=%.3f',...
				sStatsM.beta(2),sStatsF.beta(2),sStatsS.beta(2)))
			
			%%%% OUTPUT
			cellSaveCorrITD{1,1,intPopType,intPopulation} = vecDistMiss;
			cellSaveCorrITD{1,2,intPopType,intPopulation} = vecCorrMiss;
			cellSaveCorrITD{1,3,intPopType,intPopulation} = vecHetMiss;
			cellSaveCorrITD{1,4,intPopType,intPopulation} = vecActMiss;
			
			cellSaveCorrITD{2,1,intPopType,intPopulation} = vecDistHit;
			cellSaveCorrITD{2,2,intPopType,intPopulation} = vecCorrHit;
			cellSaveCorrITD{2,3,intPopType,intPopulation} = vecHetHit;
			cellSaveCorrITD{2,4,intPopType,intPopulation} = vecActHit;
			
			cellSaveCorrITD{3,1,intPopType,intPopulation} = vecDistSlow;
			cellSaveCorrITD{3,2,intPopType,intPopulation} = vecCorrSlow;
			cellSaveCorrITD{3,3,intPopType,intPopulation} = vecHetSlow;
			cellSaveCorrITD{3,4,intPopType,intPopulation} = vecActSlow;
			
			cellSaveCorrITD{4,1,intPopType,intPopulation} = vecDistFast;
			cellSaveCorrITD{4,2,intPopType,intPopulation} = vecCorrFast;
			cellSaveCorrITD{4,3,intPopType,intPopulation} = vecHetFast;
			cellSaveCorrITD{4,4,intPopType,intPopulation} = vecActFast;
		end
		%plot
		if sParams.boolSavePlots
			figure(hPopCorr);
			drawnow;
			jFig = get(handle(hPopCorr), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopCorr);
			drawnow;
			strFig = sprintf('%sagg_intertrialcorr_%s_pop%d__raw',strSes,strPopType,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
			
			figure(hPopAct);
			drawnow;
			jFig = get(handle(hPopAct), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopAct);
			drawnow;
			strFig = sprintf('%sagg_HetAct_%s_pop%d__raw',strSes,strPopType,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
			
			figure(hPopITD);
			drawnow;
			jFig = get(handle(hPopITD), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopITD);
			drawnow;
			strFig = sprintf('%sagg_intertrialcorrDistDep_%s_pop%d__raw',strSes,strPopType,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% plot dependency on reaction times
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs(indSelectRespTrials) - cellMultiSes{intPopulation}.structStim.SecsOn(indSelectRespTrials); %get RTs for hit trials and selected contrasts
		[vecRTsSorted,vecRTSortIndex] = sort(vecRTs,'ascend');
		matHitRTSorted = matRespDistSub(vecRTSortIndex,:);
		vecHitRT = nanmean(matHitRTSorted,2);
		
		%Z activity
		matActZ = matRespNormPerContrast(:,indSelectRespTrials)';

		matZRTSorted = matActZ(vecRTSortIndex,:);
		vecZRT = mean(matZRTSorted,2);
		
		%dF/F
		matActdFoF = matTrialResponse(:,indSelectRespTrials)';

		matARTSorted = matActdFoF(vecRTSortIndex,:);
		vecART = mean(matARTSorted,2);
		
		%variance
		vecVRT = var(matARTSorted,[],2);
		
		%save
		cellSaveRTDependency{1,intPopulation} = vecRTsSorted;
		cellSaveRTDependency{2,intPopulation} = vecHitRT;
		cellSaveRTDependency{3,intPopulation} = vecZRT;
		cellSaveRTDependency{4,intPopulation} = vecART;
		cellSaveRTDependency{5,intPopulation} = vecVRT;
		
		%plot
		figure
		subplot(2,2,1);
		scatter(vecRTsSorted,vecHitRT,'kx')
		
		%perform regressions
		sStatsC=regstats(vecHitRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		
		title(sprintf('Pop heterogeneity; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean population activity dissimilarity')
		
		%z-scored activity
		subplot(2,2,2);
		scatter(vecRTsSorted,vecZRT,'kx')
		
		%perform regressions
		sStatsC=regstats(vecZRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		
		title(sprintf('Z-scored act; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean z-scored population activity')
		
		%dF/F activity
		subplot(2,2,3);
		scatter(vecRTsSorted,vecART,'kx')

		%perform regressions
		sStatsC=regstats(vecART,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		
		title(sprintf('dF/F0; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean dF/F population activity')
		
		%dF/F activity
		subplot(2,2,4);
		scatter(vecRTsSorted,vecVRT,'kx')

		%perform regressions
		sStatsC=regstats(vecVRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		
		title(sprintf('Variance; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Variance of population activity')
		
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_behaviorallycorrelated_act_z_dissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% heterogeneity for different contrasts, split for hit/miss (like dF/F)
		%figure
		hHetCon = figure;
		set(hHetCon,'Color',[1 1 1]);
		figure(hHetCon);
		
		%pre-allocate
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrasts(1) = 0.2;
		intContrasts = length(vecContrasts);
		vecP = nan(1,intContrasts);
		vecMetaHitY = nan(1,intContrasts);
		vecMetaHitE = nan(1,intContrasts);
		vecMetaMissY = nan(1,intContrasts);
		vecMetaMissE = nan(1,intContrasts);
		for intContrastIndex=1:intContrasts
			%get contrast
			dblContrast = vecContrasts(intContrastIndex);

			%get resp
			vecHeteroHit = matHeterogeneity(cellSelectContrasts{intContrastIndex} & indSelectResp,:);
			vecHeteroHit = vecHeteroHit(:);
			vecHeteroMiss = matHeterogeneity(cellSelectContrasts{intContrastIndex} & ~indSelectResp,:);
			vecHeteroMiss = vecHeteroMiss(:);
			
			%perform analyses per contrast
			[h,p,ci] = ttest2(vecHeteroHit,vecHeteroMiss);
			vecY = [mean(vecHeteroHit) mean(vecHeteroMiss)];
			
			%put in output
			cellSaveMatContHetero{intPopulation}(intContrastIndex,1) = vecY(1);
			cellSaveMatContHetero{intPopulation}(intContrastIndex,2) = vecY(2);
			
			%put in meta vector
			vecP(intContrastIndex) = p;
			vecMetaHitY(intContrastIndex) = vecY(1);
			vecMetaHitE(intContrastIndex) = std(vecHeteroHit)/sqrt(length(vecHeteroHit));
			vecMetaMissY(intContrastIndex) = vecY(2);
			vecMetaMissE(intContrastIndex) = std(vecHeteroMiss)/sqrt(length(vecHeteroMiss));
		end
		
		%pre-compute variables
		vecWindow = [1 length(vecContrasts)];
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		
		%get data
		for intResp=[0 1]
			if intResp == 1
				vecMeanTrace = vecMetaHitY(vecWindowSelect);
				vecSE = vecMetaHitE(vecWindowSelect);
				vecColorFill = [0.7 1 0.7];
				vecColorLine = [0 1 0];
			else
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
				vecMeanTrace = vecMetaMissY(vecWindowSelect);
				vecSE = vecMetaMissE(vecWindowSelect);
			end
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
		set(gca,'XScale','log','YScale','linear')
		set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
		title(sprintf('Pop het [%d]; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',intPopulation,vecP))
		grid on
		xlabel('Contrast')
		ylabel('Population response heterogeneity')
		xlim(vecContrasts(vecWindow))
		%ylim([-0.01 0.06])
		legend({'SEM','Miss','SEM','Hit'},'Location','Best')
		
		drawnow;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_heterogeneity_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		
		%% dF/F over time + heterogeneity over time; also split for fast/slow trials
		return
		sParamsHet.intWindowLength = 1;
		[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		
		%get stim data
		cellFieldsC = {'Contrast'};
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation},cellFieldsC);
		cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		vecC = sTypesC.matTypes;
		
		%get resp data
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		
		%general
		vecWindowSecs = [-3 5];
		vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
		vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
		vecWindowInv = length(vecLineX):-1:1;
		vecX = [vecLineX vecLineX(vecWindowInv)];
		
		%plot
		cellSaveActTime = {};
		cellSaveHetTime = {};
		hHetTime = figure;
		hActTime = figure;
		for intC=1:length(vecC)
			indTrials = cellSelectC{intC};
			for intAct=[0 1]
				if intAct == 0,figure(hHetTime);else figure(hActTime);end
				subplot(2,3,intC);
				for intType=1:3
					if intType == 1 %miss
						vecTrials = find(indTrials&indMiss);
						strType = 'Miss';
						vecColorLine = [1 0 0];
						vecColorFill = [1 0.7 0.7];
					elseif intType == 2 %slow
						vecTrials = find(indTrials&indSlow);
						strType = 'Slow';
						vecColorLine = [1 1 0];
						vecColorFill = [1 1 0.7];
					else %fast
						vecTrials = find(indTrials&indFast);
						strType = 'Fast';
						vecColorLine = [0 1 0];
						vecColorFill = [0.7 1 0.7];
					end
					matRawData = zeros(length(vecTrials),length(vecWindow(1):vecWindow(end)));
					
					%get data
					vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn(vecTrials);
					for intTrial=1:length(vecTrials)
						intStart = vecStimOn(intTrial)+vecWindow(1);
						intStop = vecStimOn(intTrial)+vecWindow(end);
						if intAct == 0,
							matRawData(intTrial,:) = vecHeterogeneity(intStart:intStop);
						else
							matRawData(intTrial,:) = vecActivity_dFoF(intStart:intStop);
						end
						
					end
					cellType{intType} = matRawData;
					
					vecMeanTrace = mean(matRawData,1);
					vecSE = std(matRawData,[],1)/sqrt(length(vecTrials));
					vecMinTrace = vecMeanTrace-vecSE;
					vecMaxTrace = vecMeanTrace+vecSE;
					vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
					
					%plot
					hold on
					errorfill(vecLineX,vecMeanTrace,vecSE,vecColorLine,vecColorFill);
					%fill(vecX,vecY,vecColorFill,'EdgeColor','none');
					%plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
					hold off
				end
				
				%save data&labels
				if intAct == 0,
					cellSaveHetTime{intC,intPopulation} = cellType;
					title(sprintf('Pop Het [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
					ylabel('Population response heterogeneity')
				else
					cellSaveActTime{intC,intPopulation} = cellType;
					title(sprintf('Pop Act [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
					ylabel('Population dF/F')
				end
				grid on
				xlabel('Time after stim onset (s)')
				xlim(vecWindowSecs)
				%ylim([-0.01 0.06])
				drawnow;
			end
		end
		
		if sParams.boolSavePlots
			figure(hHetTime)
			drawnow;
			strFig = sprintf('%s_heterogeneity_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
			
			figure(hActTime)
			drawnow;
			strFig = sprintf('%s_activity_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% decode if trial will be fast, slow or miss based on dF/F and heterogeneity
		%{
	per animal, sum miss, slow and fast trials; trake min(m,s,f)-1 trials for
	likelihood construction and take one other random trial per condition to
	decode [make procedure work over all time points]
	=> then bootstrap procedure
		%}
		
		%get data
		if ~exist('vecHeterogeneity','var') || isempty(vecHeterogeneity)
			sParamsHet.intWindowLength = 1;%round(cellMultiSes{1}.samplingFreq);
			[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		end
		
		%get resp data
		indMiss = isnan(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		
		%general
		vecWindowSecs = [-3 5];
		vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
		vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
		vecSelectSecs = [-3 0];
		vecSelect = round(vecSelectSecs*cellMultiSes{intPopulation}.samplingFreq);
		intStart = find(vecLineX>=vecSelectSecs(1),1);
		intStop = find(vecLineX>=vecSelectSecs(2),1);
		intFramesPerStim = length(vecWindow(1):vecWindow(end));
		vecSelectFrames = intStart:intStop;
		
		%save vars
		hDecRespTypeTime = figure;
		
		%get nr of trials
		intMiss = sum(indMiss);
		vecMiss = find(indMiss);
		intFast = sum(indFast);
		vecFast = find(indFast);
		intSlow = sum(indSlow);
		vecSlow = find(indSlow);
		intNrTrials = min([intMiss intFast intSlow]);
		
		%transform data to easy-access format
		%pre-allocate
		matRawDataAct = zeros(length(indMiss),intFramesPerStim);
		matRawDataHet = zeros(length(indMiss),intFramesPerStim);
		
		%get data
		vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn;
		for intTrial=1:length(indMiss)
			intStart = vecStimOn(intTrial)+vecWindow(1);
			intStop = vecStimOn(intTrial)+vecWindow(end);
			
			matRawDataHet(intTrial,:) = vecHeterogeneity(intStart:intStop);
			matRawDataAct(intTrial,:) = vecActivity_dFoF(intStart:intStop);
		end
		
		
		%bootstrap decoding
		intIters = 1000;
		
		%pre-allocate output
		matPostAct = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
		matPostHet = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
		
		%loop
		for intIter=1:intIters
			%get random trials
			vecMissRand = randperm(intMiss,intNrTrials);
			vecMissLikelihood = vecMiss(vecMissRand(1:(end-1)));
			intDecodeTrialMiss = vecMiss(vecMissRand(end));
			
			vecFastRand = randperm(intFast,intNrTrials);
			vecFastLikelihood = vecFast(vecFastRand(1:(end-1)));
			intDecodeTrialFast = vecFast(vecFastRand(end));
			
			vecSlowRand = randperm(intSlow,intNrTrials);
			vecSlowLikelihood = vecSlow(vecSlowRand(1:(end-1)));
			intDecodeTrialSlow = vecSlow(vecSlowRand(end));
			
			%get likelihoods df/f
			dblLikeActMissMean = mean(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
			dblLikeActMissSD = std(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
			
			dblLikeActFastMean = mean(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
			dblLikeActFastSD = std(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
			
			dblLikeActSlowMean = mean(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
			dblLikeActSlowSD = std(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
			
			%get likelihoods heterogeneity
			dblLikeHetMissMean = mean(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
			dblLikeHetMissSD = std(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
			
			dblLikeHetFastMean = mean(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
			dblLikeHetFastSD = std(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
			
			dblLikeHetSlowMean = mean(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
			dblLikeHetSlowSD = std(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
			
			
			%get act probabilities for miss trial
			dblActMissTrial = mean(matRawDataAct(intDecodeTrialMiss,vecSelectFrames));
			matPostAct(1,1,intIter) = normpdf(dblActMissTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(1,2,intIter) = normpdf(dblActMissTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(1,3,intIter) = normpdf(dblActMissTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			%get act probabilities for slow trial
			dblActSlowTrial = mean(matRawDataAct(intDecodeTrialSlow,vecSelectFrames));
			matPostAct(2,1,intIter) = normpdf(dblActSlowTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(2,2,intIter) = normpdf(dblActSlowTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(2,3,intIter) = normpdf(dblActSlowTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			%get act probabilities for fast trial
			dblActFastTrial = mean(matRawDataAct(intDecodeTrialFast,vecSelectFrames));
			matPostAct(3,1,intIter) = normpdf(dblActFastTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(3,2,intIter) = normpdf(dblActFastTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(3,3,intIter) = normpdf(dblActFastTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			
			%get het probabilities for miss trial
			vecHetMissTrial = mean(matRawDataHet(intDecodeTrialMiss,vecSelectFrames));
			matPostHet(1,1,intIter) = normpdf(vecHetMissTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(1,2,intIter) = normpdf(vecHetMissTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(1,3,intIter) = normpdf(vecHetMissTrial,dblLikeHetFastMean,dblLikeHetFastSD);
			
			%get het probabilities for slow trial
			vecHetSlowTrial = mean(matRawDataHet(intDecodeTrialSlow,vecSelectFrames));
			matPostHet(2,1,intIter) = normpdf(vecHetSlowTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(2,2,intIter) = normpdf(vecHetSlowTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(2,3,intIter) = normpdf(vecHetSlowTrial,dblLikeHetFastMean,dblLikeHetFastSD);
			
			%get het probabilities for fast trial
			vecHetFastTrial = mean(matRawDataHet(intDecodeTrialFast,vecSelectFrames));
			matPostHet(3,1,intIter) = normpdf(vecHetFastTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(3,2,intIter) = normpdf(vecHetFastTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(3,3,intIter) = normpdf(vecHetFastTrial,dblLikeHetFastMean,dblLikeHetFastSD);
		end
		
		%normalize probabilities
		matPostHet(:,1,:) = matPostHet(:,1,:)/mean(reshape(matPostHet(:,1,:),[1 numel(matPostHet(:,1,:))]));
		matPostHet(:,2,:) = matPostHet(:,2,:)/mean(reshape(matPostHet(:,2,:),[1 numel(matPostHet(:,2,:))]));
		matPostHet(:,3,:) = matPostHet(:,3,:)/mean(reshape(matPostHet(:,3,:),[1 numel(matPostHet(:,3,:))]));
		matPostAct(:,1,:) = matPostAct(:,1,:)/mean(reshape(matPostAct(:,1,:),[1 numel(matPostAct(:,1,:))]));
		matPostAct(:,2,:) = matPostAct(:,2,:)/mean(reshape(matPostAct(:,2,:),[1 numel(matPostAct(:,2,:))]));
		matPostAct(:,3,:) = matPostAct(:,3,:)/mean(reshape(matPostAct(:,3,:),[1 numel(matPostAct(:,3,:))]));
		
		%get fraction correct
		%act
		[dummy,matDecodedActMiss]=max(matPostAct(1,:,:),[],2);
		vecMissActCorrect = sum(squeeze(matDecodedActMiss == 1))/intIters;
		[dummy,matDecodedActSlow]=max(matPostAct(2,:,:),[],2);
		vecSlowActCorrect = sum(squeeze(matDecodedActSlow == 2))/intIters;
		[dummy,matDecodedActFast]=max(matPostAct(3,:,:),[],2);
		vecFastActCorrect = sum(squeeze(matDecodedActFast == 3))/intIters;
		
		%het
		[dummy,matDecodedHetMiss]=max(matPostHet(1,:,:),[],2);
		vecMissHetCorrect = sum(squeeze(matDecodedHetMiss == 1))/intIters;
		[dummy,matDecodedHetSlow]=max(matPostHet(2,:,:),[],2);
		vecSlowHetCorrect = sum(squeeze(matDecodedHetSlow == 2))/intIters;
		[dummy,matDecodedHetFast]=max(matPostHet(3,:,:),[],2);
		vecFastHetCorrect = sum(squeeze(matDecodedHetFast == 3))/intIters;
		
		
		%% plot decoding performance
		%{
	make angle for m/s/f separated by 2/3pi. Calc resultant vector angle of
	triplet decoding, graph aspolar plot
		%}
		cellSaveRespDecode{intPopulation,1} = zeros(2,2,3);%vectors; [act/het] x [theta/rho] x [miss/slow/fast]
		cellSaveRespDecode{intPopulation,2} = zeros(2,3);%distance; [act/het] x [miss/slow/fast]
		for intRespType=1:3
			if intRespType == 1
				strTitle = 'Miss';
				matDecodedAct = matDecodedActMiss;
				matDecodedHet = matDecodedHetMiss;
				thetaCorr = (4/3)*pi;
				cellColor = {'r','k','k'};
			elseif intRespType == 2
				strTitle = 'Slow';
				matDecodedAct = matDecodedActSlow;
				matDecodedHet = matDecodedHetSlow;
				thetaCorr = (0/3)*pi;
				cellColor = {'k','y','k'};
			elseif intRespType == 3
				strTitle = 'Fast';
				matDecodedAct = matDecodedActFast;
				matDecodedHet = matDecodedHetFast;
				thetaCorr = (2/3)*pi;
				cellColor = {'k','k','g'};
			end
			subplot(2,2,intRespType)
			%calc corr vector
			[xCorr,yCorr] = pol2cart(thetaCorr,1);
			
			%calc resultant vector act
			[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedAct == 1))/intIters);%miss
			[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedAct == 2))/intIters);%slow
			[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedAct == 3))/intIters);%fast
			x=x1+x2+x3;
			y=y1+y2+y3;
			[thetaAct,rhoAct]=cart2pol(x,y);
			[thetaDiffAct,rhoDiffAct] = cart2pol(xCorr-x,yCorr-y);
			
			%resultant vector het
			[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedHet == 1))/intIters);%miss
			[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedHet == 2))/intIters);%slow
			[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedHet == 3))/intIters);%fast
			x=x1+x2+x3;
			y=y1+y2+y3;
			[thetaHet,rhoHet]=cart2pol(x,y);
			[thetaDiffHet,rhoDiffHet] = cart2pol(xCorr-x,yCorr-y);
			
			%save data
			cellSaveRespDecode{intPopulation,1}(1,1,intRespType) = thetaAct;%vectors;[act/het] x [theta/rho] x [miss/slow/fast]
			cellSaveRespDecode{intPopulation,1}(1,2,intRespType) = rhoAct;
			cellSaveRespDecode{intPopulation,1}(2,1,intRespType) = thetaHet;
			cellSaveRespDecode{intPopulation,1}(2,2,intRespType) = rhoHet;
			cellSaveRespDecode{intPopulation,2}(1,intRespType) = rhoDiffAct;%distance; [act/het] x [miss/slow/fast]
			cellSaveRespDecode{intPopulation,2}(2,intRespType) = rhoDiffHet;%distance; [act/het] x [miss/slow/fast]
			
			%make polar background
			[x,y] = pol2cart((4/3)*pi,0.9);%miss
			text(x,y,'Miss','Color',cellColor{1});
			hold on
			[x,y] = pol2cart((0/3)*pi,0.7);%slow
			text(x,y,'Slow','Color',cellColor{2});
			[x,y] = pol2cart((2/3)*pi,0.9);%fast
			text(x,y,'Fast','Color',cellColor{3});
			[x,y] = pol2cart((1/3)*pi,1);%slow/fast
			plot([0 x],[0 y],'k');
			[x,y] = pol2cart((3/3)*pi,1);%miss/fast
			plot([0 x],[0 y],'k');
			[x,y] = pol2cart((5/3)*pi,1);%slow/miss
			plot([0 x],[0 y],'k');
			ang=0:0.001:2*pi;
			x=cos(ang);
			y=sin(ang);
			plot(x,y, 'k');
			
			%plot data
			[x,y]=pol2cart(thetaAct,rhoAct);
			scatter(x,y,'bx');
			[x,y]=pol2cart(thetaHet,rhoHet);
			scatter(x,y,'kx');
			hold off;
			xlim([-1 1]);
			ylim([-1 1]);
			title(sprintf('%s trials',strTitle))
		end
		
		%plot summary
		subplot(2,2,4);
		matPerformance = imnorm(1-abs(cellSaveRespDecode{intPopulation,2}));
		vecActPerfomance = matPerformance(1,:);
		vecHetPerfomance = matPerformance(2,:);
		dblOffset=0.1;
		scatter((1-dblOffset)*ones(size(vecActPerfomance)),vecActPerfomance,'bx')
		hold on
		dblMeanAct=mean(vecActPerfomance);
		dblSDAct = std(vecActPerfomance);
		errorbar((1+dblOffset),dblMeanAct,dblSDAct/sqrt(length(vecActPerfomance)),'xb')
		dblMeanHet=mean(vecHetPerfomance);
		dblSDHet = std(vecHetPerfomance);
		scatter((2-dblOffset)*ones(size(vecHetPerfomance)),vecHetPerfomance,'kx')
		errorbar((2+dblOffset),dblMeanHet,dblSDHet/sqrt(length(vecHetPerfomance)),'xk')
		hold off
		title('Normalized resp type decoding')
		ylabel('Norm decod perf')
		set(gca,'XTick',[1 2],'XTickLabel',{'dF/F0', 'Heterogen'})
		
		%save
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_resptypedecoding_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% plot heterogeneity vs correlation ideal/actual pop response (pref stim)
		%matTrialResponse(intNeuron,intTrial)
		[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matTrialResponse);
		vecOriTrial = cellMultiSes{intPopulation}.structStim.Orientation;
		vecOriLookup = unique(vecOriTrial);
		vecNeuronPrefOri = vecOriLookup(vecNeuronPrefStim);

		%run through trials
		vecPrefRespCorr = nan(intTrials,1);
		for intTrial=1:intTrials
			vecIdealResp = double(vecNeuronPrefOri == vecOriTrial(intTrial))';
			vecActualResp = matTrialResponse(:,intTrial);
			vecPrefRespCorr(intTrial) = corr(vecIdealResp,vecActualResp);
		end
		
		%plot
		hHetPopCorr = figure;
		scatter(vecHeterogeneity,vecPrefRespCorr);
		
		sStats=regstats(vecPrefRespCorr,vecHeterogeneity,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStats.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		
		title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		
		xlabel('Heterogeneity')
		ylabel('Correlation ideal/actual population response')
		
		%save output
		cellSaveHetPoprespCorr{1,intPopulation} = vecHeterogeneity;
		cellSaveHetPoprespCorr{2,intPopulation} = vecPrefRespCorr;
		cellSaveHetPoprespCorr{3,intPopulation} = sStats.beta;
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_heterogen_vs_ActualIdeal_popresp%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% pairwise neuronal correlations per stimulus repetition
		
		%source matrices
		%matTrialResponse(intNeuron,intTrial);
		%matRespNormPerContrastintNeuron,intTrial);
		%matStimResponse(intStimType,intRepetition,intNeuron)
		
		sTypesAll = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation','Contrast'});
		cellSelectAll = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesAll);
		matStimResponse = doMatRespTransform(matTrialResponse,cellSelectAll);
		
		%transform response & take mean over probe trials
		matStimResponse(1,:,:) = mean(matStimResponse(1:4,:,:),1);
		matStimResponse(2:4,:,:) = [];
		
		%calculate correlations
		intNeurons = size(matStimResponse,3);
		intPairs = ((intNeurons*intNeurons)-intNeurons)/2;
		intRepetitions = size(matStimResponse,2);
		matCorrelations = nan(intPairs,intRepetitions);
		matLocationsX = nan(intPairs,intRepetitions,2);
		matLocationsY = nan(intPairs,intRepetitions,2);
		for intRepetition=1:intRepetitions
			intPair = 0;
			for intNeuron1=1:(intNeurons-1)
				for intNeuron2=(intNeuron1+1):intNeurons
					intPair = intPair + 1;
					matCorrelations(intPair,intRepetition) = corr(matStimResponse(:,intRepetition,intNeuron1),matStimResponse(:,intRepetition,intNeuron2));
					matLocationsX(intPair,intRepetition,:) = [cellMultiSes{intPopulation}.neuron(intNeuron1).x cellMultiSes{intPopulation}.neuron(intNeuron2).x];
					matLocationsY(intPair,intRepetition,:) = [cellMultiSes{intPopulation}.neuron(intNeuron1).y cellMultiSes{intPopulation}.neuron(intNeuron2).y];
				end
			end
		end
		
		%plot
		hMatCorrFig=figure;
		subplot(2,2,1)
		imagesc(matTrialResponse);
		colormap(hot(256));colorbar;drawnow;freezeColors;cbfreeze;
		title('Neuronal dF/F0 response during stimulus')
		ylabel('Neuron number')
		xlabel('Trial number')
		
		subplot(2,2,2)
		imagesc(matCorrelations(randperm(intPairs,100),:),[-1 1]);colormap(redblue(256));colorbar;drawnow;freezeColors;cbfreeze;
		xlabel('Stimulus repetition')
		ylabel('Neuronal pair')
		title('Neuronal pairwise correlations')
		
		
		%calculate across-block correlations
		intStimTypeNumber = length(cellSelectAll);
		intRepPairs = ((intRepetitions*intRepetitions)-intRepetitions)/2;
		vecRepCorrs = nan(1,intRepPairs);
		vecRepDists = nan(1,intRepPairs);
		vecRepTimeElapsed = nan(1,intRepPairs);
		matRepCorrs = zeros(intRepetitions,intRepetitions);
		intRepPair = 0;
		for intRep1=1:(intRepetitions-1)
			intRep1MeanTime = cellMultiSes{intPopulation}.structStim.SecsOn(round(intRep1*intStimTypeNumber-intStimTypeNumber/2));
			for intRep2=(intRep1+1):intRepetitions
				intRepPair = intRepPair + 1;
				dblCorr = corr(matCorrelations(:,intRep1),matCorrelations(:,intRep2));
				matRepCorrs(intRep1,intRep2) = dblCorr;
				matRepCorrs(intRep2,intRep1) = dblCorr;
				vecRepCorrs(intRepPair) = dblCorr;
				vecRepDists(intRepPair) = intRep2-intRep1;
				
				%get time
				intRep2MeanTime = cellMultiSes{intPopulation}.structStim.SecsOn(round(intRep2*intStimTypeNumber-intStimTypeNumber/2));
				vecRepTimeElapsed(intRepPair) = intRep2MeanTime-intRep1MeanTime;
			end
		end
		
		subplot(2,2,3)
		imagesc(matRepCorrs,[-1 1]);colormap(redblue(256));colorbar;drawnow;freezeColors;cbfreeze;
		title('Across-repetition correlations of population response structure')
		ylabel('Stimulus repetition')
		xlabel('Stimulus repetition')
		
		subplot(2,2,4)
		scatter(vecRepTimeElapsed/60,vecRepCorrs)
		ylabel('Across-repetition correlation')
		xlabel('Time between stimulus repetitions (minutes)')
		vecLimX = get(gca,'XLim');
		xlim([0 vecLimX(2)])
		
		%do regression
		sStats=regstats(vecRepCorrs,vecRepTimeElapsed/60,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStats.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		plot([0 vecLimX(2)],[0 0],'k--');
		hold off

		title(sprintf('Across-repetition correlation vs distance in time; slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%s_acrossRep_correlation%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		cellSaveAcrossRepCorr{1,intPopulation} = matCorrelations;
		cellSaveAcrossRepCorr{2,intPopulation} = matRepCorrs;
		cellSaveAcrossRepCorr{3,intPopulation} = vecRepCorrs;
		cellSaveAcrossRepCorr{4,intPopulation} = vecRepDists;
		cellSaveAcrossRepCorr{5,intPopulation} = vecRepTimeElapsed;
		cellSaveAcrossRepCorr{6,intPopulation} = sStats;
		
		%% make example correlation figure
		%{
		for intRepetition=1:size(matCorrelations,2)
			hAcrRep = figure;
			intRecording=ceil(intRepetition/2);

			strIm = [cellMultiSes{intPopulation}.strImPath(1:(end-3)) sprintf('%02d',intRecording) filesep 'average' filesep 'OverlayProc.tif'];
			imRaw = im2double(imread(strIm));
			imP = imnorm((circshift(imRaw,[1 0]) + imRaw)/2);
			imshow(imP)
			axis xy;
			hold on;
			for intPair=randperm(size(matCorrelations,1),size(matCorrelations,1))
				dblCorr = matCorrelations(intPair,intRepetition);
				if isnan(dblCorr),continue;end
				if dblCorr < 0,vecCol = [0 0 1];else vecCol = [1 0 0];end

				vecX = matLocationsX(intPair,intRepetition,:);
				vecY = matLocationsY(intPair,intRepetition,:);

				patchline(vecX,vecY,[],'EdgeColor',vecCol,'LineWidth',2,'edgealpha',abs(dblCorr)/2)
			end
			hold off;

			%save figure
			if sParams.boolSavePlots
				drawnow;
				strFig = sprintf('%s_acrossRep2_rep%d_correlation%d_raw',strSes,intRepetition,intPopulation);
				export_fig([strFig '.tif']);
				export_fig([strFig '.pdf']);
			end
		end
		%}
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['data_aggregate' strSes '_' strrep(strDate,'    ','_')];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end
%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCAENeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAE: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAE

%}