%only response until 400 ms after stim onset
for intMouse=1:8
	close all
	clearvars -except intMouse
	boolUseNeuropilSubtraction = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
			vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
		elseif intMouse == 2
			strSes = '20140314';
			vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 3
			strSes = '20140425';
			vecBlock = [1 1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 4
			strSes = '20140507';
			vecBlock = [1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 5
			strSes = '20140530';
			vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 6
			strSes = '20140604';
			vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 7
			strSes = '20140711';
			vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
		elseif intMouse == 8
			strSes = '20140715';
			vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
		end
		if boolUseNeuropilSubtraction
			strSes = ['NPS' strSes];
		end
		load(['C:\Data\Processed\StimDetectionAgg\dataPreProAggregate' strSes '.mat']);
	end
	
	%load separate ses files
	%{
	for intFile=1:numel(cellAggregate)
		%load data
		sLoad = load([cellAggregate{intFile} '.mat']);
		ses = sLoad.ses;
		
		%transform orientation to
		ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
		
		%assign to multi-cell
		cellSes{intFile} = ses;
	end
	%}
	
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
	strSes = ['NS_PostStim' strSes];
	
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		%get neuronal tuning
		if intMouse==8
			%remove last trial
			vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
			vecRem((end-47):end) = false;
			cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
			%recalc dfof
			%cellMultiSes{1} = doRecalcdFoF(cellMultiSes{1},3);
		end
		
		%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
		[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation});
		structStim = cellMultiSes{intPopulation}.structStim;
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.15;
		indTooFast = (structStim.vecTrialRespSecs-structStim.SecsOn)<dblRemSecs & structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		%structStim.vecTrialResponse(indTooFast) = 0;
		structStim = remel(structStim,~indTooFast);
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%remove trials that were too slow
		dblRemSecs = 3;
		indTooSlow = (structStim.vecTrialRespSecs-structStim.SecsOn)>dblRemSecs;
		fprintf('Removed %d trials with responses >%dms\n',sum(indTooSlow),round(dblRemSecs*1000));
		structStim.vecTrialResponse(indTooSlow) = 0;
		structStim.vecTrialRespSecs(indTooSlow) = nan;
		structStim.FrameOff = structStim.FrameOn + 10;
		
		%take opposite directions as the same
		structStim.Orientation = mod(structStim.Orientation,180);
		vecOrientations = unique(structStim.Orientation);
		cellMultiSes{intPopulation}.structStim = structStim;
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(structStim,{'Orientation'});
		cellSelectOri = getSelectionVectors(structStim,sTypesOri);
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(~indTuned) = [];
		vecNeuronPrefStim(~indTuned) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(structStim.Orientation);
		intOris = length(vecOrientations);
		%cellMultiSes{intPopulation}.structStim = structStim;
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		intContrasts = length(cellSelectContrasts);
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = structStim.vecTrialResponse==1;
		
		%% make sigmoid
		%make sigmoid plots
		cellHitActPrefPop = cell(1,intContrasts);
		cellMissActPrefPop = cell(1,intContrasts);
		cellHitActWholePop = cell(1,intContrasts);
		cellMissActWholePop = cell(1,intContrasts);
		cellHitActNPPop = cell(1,intContrasts);
		cellMissActNPPop = cell(1,intContrasts);
		
		%make normalized activity matrix per contrast
		for intContrast=1:intContrasts
			cellHitActPrefPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissActPrefPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellHitActWholePop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissActWholePop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellHitActNPPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissActNPPop{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
		end
		
		%loop through trials
		for intTrial = 1:intTrials
			intContrast = vecStimContrasts(intTrial);
			if intContrast == 1
				indNeurons = true(size(vecNeuronPrefStim));
				indNPNeurons = true(size(vecNeuronPrefStim));
			else
				indNeurons = vecNeuronPrefStim==vecStimOris(intTrial);
				indNPNeurons = vecNeuronPrefStim~=vecStimOris(intTrial);
			end
			
			vecAct = matTrialResponse(indNeurons,intTrial);
			vecActWP = matTrialResponse(:,intTrial);
			vecActNP = matTrialResponse(indNPNeurons,intTrial);
			
			if indStimResp(intTrial)
				cellHitActPrefPop{intContrast}(find(isnan(cellHitActPrefPop{intContrast}),1)) = mean(vecAct);
				cellHitActWholePop{intContrast}(find(isnan(cellHitActPrefPop{intContrast}),1)) = mean(vecActWP);
				cellHitActNPPop{intContrast}(find(isnan(cellHitActPrefPop{intContrast}),1)) = mean(vecActNP);
			else
				cellMissActPrefPop{intContrast}(find(isnan(cellMissActPrefPop{intContrast}),1)) = mean(vecAct);
				cellMissActWholePop{intContrast}(find(isnan(cellMissActPrefPop{intContrast}),1)) = mean(vecActWP);
				cellMissActNPPop{intContrast}(find(isnan(cellMissActPrefPop{intContrast}),1)) = mean(vecActNP);
			end
		end
		
		%remove trailing nans
		for intContrast=1:length(cellHitActPrefPop)
			cellHitActPrefPop{intContrast} = cellHitActPrefPop{intContrast}(~isnan(cellHitActPrefPop{intContrast}));
			cellMissActPrefPop{intContrast} = cellMissActPrefPop{intContrast}(~isnan(cellMissActPrefPop{intContrast}));
			cellHitActWholePop{intContrast} = cellHitActWholePop{intContrast}(~isnan(cellHitActWholePop{intContrast}));
			cellMissActWholePop{intContrast} = cellMissActWholePop{intContrast}(~isnan(cellMissActWholePop{intContrast}));
			cellHitActNPPop{intContrast} = cellHitActNPPop{intContrast}(~isnan(cellHitActNPPop{intContrast}));
			cellMissActNPPop{intContrast} = cellMissActNPPop{intContrast}(~isnan(cellMissActNPPop{intContrast}));
		end
		vecMeanHit = cellfun(@mean,cellHitActPrefPop);
		vecSEHit = cellfun(@std,cellHitActPrefPop)./cellfun(@numel,cellHitActPrefPop);
		vecMeanMiss = cellfun(@mean,cellMissActPrefPop);
		vecSEMiss = cellfun(@std,cellMissActPrefPop)./cellfun(@numel,cellMissActPrefPop);
		
		vecMeanHitWP = cellfun(@mean,cellHitActWholePop);
		vecSEHitWP = cellfun(@std,cellHitActWholePop)./cellfun(@numel,cellHitActWholePop);
		vecMeanMissWP = cellfun(@mean,cellMissActWholePop);
		vecSEMissWP = cellfun(@std,cellMissActWholePop)./cellfun(@numel,cellMissActWholePop);
		
		vecMeanHitNP = cellfun(@mean,cellHitActNPPop);
		vecSEHitNP = cellfun(@std,cellHitActNPPop)./cellfun(@numel,cellHitActNPPop);
		vecMeanMissNP = cellfun(@mean,cellMissActNPPop);
		vecSEMissNP = cellfun(@std,cellMissActNPPop)./cellfun(@numel,cellMissActNPPop);
		
		%save data
		cellSaveActSigmoid{intPopulation,1} = [vecMeanMiss;vecMeanHit];
		cellSaveActSigmoid{intPopulation,2} = [vecMeanMissWP;vecMeanHitWP];
		cellSaveActSigmoid{intPopulation,3} = [vecMeanMissNP;vecMeanHitNP];
		
		%plot
		vecContrasts = unique(structStim.Contrast)*100;
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		if sParams.boolSavePlots
			h=figure;
			subplot(2,2,1)
			errorfill(vecPlotC,vecMeanMiss,vecSEMiss,[1 0 0],[1 0.5 0.5]);
			hold on
			errorfill(vecPlotC,vecMeanHit,vecSEHit,[0 1 0],[0.5 1 0.5]);
			hold off
			set(gca,'XScale','log')
			set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
			ylabel('Mean pref pop dF/F0')
			xlabel('Stimulus contrast (%)')
			ylim([-0.01 0.06])
			%return
			%
			subplot(2,2,2)
			errorfill(vecPlotC,vecMeanMissWP,vecSEMissWP,[1 0 0],[1 0.5 0.5]);
			hold on
			errorfill(vecPlotC,vecMeanHitWP,vecSEHitWP,[0 1 0],[0.5 1 0.5]);
			hold off
			set(gca,'XScale','log')
			set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
			ylabel('Mean whole pop dF/F0')
			xlabel('Stimulus contrast (%)')
			ylim([-0.01 0.06])
			
			subplot(2,2,3)
			errorfill(vecPlotC,vecMeanMissNP,vecSEMissNP,[1 0 0],[1 0.5 0.5]);
			hold on
			errorfill(vecPlotC,vecMeanHitNP,vecSEHitNP,[0 1 0],[0.5 1 0.5]);
			hold off
			set(gca,'XScale','log')
			set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
			ylabel('Mean non-pref pop dF/F0')
			xlabel('Stimulus contrast (%)')
			ylim([-0.01 0.06])
			
			%save
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_dFoF_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% calculate relative increase per neuron from miss => hit trials [new heat maps]
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
		matTrialRelActHit = matTrialRelActHit(:,~isnan(matTrialRelActHit(1,:)) & ~isinf(matTrialRelActHit(1,:)));
		
		% NEW FIG
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
		
		%save fig
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_predictability1_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		% shuffle analysis
		%neurons
		intIters = 1000;
		vecSigNr = nan(1,intIters);
		vecNeuronR2 = nan(1,intIters);
		vecTrialR2 = nan(1,intIters);
		vecAllR2 = nan(1,intIters);
		for intIter=1:intIters
			matShuffledNeurons = matTrialRelActHit;
			for intTrial=1:intHitTrials
				matShuffledNeurons(:,intTrial) = matTrialRelActHit(randperm(intNeurons),intTrial);
			end
			%get nr of sign neurons
			[h,p]=ttest(matShuffledNeurons,0,0.05,'both',2);
			vecSigNr(intIter) = sum(h);
			
			%get explained variance
			matNeuronPredictorShuffled = repmat(mean(matShuffledNeurons,2),[1 intHitTrials]);
			dblVarRes=sum((matShuffledNeurons(:)-matNeuronPredictorShuffled(:)).^2);
			dblVarTot=sum((matShuffledNeurons(:)).^2);
			vecNeuronR2(intIter) = 1-dblVarRes/dblVarTot;
			
			%trials
			matShuffledTrials = matTrialRelActHit;
			matShuffledAll = matTrialRelActHit;
			for intNeuron=1:intNeurons
				matShuffledTrials(intNeuron,:) = matTrialRelActHit(intNeurons,randperm(intHitTrials));
				matShuffledAll(intNeuron,:) = matShuffledNeurons(intNeurons,randperm(intHitTrials));
			end
			
			%get explained variance
			matTrialPredictorShuffled = repmat(mean(matShuffledTrials,1),[intNeurons 1]);
			dblVarRes=sum((matShuffledTrials(:)-matTrialPredictorShuffled(:)).^2);
			dblVarTot=sum((matShuffledTrials(:)).^2);
			vecTrialR2(intIter) = 1-dblVarRes/dblVarTot;
			
			%get explained variance for both
			matNeuronAllPredictorShuffled = repmat(mean(matShuffledAll,2),[1 intHitTrials]);
			matTrialAllPredictorShuffled = repmat(mean(matShuffledAll,1),[intNeurons 1]);
			matAllPredictorShuffled = matTrialAllPredictorShuffled+matNeuronAllPredictorShuffled;
			dblVarRes=sum((matShuffledAll(:)-matAllPredictorShuffled(:)).^2);
			dblVarTot=sum((matShuffledAll(:)).^2);
			vecAllR2(intIter) = 1-dblVarRes/dblVarTot;
		end
		%test sig neurons
		[h,p]=ttest(matTrialRelActHit,0,0.05,'both',2);
		intNrSig = sum(h);
		
		%get sd's
		dblNrSDTrials = (dblPercExplainedTrials-mean(vecTrialR2))/std(vecTrialR2);
		dblNrSDAll = (dblPercExplainedAll-mean(vecAllR2))/std(vecAllR2);
		dblNrSDNeurons = (dblPercExplainedNeurons-mean(vecNeuronR2))/std(vecNeuronR2);
		dblNrSDSig = (intNrSig-mean(vecSigNr))/std(vecSigNr);
		
		if sParams.boolSavePlots
			% NEW FIG
			figure
			
			
			%CI explained variance neurons
			subplot(2,2,1)
			histx(vecNeuronR2);
			xlim([0 roundi(dblPercExplainedNeurons,2,'ceil')]);
			hold on
			plot([dblPercExplainedNeurons dblPercExplainedNeurons],get(gca,'ylim'),'r')
			hold off
			title(sprintf('R^2 by neuron ID when shuffled [real=%.3f; %.1fsd from mean shuffled]',dblPercExplainedNeurons,dblNrSDNeurons))
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
			title(sprintf('Consistently hit-modulated neurons when shuffled [real=%d;%.1fsd from mean]',intNrSig,dblNrSDSig))
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
			
			%figure
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_predictability2_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		cellSaveRelActMaps{1,intPopulation} = matTrialRelActHit;
		cellSaveRelActMaps{2,intPopulation} = matNeuronPredictor;
		cellSaveRelActMaps{3,intPopulation} = matTrialPredictor;
		cellSaveRelActMaps{4,intPopulation} = matAllPredictor;
		cellSaveRelActValsR2{intPopulation} = [dblPercExplainedNeurons dblPercExplainedNeuronsAdj dblPercExplainedTrials dblPercExplainedTrialsAdj dblPercExplainedAll intNrSig/intNeurons];
		cellSaveRelActSDR2{intPopulation} = [dblNrSDNeurons dblNrSDTrials dblNrSDAll dblNrSDSig];
		
		%% calculate noise + signal correlations FOR HIT VERSUS MISS
		%get miss/hit/slow
		indMiss = structStim.vecTrialResponse == 0;
		indHit = structStim.vecTrialResponse == 1;
		indFast = structStim.vecTrialRespSecs < nanmean(structStim.vecTrialRespSecs);
		indSlow = structStim.vecTrialRespSecs >= nanmean(structStim.vecTrialRespSecs);
		
		structStimMiss = remel(structStim,indMiss);
		structStimHit = remel(structStim,indHit);
		structStimFast = remel(structStim,indFast);
		structStimSlow = remel(structStim,indSlow);
		
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
		if sParams.boolSavePlots
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
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_HCAEneurons_correlations_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%put in output
		cellSaveStimCorrs{1,intPopulation} = matSC_Miss;
		cellSaveStimCorrs{2,intPopulation} = matNC_Miss;
		cellSaveStimCorrs{3,intPopulation} = matSC_Hit;
		cellSaveStimCorrs{4,intPopulation} = matNC_Hit;
		cellSaveStimCorrs{5,intPopulation} = matSC_Fast;
		cellSaveStimCorrs{6,intPopulation} = matNC_Fast;
		cellSaveStimCorrs{7,intPopulation} = matSC_Slow;
		cellSaveStimCorrs{8,intPopulation} = matNC_Slow;
		
		
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
		
		%normalize
		matRespNormPerNeuron = zscore(matTrialResponse,[],2);
		
		%plot
		if sParams.boolSavePlots
			figure
			colormap(hot(256))
			%imagesc(matRespNorm)
			imagesc(matTrialResponse)
			%set(gca,'XTick',1:sum(vecSelectContrastTrials):length(vecTrialContrastVector))
			colorbar
			title('Z-scored activation level per neuron per trial')
			ylabel('Neuron number')
			xlabel('Trial number, sorted by stimulus ori&contrast')
			%save plot
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%high/low DCAE
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intValuesTotal = (intNeurons*intNeurons-intNeurons)/2;
		matHeterogeneity = nan(intTrials,intValuesTotal);
		vecHeterogenPrefPop = nan(intTrials,1);
		vecHeterogenNonPrefPop = nan(intTrials,1);
		matRespDistNoRem = nan(intTrials,intValuesTotal);
		
		matTrialIndex = nan(intTrials,intValuesTotal);
		%calculate inverse correlation matrix per trial
		for intTrial=1:intTrials
			%general
			matTrialIndex(intTrial,:) = intTrial;
			
			%% STANDARD
			%perform calculation for all neurons
			vecActivity = matRespNormPerNeuron(:,intTrial);
			matDistAll = abs(bsxfun(@minus,vecActivity,vecActivity'));
			
			%save data as vectors in matrix
			matSelectAll = tril(true(size(matDistAll)),-1);
			vecRespDistAll = matDistAll(matSelectAll);
			
			matRespDistNoRem(intTrial,:) = vecRespDistAll;
			matHeterogeneity(intTrial,:) = vecRespDistAll;
			
			%% only pref pop
			indNeurons = vecNeuronPrefStim==vecStimOris(intTrial);
			vecActivity = matRespNormPerNeuron(indNeurons,intTrial);
			matDistPP=abs(bsxfun(@minus,vecActivity,vecActivity'));
			
			%save data as vectors in matrix
			matSelectPP = tril(true(size(matDistPP)),-1);
			vecRespDistPP = matDistPP(matSelectPP);
			
			vecHeterogenPrefPop(intTrial) = mean(vecRespDistPP);
			
			%% only non-pref pop
			indNeurons = vecNeuronPrefStim~=vecStimOris(intTrial);
			vecActivity = matRespNormPerNeuron(indNeurons,intTrial);
			matDistNPP=abs(bsxfun(@minus,vecActivity,vecActivity'));
			
			%save data as vectors in matrix
			matSelectNPP = tril(true(size(matDistNPP)),-1);
			vecRespDistNPP = matDistNPP(matSelectNPP);
			
			vecHeterogenNonPrefPop(intTrial) = mean(vecRespDistNPP);
		end
		
		%no rem
		matRespDist = matRespDistNoRem;
		%split data set
		indSelectContrasts = structStim.Contrast == 0.005 | structStim.Contrast == 0.02  | structStim.Contrast == 0.08  | structStim.Contrast == 0.32;
		indSelectResp = structStim.vecTrialResponse;
		indSelectNoRespTrials = indSelectContrasts & ~indSelectResp;
		
		matNoResp = matRespDist(indSelectNoRespTrials,:);
		
		%hit trials only
		indSelectRespTrials = indSelectContrasts & indSelectResp;
		matRespDistSub = matRespDist(indSelectRespTrials,:);
		
		%put in output
		cellSaveNormActDissim{1,1,intPopulation} = matRespDistSub(:); %within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals) with every cell = vector of dissimilarity values
		cellSaveNormActDissim{1,2,intPopulation} = matNoResp(:);
		
		%no rem
		matRespDist = matRespDistNoRem;
		
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
		
		%source matrices
		%matTrialResponse(intNeuron,intTrial);
		%matRespNormPerContrastintNeuron,intTrial);
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		%for whole pop or only preferred pop
		hPopCorr = figure;
		for intPopType=1:2
			if intPopType == 1
				strPopType = 'PrefNeurons';
			elseif intPopType == 2
				strPopType = 'NonPrefNeurons';
			end
			for intPrefType=1:4
				indSelect = vecNeuronPrefStim' == intPrefType;
				if intPopType == 2
					indSelect = ~indSelect;
				end
				%strPopType
				%sum(indSelect)
				
				%calculate trial-by-trial correlations of population activity
				if sum(indSelect) > 0
					cellCorr{intPrefType}=corr(matTrialResponse(indSelect,:));
				else
					cellCorr{intPrefType} = nan(intTrials,intTrials);
				end
			end
			matCorr = nan(size(cellCorr{1}));
			for intOriType=vecOriType
				if isnan(intOriType),continue;end
				%get trials
				indSelectTrials = vecOriType == intOriType;
				%assign
				matCorr(indSelectTrials,indSelectTrials) = cellCorr{intPrefType}(indSelectTrials,indSelectTrials);
			end
			
			
			%compute absolute trial distance
			matTrialDist = abs(repmat(1:intTrials,[intTrials 1])-repmat((1:intTrials)',[1 intTrials]));
			%matSelect = tril(true(size(matTrialDist)),-1);
			indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
			indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
			indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs <= nanmedian(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs > nanmedian(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			
			%separate across-trial correlation correlations for different trial types (1-4)
			%also split for hit/miss trials
			intHitTrials = sum(indHit);
			intMissTrials = sum(indMiss);
			intFastTrials = sum(indFast);
			intSlowTrials = sum(indSlow);
			vecCorrHit = [];
			vecCorrMiss = [];
			vecCorrFast = [];
			vecCorrSlow = [];
			for intOriType=vecOriType
				%% hit; select trials
				if isnan(intOriType),continue;end
				indSelectHitTrials = vecOriType == intOriType & indHit;
				if sum(indSelectHitTrials) > 2
					%get trials
					matCorrSubHit = matCorr(indSelectHitTrials,indSelectHitTrials);
					
					%subselection matrix hit
					matSubSelectHit = tril(true(size(matCorrSubHit)),-1);
					
					%subselect
					vecCorrHit = [vecCorrHit matCorrSubHit(matSubSelectHit)'];
				end
				
				%% miss; select trials
				indSelectMissTrials = vecOriType == intOriType & indMiss;
				if sum(indSelectMissTrials) > 2
					%get trials
					matCorrSubMiss = matCorr(indSelectMissTrials,indSelectMissTrials);
					
					%subselection matrix hit
					matSubSelectMiss = tril(true(size(matCorrSubMiss)),-1);
					
					%subselect
					vecCorrMiss = [vecCorrMiss matCorrSubMiss(matSubSelectMiss)'];
				end
				
				%% fast; select trials
				indSelectFastTrials = vecOriType == intOriType & indFast;
				if sum(indSelectFastTrials) > 2
					%get trials
					matCorrSubFast = matCorr(indSelectFastTrials,indSelectFastTrials);
					
					%subselection matrix hit
					matSubSelectFast = tril(true(size(matCorrSubFast)),-1);
					
					%subselect
					vecCorrFast = [vecCorrFast matCorrSubFast(matSubSelectFast)'];
				end
				
				%% slow; select trials
				indSelectSlowTrials = vecOriType == intOriType & indSlow;
				if sum(indSelectSlowTrials) > 2
					%get trials
					matCorrSubSlow = matCorr(indSelectSlowTrials,indSelectSlowTrials);
					
					%subselection matrix hit
					matSubSelectSlow = tril(true(size(matCorrSubSlow)),-1);
					
					%subselect
					vecCorrSlow = [vecCorrSlow matCorrSubSlow(matSubSelectSlow)'];
				end
			end
			
			%plot
			if sParams.boolSavePlots
				subplot(2,2,1+(intPopType-1)*2)
				matPlotCorr = matCorr(~isnan(vecOriType),~isnan(vecOriType));
				matPlotCorr(diag(true(1,length(matPlotCorr)))) = nan;
				vecLim = [-nanmax(abs(matPlotCorr(:))) nanmax(abs(matPlotCorr(:)))];
				imagesc(matPlotCorr,vecLim);
				matC = redbluepurple(128);
				nancolorbar(matPlotCorr,vecLim,matC);
				title(sprintf('Inter-trial correlation of %s population activity (block %d)',strPopType,intPopulation))
				xlabel('Trial')
				ylabel('Trial')
				
				subplot(2,2,2+(intPopType-1)*2)
				dblMeanCorrHit = mean(vecCorrHit);
				dblMeanCorrMiss = mean(vecCorrMiss);
				dblMeanCorrFast = mean(vecCorrFast);
				dblMeanCorrSlow = mean(vecCorrSlow);
				dblErrCorrHit = std(vecCorrHit)/sqrt(length(vecCorrHit));
				dblErrCorrMiss = mean(vecCorrMiss)/sqrt(length(vecCorrMiss));
				dblErrCorrFast = mean(vecCorrFast)/sqrt(length(vecCorrFast));
				dblErrCorrSlow = mean(vecCorrSlow)/sqrt(length(vecCorrSlow));
				
				if intPopType==1 %pref
					vecMeanCorrPref = [dblMeanCorrMiss dblMeanCorrSlow dblMeanCorrFast];
				else
					vecMeanCorrNonPref = [dblMeanCorrMiss dblMeanCorrSlow dblMeanCorrFast];
				end
				[h,pHM] = ttest2(vecCorrHit,vecCorrMiss);
				[h,pMF] = ttest2(vecCorrMiss,vecCorrFast);
				[h,pMS] = ttest2(vecCorrMiss,vecCorrSlow);
				[h,pFS] = ttest2(vecCorrFast,vecCorrSlow);
				
				errorbar(1:4,[dblMeanCorrMiss dblMeanCorrHit dblMeanCorrFast dblMeanCorrSlow],[dblErrCorrMiss dblErrCorrHit dblErrCorrFast dblErrCorrSlow],'Linestyle','none','Marker','x','Color','k');
				set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
				title(sprintf('Mean inter-trial correlation of %s population activity\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
				ylabel('Pearson correlation of population activity')
			end
		end
		%%%% OUTPUT
		cellSaveITC{intPopulation} = [vecMeanCorrNonPref;vecMeanCorrPref];
		
		%plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(hPopCorr), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopCorr);
			drawnow;
			strFig = sprintf('%sagg_intertrialcorr_%s_pop%d__raw',strSes,strPopType,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%{
		%% inter-trial correlations
		%for whole pop or only preferred pop
		matRespNormPerTrial = zscore(matTrialResponse,[],1);
		matCorrAll = nan(intTrials,intTrials);
		matCorrPref = nan(intTrials,intTrials);
		matCorrNonPref = nan(intTrials,intTrials);
		for intContrast=2:(intContrasts-1)
			indContrastTrials = vecStimContrasts==intContrast;
			for intOri=1:intOris
				%get general info about this ori
				indOriTrials = vecStimOris==intOri;
				indPrefPop = vecNeuronPrefStim == intOri;
				indNonPrefPop = vecNeuronPrefStim ~= intOri;
				indSelectTrials = indOriTrials & indContrastTrials;
				intOriConReps = sum(indSelectTrials);
				vecOriConTrials = find(indSelectTrials);
				
				%all neurons
				matCorrOri = nan(intOriConReps,intOriConReps);
				for intTrialCounter=1:intOriConReps
					intOriConTrial = vecOriConTrials(intTrialCounter);
					matCorrOri(intTrialCounter,:) = mean(bsxfun(@times,matRespNormPerTrial(:,intOriConTrial),matRespNormPerTrial(:,vecOriConTrials))) * (intOriConReps/(intOriConReps-1));
				end
				matCorrAll(vecOriConTrials,vecOriConTrials) = matCorrOri;
				
				%pref neurons
				matCorrOri = nan(intOriConReps,intOriConReps);
				for intTrialCounter=1:intOriConReps
					intOriConTrial = vecOriConTrials(intTrialCounter);
					matCorrOri(intTrialCounter,:) = mean(bsxfun(@times,matRespNormPerTrial(indPrefPop,intOriConTrial),matRespNormPerTrial(indPrefPop,vecOriConTrials))) * (intOriConReps/(intOriConReps-1));
				end
				matCorrPref(vecOriConTrials,vecOriConTrials) = matCorrOri;
				
				%pref neurons
				matCorrOri = nan(intOriConReps,intOriConReps);
				for intTrialCounter=1:intOriConReps
					intOriConTrial = vecOriConTrials(intTrialCounter);
					matCorrOri(intTrialCounter,:) = mean(bsxfun(@times,matRespNormPerTrial(indNonPrefPop,intOriConTrial),matRespNormPerTrial(indNonPrefPop,vecOriConTrials))) * (intOriConReps/(intOriConReps-1));
				end
				matCorrNonPref(vecOriConTrials,vecOriConTrials) = matCorrOri;
			end
		end
		hPopCorr = figure;
		subplot(2,2,1)
		imagesc(matCorrAll,[-1 1]);colormap(redblue);nancolorbar(matCorrAll,[-1 1],redblue);
		title('Whole population');
		xlabel('Trial')
		ylabel('Trial')
		
		subplot(2,2,3)
		imagesc(matCorrNonPref,[-1 1]);colormap(redblue);nancolorbar(matCorrNonPref,[-1 1],redblue);
		title('Non-preferred population');
		xlabel('Trial')
		ylabel('Trial')
		
		subplot(2,2,4)
		imagesc(matCorrPref,[-1 1]);colormap(redblue);nancolorbar(matCorrPref,[-1 1],redblue);
		title('Preferred population');
		xlabel('Trial')
		ylabel('Trial')
		
		%split for miss/slow/fast
		vecMiss = find(indMiss);
		vecFast = find(indFast);
		vecSlow = find(indSlow);
		vecMeanCorrPref = nan(1,3);
		vecErCorrPref = nan(1,3);
		vecMeanCorrNonPref = nan(1,3);
		vecErCorrNonPref = nan(1,3);
		for intResp = 1:3
			%get data
			if intResp==1,vecResp=vecMiss;
			elseif intResp==2,vecResp=vecSlow;
			elseif intResp==3,vecResp=vecFast;
			end
			%calc
			matThisCorrP = matCorrPref(vecResp,vecResp);
			matThisCorrNP = matCorrNonPref(vecResp,vecResp);
			matSelect = tril(true(size(matThisCorrP)),-1);
			vecCorrP = matThisCorrP(matSelect);
			vecCorrNP = matThisCorrNP(matSelect);
			dblErrF = sum(matSelect(:));
			%dblErrF = length(vecResp)^2;
			
			%cellRespP{intResp} = [cellRespP{intResp};vecCorrP(~isnan(vecCorrP))];
			%cellRespNP{intResp} = [cellRespP{intResp};vecCorrNP(~isnan(vecCorrNP))];
			
			%save
			vecMeanCorrPref(intResp) = nanmean(vecCorrP(:));
			vecErCorrPref(intResp) = nanstd(vecCorrP(:))/sqrt(dblErrF);
			vecMeanCorrNonPref(intResp) = nanmean(vecCorrNP(:));
			vecErCorrNonPref(intResp) = nanstd(vecCorrNP(:))/sqrt(dblErrF);
		end
		hPopCorrPlot = figure;
		subplot(2,2,1)
		cellColor={'r','m','g'};
		hold on;
		for intResp=1:3
			errorbar(intResp,vecMeanCorrNonPref(intResp),vecErCorrNonPref(intResp),['x' cellColor{intResp}]);
		end
		hold off
		title('Non-preferred population');
		cellLabels = {'Miss','Slow','Fast'};
		set(gca,'xtick',1:3,'xticklabel',cellLabels)
		ylabel('Inter-trial correlation');
		ylim([0 1]);
		
		subplot(2,2,2)
		hold on;
		for intResp=1:3
			errorbar(intResp,vecMeanCorrPref(intResp),vecErCorrPref(intResp),['x' cellColor{intResp}]);
		end
		hold off
		title('Preferred population');
		set(gca,'xtick',1:3,'xticklabel',cellLabels)
		ylabel('Inter-trial correlation');
		ylim([0 1]);
		
		
		%save data
		cellSaveITC{intPopulation} = [vecMeanCorrNonPref;vecMeanCorrPref];
		
		%% plot
		if sParams.boolSavePlots
			figure(hPopCorr);
			drawnow;
			jFig = get(handle(hPopCorr), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopCorr);
			drawnow;
			strFig = sprintf('%sagg_intertrialcorrmaps_pop%d__raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
			
			figure(hPopCorrPlot);
			drawnow;
			jFig = get(handle(hPopCorrPlot), 'JavaFrame');
			jFig.setMaximized(true);
			figure(hPopCorrPlot);
			drawnow;
			strFig = sprintf('%sagg_intertrialcorrplot_pop%d__raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%}
		
		%% plot dependency on reaction times
		vecRTs = structStim.vecTrialRespSecs(indSelectRespTrials) - structStim.SecsOn(indSelectRespTrials); %get RTs for hit trials and selected contrasts
		[vecRTsSorted,vecRTSortIndex] = sort(vecRTs,'ascend');
		matHitRTSorted = matRespDistSub(vecRTSortIndex,:);
		vecHitRT = nanmean(matHitRTSorted,2);
		
		%Z activity
		matActZ = matRespNormPerNeuron(:,indSelectRespTrials)';
		
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
		if sParams.boolSavePlots
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
		sParamsHet.intWindowLength = 1;
		%[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		
		%pre-allocate
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrasts(1) = 0.2;
		intContrasts = length(vecContrasts);
		vecP = nan(1,intContrasts);
		vecMetaHitY = nan(1,intContrasts);
		vecMetaHitE = nan(1,intContrasts);
		vecMetaMissY = nan(1,intContrasts);
		vecMetaMissE = nan(1,intContrasts);
		vecPPP = nan(1,intContrasts);
		vecMetaHitPPY = nan(1,intContrasts);
		vecMetaHitPPE = nan(1,intContrasts);
		vecMetaMissPPY = nan(1,intContrasts);
		vecMetaMissPPE = nan(1,intContrasts);
		vecPNPP = nan(1,intContrasts);
		vecMetaHitNPPY = nan(1,intContrasts);
		vecMetaHitNPPE = nan(1,intContrasts);
		vecMetaMissNPPY = nan(1,intContrasts);
		vecMetaMissNPPE = nan(1,intContrasts);
		vecActCohensD = nan(1,intContrasts);
		vecActPPCohensD = nan(1,intContrasts);
		vecActNPCohensD = nan(1,intContrasts);
		vecHetCohensD = nan(1,intContrasts);
		vecHetPPCohensD = nan(1,intContrasts);
		vecHetNPPCohensD = nan(1,intContrasts);
		for intContrastIndex=1:intContrasts
			%get contrast
			dblContrast = vecContrasts(intContrastIndex);
			
			%get resp
			vecHeteroHit = mean(matHeterogeneity(cellSelectContrasts{intContrastIndex} & indSelectResp,:),2);
			%vecHeteroHit = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & indSelectResp);
			vecHeteroHit = vecHeteroHit(:);
			vecHeteroMiss = mean(matHeterogeneity(cellSelectContrasts{intContrastIndex} & ~indSelectResp,:),2);
			%vecHeteroMiss = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & ~indSelectResp);
			vecHeteroMiss = vecHeteroMiss(:);
			
			%get pref resp
			vecHeteroHitPP = vecHeterogenPrefPop(cellSelectContrasts{intContrastIndex} & indSelectResp);
			vecHeteroHitPP = vecHeteroHitPP(:);
			vecHeteroMissPP = vecHeterogenPrefPop(cellSelectContrasts{intContrastIndex} & ~indSelectResp);
			vecHeteroMissPP = vecHeteroMissPP(:);
			
			%get pref resp
			vecHeteroHitNPP = vecHeterogenNonPrefPop(cellSelectContrasts{intContrastIndex} & indSelectResp);
			vecHeteroHitNPP = vecHeteroHitNPP(:);
			vecHeteroMissNPP = vecHeterogenNonPrefPop(cellSelectContrasts{intContrastIndex} & ~indSelectResp);
			vecHeteroMissNPP = vecHeteroMissNPP(:);
			
			%perform analyses per contrast
			[h,p,ci] = ttest2(vecHeteroHit,vecHeteroMiss);
			vecY = [mean(vecHeteroHit) mean(vecHeteroMiss)];
			[h,pPP,ci] = ttest2(vecHeteroHitPP,vecHeteroMissPP);
			vecYPP = [mean(vecHeteroHitPP) mean(vecHeteroMissPP)];
			[h,pNPP,ci] = ttest2(vecHeteroHitNPP,vecHeteroMissNPP);
			vecYNPP = [mean(vecHeteroHitNPP) mean(vecHeteroMissNPP)];
			
			%put in output
			cellSaveMatContHetero{intPopulation}(intContrastIndex,1) = vecY(1);
			cellSaveMatContHetero{intPopulation}(intContrastIndex,2) = vecY(2);
			cellSaveMatContHetero{intPopulation}(intContrastIndex,3) = vecYPP(1); %pref pop
			cellSaveMatContHetero{intPopulation}(intContrastIndex,4) = vecYPP(2); %pref pop
			cellSaveMatContHetero{intPopulation}(intContrastIndex,5) = vecYNPP(1); %non-pref pop
			cellSaveMatContHetero{intPopulation}(intContrastIndex,6) = vecYNPP(2); %non-pref pop
			
			%put in meta vector
			vecP(intContrastIndex) = p;
			vecMetaHitY(intContrastIndex) = vecY(1);
			vecMetaHitE(intContrastIndex) = std(vecHeteroHit)/sqrt(intTrials);
			vecMetaMissY(intContrastIndex) = vecY(2);
			vecMetaMissE(intContrastIndex) = std(vecHeteroMissPP)/sqrt(intTrials);
			vecPPP(intContrastIndex) = pPP;
			vecMetaHitPPY(intContrastIndex) = vecYPP(1);
			vecMetaHitPPE(intContrastIndex) = std(vecHeteroHitPP)/sqrt(intTrials);
			vecMetaMissPPY(intContrastIndex) = vecYPP(2);
			vecMetaMissPPE(intContrastIndex) = std(vecHeteroMissPP)/sqrt(intTrials);
			vecPNPP(intContrastIndex) = pNPP;
			vecMetaHitNPPY(intContrastIndex) = vecYNPP(1);
			vecMetaHitNPPE(intContrastIndex) = std(vecHeteroHitNPP)/sqrt(intTrials);
			vecMetaMissNPPY(intContrastIndex) = vecYNPP(2);
			vecMetaMissNPPE(intContrastIndex) = std(vecHeteroMissNPP)/sqrt(intTrials);
			
			%cohen's d
			vecActCohensD(intContrastIndex) = getCohensD(cellHitActWholePop{intContrastIndex},cellMissActWholePop{intContrastIndex});
			vecActPPCohensD(intContrastIndex) = getCohensD(cellHitActPrefPop{intContrastIndex},cellMissActPrefPop{intContrastIndex});
			vecActNPCohensD(intContrastIndex) = getCohensD(cellHitActNPPop{intContrastIndex},cellMissActNPPop{intContrastIndex});
			
			vecHetCohensD(intContrastIndex) = getCohensD(vecHeteroHit,vecHeteroMiss);
			vecHetPPCohensD(intContrastIndex) = getCohensD(vecHeteroHitPP,vecHeteroMissPP);
			vecHetNPPCohensD(intContrastIndex) = getCohensD(vecHeteroHitNPP,vecHeteroMissNPP);
		end
		%save cohen's d
		cellSaveActCohensD{intPopulation} = vecActCohensD;
		cellSaveActPPCohensD{intPopulation} = vecActPPCohensD;
		cellSaveActNPPCohensD{intPopulation} = vecActNPCohensD;
		cellSaveHetCohensD{intPopulation} = vecHetCohensD;
		cellSaveHetPPCohensD{intPopulation} = vecHetPPCohensD;
		cellSaveHetNPPCohensD{intPopulation} = vecHetNPPCohensD;
		
		%pre-compute variables
		%if sParams.boolSavePlots
		%figure
		hHetCon = figure;
		set(hHetCon,'Color',[1 1 1]);
		figure(hHetCon);
		
		subplot(2,2,1)
		vecWindow = [1 length(vecContrasts)];
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		
		%plot heterogeneity
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
		
		%plot pref-pop heterogeneity
		subplot(2,2,2)
		for intResp=[0 1]
			if intResp == 1
				vecMeanTrace = vecMetaHitPPY(vecWindowSelect);
				vecSE = vecMetaHitPPE(vecWindowSelect);
				vecColorFill = [0.7 1 0.7];
				vecColorLine = [0 1 0];
			else
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
				vecMeanTrace = vecMetaMissPPY(vecWindowSelect);
				vecSE = vecMetaMissPPE(vecWindowSelect);
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
		title(sprintf('Pref-Pop het [%d]; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',intPopulation,vecPPP))
		grid on
		xlabel('Contrast')
		ylabel('Pref-Pop response heterogeneity')
		xlim(vecContrasts(vecWindow))
		%ylim([-0.01 0.06])
		legend({'SEM','Miss','SEM','Hit'},'Location','Best')
		
		%plot non-pref-pop heterogeneity
		subplot(2,2,3)
		for intResp=[0 1]
			if intResp == 1
				vecMeanTrace = vecMetaHitNPPY(vecWindowSelect);
				vecSE = vecMetaHitNPPE(vecWindowSelect);
				vecColorFill = [0.7 1 0.7];
				vecColorLine = [0 1 0];
			else
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
				vecMeanTrace = vecMetaMissNPPY(vecWindowSelect);
				vecSE = vecMetaMissNPPE(vecWindowSelect);
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
		title(sprintf('Non-pref-Pop het [%d]; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',intPopulation,vecPNPP))
		grid on
		xlabel('Contrast')
		ylabel('Non pref-Pop response heterogeneity')
		xlim(vecContrasts(vecWindow))
		%ylim([-0.01 0.06])
		legend({'SEM','Miss','SEM','Hit'},'Location','Best')
		
		%save fig
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%s_heterogeneity_over_contrasts_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
		%end
		
		%% calculate Cohen's d vs. behavioural d' to check dependency of the effects
		%[phat,pci] = binofit(x,n)
		
		%% dF/F over time + heterogeneity over time; also split for fast/slow trials
		sParamsHet.intWindowLength = 1;
		[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		
		%get stim data
		cellFieldsC = {'Contrast'};
		sTypesC = getStimulusTypes(structStim,cellFieldsC);
		cellSelectC = getSelectionVectors(structStim,sTypesC);
		vecC = sTypesC.matTypes;
		
		%get resp data
		indMiss = structStim.vecTrialResponse == 0;
		indFast = structStim.vecTrialRespSecs < nanmean(structStim.vecTrialRespSecs);
		indSlow = structStim.vecTrialRespSecs >= nanmean(structStim.vecTrialRespSecs);
		
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
					vecStimOn = structStim.FrameOn(vecTrials);
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
		indMiss = isnan(structStim.vecTrialRespSecs);
		indFast = structStim.vecTrialRespSecs < nanmean(structStim.vecTrialRespSecs);
		indSlow = structStim.vecTrialRespSecs >= nanmean(structStim.vecTrialRespSecs);
		
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
		vecStimOn = structStim.FrameOn;
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
		if sParams.boolSavePlots
			%fig
			hDecRespTypeTime = figure;
		end
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
			if sParams.boolSavePlots
				subplot(2,2,intRespType)
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
		end
		
		%plot summary
		if sParams.boolSavePlots
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
			drawnow;
			strFig = sprintf('%s_resptypedecoding_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% plot heterogeneity vs correlation ideal/actual pop response (pref stim)
		%matTrialResponse(intNeuron,intTrial)
		[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matTrialResponse);
		vecOriTrial = structStim.Orientation;
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
		sStats=regstats(vecPrefRespCorr,vecHeterogeneity,'linear',{'beta','rsquare','tstat'});
		if sParams.boolSavePlots
			hHetPopCorr = figure;
			scatter(vecHeterogeneity,vecPrefRespCorr);
			
			vecX = get(gca,'XLim');
			vecY = polyval(sStats.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			
			title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
			
			xlabel('Heterogeneity')
			ylabel('Correlation ideal/actual population response')
		end
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