%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
for intMouse=1:8%3%5%7
	close all
	clearvars -except intMouse
	boolUseNeuropilSubtraction = true;
	
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
		load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat']);
	end
	
	%load separate ses files
	vecRecordings = [];
	for intFile=1:numel(cellAggregate)
		%define data location
		if boolUseNeuropilSubtraction
			strSes2 = strSes(4:end);
		end
		strFile = cellAggregate{intFile};
		vecLast = find(strFile==filesep,2,'last');
		intRec = str2double(getFlankedBy(cellAggregate{intFile},[strSes2 'xyt'],'_'));
		if ~strcmp(strFile((vecLast(1)+1):(vecLast(2)-3)),'xyt')
			strFile = [strFile(1:vecLast(2)) sprintf('xyt%02d',intRec) filesep strFile((vecLast(2)+1):end)];
		end
		
		%load data
		sLoad = load([strFile '.mat']);
		ses = sLoad.ses;
		
		%transform orientation to
		ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
		
		%assign to multi-cell
		cellSes{intFile} = ses;
		%xyt01
		%get recording
		vecRecordings = [vecRecordings intRec];
	end
	
	
	%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	vecBlockTypes = unique(vecBlock);
	intNumBlocks = length(vecBlockTypes);
	%vecNeuronNum = zeros(1,intNumBlocks);
	%cellKeepList = cell(1,intNumBlocks);
	%#ok<*ASGLU>
	%#ok<*AGROW>
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%% hit rate + RT over time
	%put in structure
	if boolUseNeuropilSubtraction
		sIn.strSes = strSes(4:end);
	else
		sIn.strSes = strSes;
	end
	sIn.strMasterPath = 'D:\Data\Processed\imagingdata';
	sIn.vecRecordings = vecRecordings;
	%%{
	fprintf('Building aggregate stimulus structure...\n')
	sStimAggregate = buildStimAggregate(sIn);
	
	%change name for no split
	strSes = ['NS_Supp' strSes];
	
	%get data
	sStimAggregate.Orientation = mod(sStimAggregate.Orientation,180);
	sTypesAll = getStimulusTypes(sStimAggregate,{'Orientation','Contrast'});
	cellSelectAll = getSelectionVectors(sStimAggregate,sTypesAll);
	matStimResponse = doMatRespTransform(sStimAggregate.vecTrialResponse,cellSelectAll);
	vecRT = sStimAggregate.vecTrialRespSecs-sStimAggregate.SecsOn;
	vecRT(sStimAggregate.vecTrialResponse==0) = nan;
	matStimRespTime = doMatRespTransform(vecRT,cellSelectAll);
	
	%% calculate psychometric curve on % responded and RT
	fprintf('Calculating psychometric curve...\n')
	sPC = struct;
	sPC = calcPsychometricCurve(sStimAggregate);
	%matChiSquareP = sPC.matP;
	cellSaveBehavDetect = shiftdim(sPC.cellResponses,1); %behavioral detection [6 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
	cellSaveBehavRT = cell(6,1); %behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
	for intC = 1:6
		vecRT = sPC.matRTs(:,intC);
		cellSaveBehavRT{intC} = vecRT(~isnan(vecRT));
	end

	%resp over C
	matRTC = nan(6,size(matStimResponse,2));
	matRespC = nan(6,size(matStimResponse,2));
	for intC=1:6
		intStart=(intC-1)*4+1;
		intStop=intC*4;
		matRespC(intC,:) = mean(matStimResponse(intStart:intStop,:));
		matRTC(intC,:) = nanmean(matStimRespTime(intStart:intStop,:));
	end
	
	%resp over O
	matRespO = nan(4,size(matStimResponse,2));
	matRTO = nan(4,size(matStimResponse,2));
	for intO=1:4
		intStart=(intO-1)*6+1;
		intStop=intO*6;
		matRespO(intO,:) = mean(matStimResponse(intStart:intStop,:));
		matRTO(intO,:) = nanmean(matStimRespTime(intStart:intStop,:));
	end
	
	%plot
	hBehavOverTime = figure;
	
	subplot(2,2,1)
	plot(matRespC')
	hold on
	scatter(1:size(matRespC,2),nanmean(matRespC),[],[0 0 0],'filled')
	hold off
	ylabel('Response Proportion')
	xlabel('Repetition Block')
	title([strSes '; Response proportion per contrast (colors)'])
	
	subplot(2,2,2)
	plot(matRespO')
	hold on
	scatter(1:size(matRespO,2),nanmean(matRespO),[],[0 0 0],'filled')
	hold off
	ylabel('Response Proportion')
	xlabel('Repetition Block')
	title('Response proportion per orientation (colors)')
	
	subplot(2,2,3)
	plot(matRTC')
	hold on
	scatter(1:size(matRTC,2),nanmean(matRTC),[],[0 0 0],'filled')
	hold off
	ylabel('Reaction time')
	xlabel('Repetition Block')
	title('Reaction time per contrast (colors)')
	ylim([0 3])
	
	subplot(2,2,4)
	plot(matRTO')
	hold on
	scatter(1:size(matRTO,2),nanmean(matRTO),[],[0 0 0],'filled')
	hold off
	ylabel('Reaction time')
	xlabel('Repetition Block')
	title('Reaction time per orientation (colors)')
	ylim([0 3])
	
	%save data
	cellSaveHitOverTime{1} = matRespC;
	cellSaveHitOverTime{2} = matRTC;
	cellSaveHitOverTime{3} = matRespO;
	cellSaveHitOverTime{4} = matRTO;
	
	%save figure
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%s_behavOverTime%d_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	%}
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
		%structStim.FrameOff = structStim.FrameOn+1;
		
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
		
		s=semblance(1:intTrials,matTrialResponse(1,:),matTrialResponse(2,:),intTrials)
		
		
		%normalize per contrast
		%matRespNormPerContrast = nan(size(matTrialResponse));
		%for intContrastIndex=1:length(cellSelectContrasts)
		%	vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
		%	matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		%end
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		
		%% calculate response to preferred stimulus over time
		%calculate response to preferred orientation over time
		matSelectC = repmat(cellSelectContrasts{end},[intNeurons 1]);
		vecLookupO = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		matPrefO = vecLookupO(repmat(vecNeuronPrefStim',[1 intTrials]));
		matStimO = repmat(cellMultiSes{intPopulation}.structStim.Orientation,[intNeurons 1]);
		matSelectO = matPrefO == matStimO;
		matSelect = matSelectC & matSelectO;
		vecReps = sum(matSelect,2);
		%{
		if any(vecReps ~= vecReps(1)),error([mfilename ':UnevenRepetitions'],'Number of repetitions unequal; %s pop%d',strSes,intPopulation);end
		matPrefResp = reshape(matTrialResponse(matSelect),size(matSelect,1),vecReps(1));
		
		vecTrialOri = nan(1,intTrials);
		for intOri=1:length(vecLookupO)
			vecTrialOri(vecLookupO(intOri)==cellMultiSes{intPopulation}.structStim.Orientation) = intOri;
		end
		
		%calculate preferred orientation over time
		vecOrisFullC = cellMultiSes{intPopulation}.structStim.Orientation(cellSelectContrasts{end});
		matRespFullC = matTrialResponse(:,cellSelectContrasts{end});
		intStepSize = 8;
		intTotalSteps = size(matRespFullC,2)/intStepSize;
		matPrefOri = nan(intNeurons,intTotalSteps);
		for intStep=1:intTotalSteps
			intStart = (intStep-1)*intStepSize+1;
			vecSelect = intStart:(intStart+intStepSize-1);
			vecOri = vecOrisFullC(vecSelect);
			[dummy,vecSort] = sort(vecOri,'ascend');
			matThisResp = matRespFullC(:,vecSelect(vecSort));
			matRepByOri = reshape(matThisResp,intNeurons,2,4);
			matOriResp = squeeze(mean(matRepByOri,2));
			[dummy,vecPrefOri] = max(matOriResp,[],2);
			matPrefOri(:,intStep) = vecPrefOri;
		end
		
		%calculate max resp per repetition
		intRepSize = 48;
		intReps = size(matTrialResponse,2)/intRepSize;
		matMaxResp = nan(intNeurons,intReps);
		for intRep=1:intReps
			intStart = (intRep-1)*intRepSize+1;
			vecSelect = intStart:(intStart+intRepSize-1);
			matMaxResp(:,intRep) = max(matTrialResponse(:,vecSelect),[],2);
		end
		matNormMaxResp = matMaxResp ./ repmat(mean(matMaxResp,2),[1 intReps]);
		
		%%{
		%calculate relative fluorescence of soma vs neuropil
		intDownsample = 25;
		matRatioF = nan(intNeurons,ceil(length(cellMultiSes{intPopulation}.neuron(1).F)/intDownsample));
		vecFilter = buildGaussian(-2:0.5:2,0,1,0);
		vecMaxEpochDurBelowThresh = nan(1,intNeurons);
		for intNeuron=1:intNeurons
			vecRatio = cellMultiSes{intPopulation}.neuron(intNeuron).F ./ (cellMultiSes{intPopulation}.neuron(intNeuron).F+cellMultiSes{intPopulation}.neuron(intNeuron).npF);
			vecBelowThresh = [false diff(vecRatio > 0.5)];
			vecUp = find(vecBelowThresh==1);
			vecDown = find(vecBelowThresh==-1);
			if isempty(vecUp) || isempty(vecDown),vecUp=0;vecDown=0;
			elseif length(vecDown) < length(vecUp),vecDown(end+1) = vecUp(end);
			elseif length(vecUp) < length(vecDown),vecUp(end+1) = length(vecRatio);
			end
			
			vecMaxEpochDurBelowThresh(intNeuron) = max(vecUp-vecDown);
			matRatioF(intNeuron,:) = conv(decimate(vecRatio,intDownsample)-0.5,vecFilter,'same')+0.5;
		end
		
		%plot traces
		h=figure;
		intPlottedNeurons = 10;
		vecPlottedNeurons = randperm(intNeurons,intPlottedNeurons);
		vecPlottedNeurons = sort(vecPlottedNeurons);
		hold on
		plot(matRatioF(vecPlottedNeurons,5:(end-5))');
		plot(get(gca,'XLim'),[0.5 0.5],'k--')
		hold off;
		ylim([0.35 0.65+eps]);
		ylabel('Ratio fluorescence Fsoma/(Fsoma+Fneuropil)')
		xlabel('Time (s)')
		title(sprintf('%s pop%d',strSes,intPopulation))
		legend(arrayfun(@num2str,vecPlottedNeurons', 'UniformOutput', false))

		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_neuron_stability1_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%plot max epoch dur
		h=figure;
		vecMaxEpochs = vecMaxEpochDurBelowThresh/cellMultiSes{intPopulation}.samplingFreq;
		vecBins = 0:0.04:1.1;
		vecPlotBins = 0.02:0.04:1.12;
		vecV = hist(vecMaxEpochs,vecBins);
		bar(vecPlotBins,vecV,1)
		xlabel('Maximum below threshold epoch duration (s)');
		ylabel('Number of neurons');
		xlim([0 1])

		title(sprintf('%s pop%d',strSes,intPopulation))
		cellSaveMaxEpochs{intPopulation} = vecMaxEpochs;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_neuron_stability2_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%}
		%{
		%% calculate hit-correlated enhancement
		%pre-allocate
		matHCAE = nan(intNeurons,length(cellSelectContrasts));
		cellHitResp = cell(6,1);
		cellMissResp = cell(6,1);
		indHitTrials = logical(cellMultiSes{intPopulation}.structStim.vecTrialResponse);
		for intContrast=1:intContrasts
			indSelectContrast = cellSelectContrasts{intContrast};
			intCounter=intContrast;
			intHitTrialCounter = 1;
			intMissTrialCounter = 1;
			cellHitResp{intContrast} = nan(1,intNeurons*intTrials);
			cellMissResp{intContrast} = nan(1,intNeurons*intTrials);
			%select pref stim trials per neuron
			for intNeuron=1:intNeurons
				intPrefStim = vecNeuronPrefStim(intNeuron);
				indPrefStimTrials = cellSelectOri{intPrefStim};
				indPrefHitTrials = indPrefStimTrials & indHitTrials & indSelectContrast;
				indPrefMissTrials = indPrefStimTrials & ~indHitTrials & indSelectContrast;
				
				%get data
				vecHitResp = matTrialResponse(intNeuron,indPrefHitTrials);
				vecMissResp = matTrialResponse(intNeuron,indPrefMissTrials);
				cellHitResp{intContrast}(intHitTrialCounter:(intHitTrialCounter+sum(indPrefHitTrials)-1)) = vecHitResp;
				cellMissResp{intContrast}(intMissTrialCounter:(intMissTrialCounter+sum(indPrefMissTrials)-1)) = vecMissResp;
				intHitTrialCounter = intHitTrialCounter + sum(indPrefHitTrials);
				intMissTrialCounter = intMissTrialCounter + sum(indPrefMissTrials);
				
				%get normalized response to calculate hit correlated activity enhancement (range [-1 1])
				vecHitNormResp = matTrialNormResponse(intNeuron,indPrefHitTrials);
				vecMissNormResp = matTrialNormResponse(intNeuron,indPrefMissTrials);
				dblMeanHitResp = mean(vecHitNormResp);if dblMeanHitResp < 0,dblMeanHitResp = 0;end
				dblMeanMissResp = mean(vecMissNormResp);if dblMeanMissResp < 0,dblMeanMissResp = 0;end
				matHCAE(intNeuron,intContrast) = (dblMeanHitResp-dblMeanMissResp)/(dblMeanHitResp+dblMeanMissResp);
			end
			%remove nans
			cellHitResp{intContrast}(isnan(cellHitResp{intContrast})) = [];
			cellMissResp{intContrast}(isnan(cellMissResp{intContrast})) = [];
		end
		%set nans to 0
		matHCAE(isnan(matHCAE)) = 0;
		vecStimDetectActInc = mean(matHCAE(:,2:5),2);
		
		%define correlated/uncorrelated neurons
		dblFrac = 1/3;
		intNumNeurons = round(intNeurons * dblFrac);
		[vecHCAE,vecNeuronsHCAE] = findmax(vecStimDetectActInc,intNumNeurons);
		[vecHUAE,vecNeuronsHAAE] = findmin(vecStimDetectActInc,intNumNeurons);
		indNeuronsHCAE = false(size(vecStimDetectActInc));
		indNeuronsHCAE(vecNeuronsHCAE) = true;
		indNeuronsHAAE = false(size(vecStimDetectActInc));
		indNeuronsHAAE(vecNeuronsHAAE) = true;
		
		indSelectHitCorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitCorrelatedNeurons(vecNeuronsHCAE) = true;
		indSelectHitUncorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitUncorrelatedNeurons(~indNeuronsHCAE & ~indNeuronsHAAE) = true;
		
		%get which neurons are most correlated for each contrast
		matMostCorrelated = false(size(matHCAE));
		for intContrast=1:intContrasts
			indMembers = false(1,intNeurons);
			[vecHCAE_Temp,vecNeuronsHCAE_Temp] = findmax(matHCAE(:,intContrast),intNumNeurons);
			indMembers(vecNeuronsHCAE_Temp) = true;
			matMostCorrelated(:,intContrast) = indMembers;
		end
		matHCN_Consistency = corr(matMostCorrelated);
		%set auto correlation to 0
		matHCN_Consistency(diag(diag(true(size(matHCN_Consistency))))) = 0;
		
		% plot highest hit-correlated neuron consistency
		figure
		imagesc(matHCN_Consistency,[-0.4 0.4]);
		colormap(redblue(128));colorbar;drawnow;freezeColors;cbfreeze;
		vecContrast = unique(cellMultiSes{1}.structStim.Contrast)*100;
		set(gca,'XTick',1:length(vecContrast),'XTickLabel',vecContrast)
		set(gca,'YTick',1:length(vecContrast),'YTickLabel',vecContrast)
		ylabel('Contrast (%)')
		xlabel('Contrast (%)')
		title(sprintf('Correlation of highest %.0f%% of hit-correlated neurons',dblFrac*100))
		
		%save data
		cellSaveNeuronConsistency{intPopulation} = matHCN_Consistency;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_HCN_consistency_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%}
		%% calc heterogeneity
		%normalize
		matRespNormPerNeuron = zscore(matTrialResponse,[],2);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vecHeterogeneity = nan(intTrials,1);
		vecActPrefPop = nan(intTrials,1);
		vecActZPrefPop = nan(intTrials,1);
		
		matR = nan(intNeurons,intNeurons,intTrials);
		for intNeuron1=1:intNeurons
			for intNeuron2=intNeuron1:intNeurons
				matR(intNeuron1,intNeuron2,:) = matRespNormPerNeuron(intNeuron1,:).*matRespNormPerNeuron(intNeuron2,:);
			end
		end
		
		matSelectAllR = false(size(matR));
		vecR = nan(intTrials,1);
		vecR_SD = nan(intTrials,1);
		for intTrial=1:intTrials
			vecZ = matRespNormPerNeuron(:,intTrial);
			matZ1 = repmat(vecZ,[1 intNeurons]);
			matZ2 = repmat(vecZ',[intNeurons 1]);
			matDist = abs(matZ1 - matZ2);
			vecHeterogeneity(intTrial) = mean(matDist(matSelect));
			
			matSelectR = matSelectAllR;
			matSelectR(:,:,intTrial) = triu(true(intNeurons),1);
			vecR(intTrial) = mean(matR(matSelectR));
			vecR_SD(intTrial) = std(matR(matSelectR));
			
			
			%get pref pop
			%intTrialOri = [];
			vecPrefNeurons = vecNeuronPrefStim == vecStimOris(intTrial);
			vecActZPrefPop(intTrial) = mean(matRespNormPerNeuron(vecPrefNeurons,intTrial));
			vecActPrefPop(intTrial) = mean(matTrialResponse(vecPrefNeurons,intTrial));
		end
		
		%% plot sigmoids
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		vecContrasts=[0 0.5 2 8 32 100];
		vecCX=[0.2 0.5 2 8 32 100];
		if intPopulation == 1
			matMeanAct = nan(length(cellSelectContrasts),2);%[c x resp/no-resp]
			matErrAct = nan(length(cellSelectContrasts),2); %[c x resp/no-resp]
			matMeanHet = nan(length(cellSelectContrasts),2); %[c x resp/no-resp]
			matErrHet = nan(length(cellSelectContrasts),2); %[c x resp/no-resp]
		end
		for intResp=[1 2]
			indSelectResp = indSelectRespTrials == intResp-1;
			for intContrast=1:length(cellSelectContrasts)
				indSelectContrastTrials = cellSelectContrasts{intContrast};
				indSelectTrials = indSelectResp & indSelectContrastTrials;
				
				vecThisAct = vecActPrefPop(indSelectTrials);
				vecThisHet = vecHeterogeneity(indSelectTrials);
				
				matMeanAct(intContrast,intResp,intPopulation) = mean(vecThisAct);
				matErrAct(intContrast,intResp,intPopulation) = std(vecThisAct)/sqrt(length(vecThisAct));
				
				matMeanHet(intContrast,intResp,intPopulation) = mean(vecThisHet);
				matErrHet(intContrast,intResp,intPopulation) = std(vecThisHet)/sqrt(length(vecThisHet));
			end
		end
		
		h=figure;
		subplot(2,2,1)
		errorfill(vecCX,matMeanAct(:,1,intPopulation),matErrAct(:,1,intPopulation),[1 0 0],[1 0.5 0.5]);
		hold on;
		errorfill(vecCX,matMeanAct(:,2,intPopulation),matErrAct(:,2,intPopulation),[0 1 0],[0.5 1.0 0.5]);
		hold off;
		set(gca,'XTick',vecCX,'XTickLabel',vecContrasts, 'XScale','log')
		xlim([min(vecCX) max(vecCX)])
		ylabel('Pref-pop dF/F0')
		xlabel('Stimulus contrast (%)')
		
		subplot(2,2,2)
		errorfill(vecCX,matMeanHet(:,1,intPopulation),matErrHet(:,1,intPopulation),[1 0 0],[1 0.5 0.5]);
		hold on;
		errorfill(vecCX,matMeanHet(:,2,intPopulation),matErrHet(:,2,intPopulation),[0 1 0],[0.5 1.0 0.5]);
		hold off;
		set(gca,'XTick',vecCX,'XTickLabel',vecContrasts, 'XScale','log')
		xlim([min(vecCX) max(vecCX)])
		ylabel('Heterogeneity')
		xlabel('Stimulus contrast (%)')
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_resp_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end

		%save data
		cellSaveContrastResp{1} = matMeanAct;
		cellSaveContrastResp{2} = matErrAct;
		cellSaveContrastResp{3} = matMeanHet;
		cellSaveContrastResp{4} = matErrHet;
		
		%% RT dependence
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		vecRT = vecRTs(indSelectRespTrials);
		[vecRTsSorted,vecRTSortIndex] = sort(vecRT,'ascend');
		
		%heterogeneity
		vecHetHit = vecHeterogeneity(indSelectRespTrials);
		vecHetRT = vecHetHit(vecRTSortIndex);
		
		%Z-scored activity
		matActZ = matRespNormPerNeuron(:,indSelectRespTrials)';

		matZRTSorted = matActZ(vecRTSortIndex,:);
		vecZRT = mean(matZRTSorted,2);
		
		%dF/F
		matActdFoF = matTrialResponse(:,indSelectRespTrials)';

		matARTSorted = matActdFoF(vecRTSortIndex,:);
		vecART = mean(matARTSorted,2);
		
		%variance
		vecVRT = var(matARTSorted,[],2);
		
		%sparseness
		vecSparseness = kurtosis(matARTSorted,[],2)-3;
		
		%pearson-like similarity metric
		vecRHit = vecR(indSelectRespTrials);
		vecCorrLike = vecRHit(vecRTSortIndex);
		
		%spread in pearson-like similarity metric
		vecRSDHit = vecR_SD(indSelectRespTrials);
		vecCorrLikeSpread = vecRSDHit(vecRTSortIndex);
		
		%dF/F pref pop
		vecActPPHit = vecActPrefPop(indSelectRespTrials);
		vecAct_PP = vecActPPHit(vecRTSortIndex);
		
		%z-scored dF/F pref pop
		vecActZPPHit = vecActZPrefPop(indSelectRespTrials);
		vecActZ_PP = vecActZPPHit(vecRTSortIndex);
		
		
		%save
		cellSaveRTDependency{1,intPopulation} = vecRTsSorted';
		cellSaveRTDependency{2,intPopulation} = vecHetRT;
		cellSaveRTDependency{3,intPopulation} = vecZRT;
		cellSaveRTDependency{4,intPopulation} = vecART;
		cellSaveRTDependency{5,intPopulation} = vecVRT;
		cellSaveRTDependency{6,intPopulation} = vecSparseness;
		cellSaveRTDependency{7,intPopulation} = vecCorrLike;
		cellSaveRTDependency{8,intPopulation} = vecCorrLikeSpread;
		cellSaveRTDependency{9,intPopulation} = vecAct_PP;
		cellSaveRTDependency{10,intPopulation} = vecActZ_PP;
		
		%plot
		figure
		subplot(3,3,1);
		scatter(vecRTsSorted,vecHetRT,'kx')
		
		%perform regressions
		sStatsC=regstats(vecHetRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('Pop heterogeneity; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean population activity dissimilarity')
		
		%z-scored activity
		subplot(3,3,2);
		scatter(vecRTsSorted,vecZRT,'kx')
		
		%perform regressions
		sStatsC=regstats(vecZRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('Z-scored act; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean z-scored population activity')
		
		%dF/F activity
		subplot(3,3,3);
		scatter(vecRTsSorted,vecART,'kx')

		%perform regressions
		sStatsC=regstats(vecART,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('dF/F0; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean dF/F population activity')
		
		%dF/F activity
		subplot(3,3,4);
		scatter(vecRTsSorted,vecVRT,'kx')

		%perform regressions
		sStatsC=regstats(vecVRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('Variance; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Variance of population activity')
		
		
		% sparseness
		subplot(3,3,5);
		scatter(vecRTsSorted,vecSparseness,'kx')

		%perform regressions
		sStatsC=regstats(vecSparseness,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('Sparseness; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Population Sparseness')
		
		% corr-like
		subplot(3,3,6);
		scatter(vecRTsSorted,vecCorrLike,'kx')

		%perform regressions
		sStatsC=regstats(vecCorrLike,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('Pearson-like; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Pearson-like')
		
		% corr-like spread
		subplot(3,3,7);
		scatter(vecRTsSorted,vecCorrLikeSpread,'kx')

		%perform regressions
		sStatsC=regstats(vecCorrLikeSpread,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('CorrLikeSpread; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Pearson-like sd')
		
		% pref pop dF/F
		subplot(3,3,8);
		scatter(vecRTsSorted,vecAct_PP,'kx')

		%perform regressions
		sStatsC=regstats(vecAct_PP,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('PrefPop Act; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('dF/F pref-pop')

		% pref pop actZ
		subplot(3,3,9);
		scatter(vecRTsSorted,vecActZ_PP,'kx')

		%perform regressions
		sStatsC=regstats(vecActZ_PP,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		hold off
		xlim([0 3])
		title(sprintf('PrefPop ActZ; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Z-scored dF/F pref-pop')

		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_RT_dependence_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%}
		%% for same het, increase in het for fast/slow
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Contrast,unique(cellMultiSes{intPopulation}.structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn<3;
		vecStimRT = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn;
	
		%get heterogeneity
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		matPreTrialResponse = getPreTrialResponseData(cellMultiSes{intPopulation},3);

		[vecHetStim,vecActStim] = calcMatRespHeteroGen(matTrialResponse);
		[vecHetPreStim,vecActPreStim] = calcMatRespHeteroGen(matPreTrialResponse);
		
		%calculate pre-stim quintile edges
		intSegments = 3;
		[vecMeans,vecBinEdges,cellValues] = getQuintiles(vecHetPreStim,[],intSegments);
		vecBinEdges = vecBinEdges';
		
		%pre-allocate & calculate
		cellTrialIndices = cell(6,3);
		cellHetStim = cell(6,3);
		cellHetPreStim = cell(6,3);
		cellActStim = cell(6,3);
		cellActPreStim = cell(6,3);
		cellSaveChange{intPopulation} = ones(6,3,intSegments,2); %[contrast x behav x quintile x pre/post]
		
		%%
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		
		h=figure;
		for intC=1:6
			%get trials
			indC = vecStimContrasts==intC;
			dblCutOffRT = nanmedian(vecStimRT(indC));
			cellTrialIndices{intC,1} = ~indStimResp & indC;%miss
			cellTrialIndices{intC,2} = vecStimRT<dblCutOffRT & indC;%slow
			cellTrialIndices{intC,3} = vecStimRT>=dblCutOffRT & indC;%fast
			
			%get stim het
			cellHetStim{intC,1} = vecHetStim(~indStimResp & indC);%miss
			cellHetStim{intC,2} = vecHetStim(vecStimRT<dblCutOffRT & indC);%slow
			cellHetStim{intC,3} = vecHetStim(vecStimRT>=dblCutOffRT & indC);%fast
			
			%get pre-stim het
			cellHetPreStim{intC,1} = vecHetPreStim(~indStimResp & indC);%miss
			cellHetPreStim{intC,2} = vecHetPreStim(vecStimRT<dblCutOffRT & indC);%slow
			cellHetPreStim{intC,3} = vecHetPreStim(vecStimRT>=dblCutOffRT & indC);%fast
			
			%get stim act
			cellActStim{intC,1} = vecActStim(~indStimResp & indC);%miss
			cellActStim{intC,2} = vecActStim(vecStimRT<dblCutOffRT & indC);%slow
			cellActStim{intC,3} = vecActStim(vecStimRT>=dblCutOffRT & indC);%fast
			
			%get pre-stim act
			cellActPreStim{intC,1} = vecActPreStim(~indStimResp & indC);%miss
			cellActPreStim{intC,2} = vecActPreStim(vecStimRT<dblCutOffRT & indC);%slow
			cellActPreStim{intC,3} = vecActPreStim(vecStimRT>=dblCutOffRT & indC);%fast
			
			
			%plot increase
			subplot(2,3,intC)
			hold on;
			
			%miss
			[nVec,meanVecPre] = makeBins(cellHetPreStim{intC,1},cellHetPreStim{intC,1},[-inf vecBinEdges inf]);
			[nVecM,meanVecStim] = makeBins(cellHetPreStim{intC,1},cellHetStim{intC,1},[-inf vecBinEdges inf]);
			plot(repmat([1 2],[intSegments 1])',[meanVecPre' meanVecStim']','rx-')
			cellSaveChange{intPopulation}(intC,1,:,1) = meanVecPre;
			cellSaveChange{intPopulation}(intC,1,:,2) = meanVecStim;
			
			%slow
			[nVec,meanVecPre] = makeBins(cellHetPreStim{intC,2},cellHetPreStim{intC,2},[-inf vecBinEdges inf]);
			[nVecS,meanVecStim] = makeBins(cellHetPreStim{intC,2},cellHetStim{intC,2},[-inf vecBinEdges inf]);
			plot(repmat([3 4],[intSegments 1])',[meanVecPre' meanVecStim']','bx-')
			cellSaveChange{intPopulation}(intC,2,:,1) = meanVecPre;
			cellSaveChange{intPopulation}(intC,2,:,2) = meanVecStim;
			
			%fast
			[nVec,meanVecPre] = makeBins(cellHetPreStim{intC,3},cellHetPreStim{intC,3},[-inf vecBinEdges inf]);
			[nVecF,meanVecStim] = makeBins(cellHetPreStim{intC,3},cellHetStim{intC,3},[-inf vecBinEdges inf]);
			plot(repmat([5 6],[intSegments 1])',[meanVecPre' meanVecStim']','gx-')
			cellSaveChange{intPopulation}(intC,3,:,1) = meanVecPre;
			cellSaveChange{intPopulation}(intC,3,:,2) = meanVecStim;
			
			hold off;
			ylabel('Heterogeneity');
			ylim([0.5 1.5]);
			title(sprintf('Contrast %.1f%%; n-trials, M,T1=%d,T2=%d,T3=%d; \nS,T1=%d,T2=%d,T3=%d; F,T1=%d,T2=%d,T3=%d',vecContrasts(intC),nVecM,nVecS,nVecF));
		end
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sprestim_het_by_quintile__pop%draw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%%{
		%%
		%calc values from cell array
		cellCalc = cellActStim;
		matMeanHet = cellfun(@mean,cellCalc);
		matSDHet = cellfun(@std,cellCalc);
		matN = cellfun(@numel,cellCalc);
		
		%plot
		h=figure;
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrastPlot = vecContrasts;
		vecContrastPlot(1) = 0.2;
		hold on;
		for intResp=[1 3 2]
			vecColorLine = [0 0 0];
			vecColorLine(intResp) = 1;
			vecColorFill = [0.5 0.5 0.5];
			vecColorFill(intResp) = 1;
			errorfill(vecContrastPlot,matMeanHet(:,intResp)',matSDHet(:,intResp)'./sqrt(matN(:,intResp)'),vecColorLine,vecColorFill);
		end
		hold off;
		set(gca,'XScale','log')
		%}
		%%{
		%% running
		%pref pop act miss still/move + hit still/move
		%get other data
		vecOriTrials = (cellMultiSes{intPopulation}.structStim.Orientation/45)+1;
		indSelectResp = cellMultiSes{intPopulation}.structStim.vecTrialResponse;
		
		indMoved = false(1,intTrials);
		for intTrial=1:intTrials
			vecMoveSecs = cellMultiSes{intPopulation}.structStim.cellMoveSecs{intTrial};
			dblSecsOn = cellMultiSes{intPopulation}.structStim.SecsOn(intTrial);
			dblSecsOff = cellMultiSes{intPopulation}.structStim.SecsOff(intTrial);
			
			indMoved(intTrial) = any(vecMoveSecs > dblSecsOn & vecMoveSecs < dblSecsOff);
		end
		
		%figure
		hActCon = figure;
		set(hActCon,'Color',[1 1 1]);
		figure(hActCon);
		
		
		%pre-allocate
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrasts(1) = 0.2;
		intContrasts = length(vecContrasts);
		vecP = nan(1,intContrasts);
		vecMetaMoveHitY = nan(1,intContrasts);
		vecMetaMoveHitE = nan(1,intContrasts);
		vecMetaMoveMissY = nan(1,intContrasts);
		vecMetaMoveMissE = nan(1,intContrasts);
		vecMetaStillHitY = nan(1,intContrasts);
		vecMetaStillHitE = nan(1,intContrasts);
		vecMetaStillMissY = nan(1,intContrasts);
		vecMetaStillMissE = nan(1,intContrasts);
		vecMetaMoveStillY = nan(1,intContrasts);
		vecMetaMoveStillE = nan(1,intContrasts);
		intContrastCounter = 0;
		for intContrastIndex=1:intContrasts
			%get contrast
			dblContrast = vecContrasts(intContrastIndex);
			intContrastCounter = intContrastCounter + 1;
			
			%act hit still & het hit still
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & indSelectResp & ~indMoved);
			vecActHitStill = nan(1,length(vecTrials));
			vecHetHitStill = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActHitStill(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
				vecHetHitStill(intTrialCounter) = vecHeterogeneity(intTrial);
			end
			
			%act miss still & het miss still
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & ~indSelectResp & ~indMoved);
			vecActMissStill = nan(1,length(vecTrials));
			vecHetMissStill = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActMissStill(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
				vecHetMissStill(intTrialCounter) = vecHeterogeneity(intTrial);
			end
			
			%act hit move still & het hit move
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & indSelectResp & indMoved);
			vecActHitMove = nan(1,length(vecTrials));
			vecHetHitMove = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActHitMove(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
				vecHetHitMove(intTrialCounter) = vecHeterogeneity(intTrial);
			end
			
			%act miss move & het miss move
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & ~indSelectResp & indMoved);
			vecActMissMove = nan(1,length(vecTrials));
			vecHetMissMove = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActMissMove(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
				vecHetMissMove(intTrialCounter) = vecHeterogeneity(intTrial);
			end
			
			%perform analyses per contrast
			vecY = [mean(vecActHitStill) mean(vecActMissStill) mean(vecActHitMove) mean(vecActMissMove) (mean(vecActHitMove)-mean(vecActMissMove))-(mean(vecActHitStill)-mean(vecActMissStill))];
			
			%overall
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & indSelectResp);
			vecActHit = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActHit(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
			end
			
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & ~indSelectResp);
			vecActMiss = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActMiss(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
			end
			
			%heterogen
			vecHetHit = vecHeterogeneity(cellSelectContrasts{intContrastCounter} & indSelectResp)';
			vecHetMiss = vecHeterogeneity(cellSelectContrasts{intContrastCounter} & ~indSelectResp)';
			
			%put in output
			cellSaveMatContAct{intPopulation}(intContrastCounter,1) = vecY(1); %hit still
			cellSaveMatContAct{intPopulation}(intContrastCounter,2) = vecY(2); %miss still
			cellSaveMatContAct{intPopulation}(intContrastCounter,3) = vecY(3); %hit move
			cellSaveMatContAct{intPopulation}(intContrastCounter,4) = vecY(4); %miss move
			
			cellSaveMatContAct{intPopulation}(intContrastCounter,5) = vecY(5); %dMove - dStill
			cellSaveMatContAct{intPopulation}(intContrastCounter,6) = sum(indMoved)/length(indMoved); %dMove - dStill
			
			cellSaveMatContAct{intPopulation}(intContrastCounter,7) = mean(vecActHit); %overall hit
			cellSaveMatContAct{intPopulation}(intContrastCounter,8) = mean(vecActMiss); %overall miss
			
			cellSaveMatContAct{intPopulation}(intContrastCounter,9) = mean(vecHetHitStill); %hit still
			cellSaveMatContAct{intPopulation}(intContrastCounter,10) = mean(vecHetMissStill); %miss still
			cellSaveMatContAct{intPopulation}(intContrastCounter,11) = mean(vecHetHitMove); %hit move
			cellSaveMatContAct{intPopulation}(intContrastCounter,12) = mean(vecHetMissMove); %miss move
			
			cellSaveMatContAct{intPopulation}(intContrastCounter,13) = mean(vecHetHit); %overall het hit
			cellSaveMatContAct{intPopulation}(intContrastCounter,14) = mean(vecHetMiss); %overall het miss
			
			cellSaveMatContAct{intPopulation}(intContrastCounter,15) = getCohensD(vecActHitStill,vecActMissStill); %overall MES dFoF Still
			cellSaveMatContAct{intPopulation}(intContrastCounter,16) = getCohensD(vecHetHitStill,vecHetMissStill); %overall MES Het Still
			cellSaveMatContAct{intPopulation}(intContrastCounter,17) = getCohensD(vecActHitMove,vecActMissMove); %overall MES dFoF Move
			cellSaveMatContAct{intPopulation}(intContrastCounter,18) = getCohensD(vecHetHitMove,vecHetMissMove); %overall MES Het Move
			
			%put in meta vector
			vecMetaStillHitY(intContrastIndex) = vecY(1);
			vecMetaStillHitE(intContrastIndex) = std(vecActHitStill)/sqrt(length(vecActHitStill));
			vecMetaStillMissY(intContrastIndex) = vecY(2);
			vecMetaStillMissE(intContrastIndex) = std(vecActMissStill)/sqrt(length(vecActMissStill));
			vecMetaMoveHitY(intContrastIndex) = vecY(3);
			vecMetaMoveHitE(intContrastIndex) = std(vecActHitMove)/sqrt(length(vecActHitMove));
			vecMetaMoveMissY(intContrastIndex) = vecY(4);
			vecMetaMoveMissE(intContrastIndex) = std(vecActMissMove)/sqrt(length(vecActMissMove));
			vecMetaMoveStillY(intContrastIndex) = vecY(5);
			vecMetaMoveStillE(intContrastIndex) = 0;
		end
		
		%pre-compute variables
		vecWindow = 1:length(vecContrasts);
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		hold on;
		%get data
		for intResp=1:5
			if intResp == 1
				vecMeanTrace = vecMetaStillHitY(vecWindowSelect);
				vecSE = vecMetaStillHitE(vecWindowSelect);
				vecColorFill = [0 0.5 0];
				vecColorLine = [0 1 0];
			elseif intResp == 2
				vecMeanTrace = vecMetaStillMissY(vecWindowSelect);
				vecSE = vecMetaStillMissE(vecWindowSelect);
				vecColorLine = [1 0 0];
				vecColorFill = [0.5 0 0];
			elseif intResp == 3
				vecMeanTrace = vecMetaMoveHitY(vecWindowSelect);
				vecSE = vecMetaMoveHitE(vecWindowSelect);
				vecColorLine = [0 1 1];
				vecColorFill = [0 1 0.5];
			elseif intResp == 4
				vecMeanTrace = vecMetaMoveMissY(vecWindowSelect);
				vecSE = vecMetaMoveMissE(vecWindowSelect);
				vecColorLine = [1 0 1];
				vecColorFill = [1 0 0.5];
			elseif intResp == 5
				vecMeanTrace = vecMetaMoveStillY(vecWindowSelect);
				vecSE = vecMetaMoveStillE(vecWindowSelect);
				vecColorLine = [0 0 0];
				vecColorFill = [0.5 0.5 0.5];
			end
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			%fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			scatter(vecLineX,vecMeanTrace,30,vecColorLine);
		end
		hold off;
		set(gca,'XScale','log','YScale','linear')
		set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
		title(sprintf('dMove-dMiss [%d]',intPopulation))
		grid on
		xlabel('Contrast')
		ylabel('Population response dF/F0')
		xlim([min(vecContrasts(vecWindow)) max(vecContrasts(vecWindow))]);
		%ylim([-0.01 0.06])
		
		drawnow;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_dfof_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% remove different quintiles
		vecHetAll = nan(1,intTrials);
		matHetQ = nan(5,intTrials);
		for intTrial=1:intTrials
			%% STANDARD
			%perform calculation for all neurons
			vecActivity = sort(matRespNormPerNeuron(:,intTrial),'ascend');
			matZ1 = repmat(vecActivity,[1 intNeurons]);
			matZ2 = repmat(vecActivity',[intNeurons 1]);
			matDistAll = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectAll = tril(true(size(matDistAll)),-1);
			vecHetAll(intTrial) = mean(matDistAll(matSelectAll));
			
			%% REMOVE quintiles
			%do calculation
			dblSize = 0.2;
			intRemNeurons = round(intNeurons*dblSize);
			vecNeurons = 1:intNeurons;
			for intQuintile = 1:5
				intStartNeuron = floor(intNeurons*((intQuintile-1)*dblSize))+1;
				vecNeuronsQ = vecNeurons;
				vecNeuronsQ(intStartNeuron:(intStartNeuron+intRemNeurons-1)) = [];
				
				vecActivityQ = vecActivity(vecNeuronsQ);
				matZ1 = repmat(vecActivityQ,[1 length(vecActivityQ)]);
				matZ2 = repmat(vecActivityQ',[length(vecActivityQ) 1]);
				matDist = abs(matZ1 - matZ2);
				
				%save data as vectors in matrix
				matSelect = tril(true(size(matDist)),-1);
				matHetQ(intQuintile,intTrial) = mean(matDist(matSelect));
			end
		end
		cellSaveHetQs{intPopulation} = matHetQ; 
		cellSaveHetAll{intPopulation} = vecHetAll;
		cellSaveRespTimeTrials{intPopulation} = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs-cellMultiSes{intPopulation}.structStim.SecsOn; 
		cellSaveContrastTrials{intPopulation} = cellMultiSes{intPopulation}.structStim.Contrast; 
		
		%% stimulus presence decoding
		%decode 0% and 100% stimuli for presence
		%afterwards, split for hit/miss trials
		
		%define variables
		indNoStim=cellSelectContrasts{1};
		indStim=cellSelectContrasts{end};
		intNumTrials = length(indNoStim);
		indTrials = true(1,intNumTrials);
		vecStimType = cellMultiSes{intPopulation}.structStim.Contrast;
		matDetect = nan(3,6,2); %[upper/mean/lower] x [contrasts] x [hit/miss]
		vecStimDecoded = zeros(1,intNumTrials);
		matPostProb = zeros(2,intNumTrials);
		vecRatioProb = zeros(1,intNumTrials);
		vecOriTrials = (cellMultiSes{intPopulation}.structStim.Orientation/45)+1;
		for intTrial = 1:intNumTrials
			%make likelihood indices for trial
			indLikelihoodStim = indStim;
			indLikelihoodNoStim = indNoStim;
			indLikelihoodStim(intTrial) = false;
			indLikelihoodNoStim(intTrial) = false;
			
			%loop through neurons
			intNeuronCounter = 0;
			vecPrefPop = 1:length(vecNeuronPrefStim);%find(vecNeuronPrefStim == vecOriTrials(intTrial));
			vecStimP = nan(1,length(vecPrefPop));
			vecNoStimP = nan(1,length(vecPrefPop));
			for intNeuron=vecPrefPop %take only pref pop
				intNeuronCounter = intNeuronCounter + 1;
				%build likelihood
				vecLikeStimAct = matTrialResponse(intNeuron,indLikelihoodStim);
				vecLikeNoStimAct = matTrialResponse(intNeuron,indLikelihoodNoStim);
				
				%get data
				dblAct = matTrialResponse(intNeuron,intTrial);
				
				%do decoding
				dblPostStimTemp = normpdf(dblAct,mean(vecLikeStimAct),std(vecLikeStimAct));
				dblPostNoStimTemp = normpdf(dblAct,mean(vecLikeNoStimAct),std(vecLikeNoStimAct));
				
				%put in temp output
				vecStimP(intNeuronCounter) = dblPostStimTemp;
				vecNoStimP(intNeuronCounter) = dblPostNoStimTemp;
			end
			
			%get population decoding
			dblPostNoStim = prod(vecNoStimP);
			dblPostStim = prod(vecStimP);
			
			%put in output
			matPostProb(:,intTrial) = [dblPostNoStim dblPostStim];
			[dummy,vecStimDecoded(intTrial)] = max([dblPostNoStim dblPostStim]);
			vecRatioProb(:,intTrial) = dblPostStim/(dblPostNoStim+dblPostStim);
		end
		
		%decoded stim presence probability + 95% CI
		vecHeterogeneity = calcMatRespHeteroGen(matTrialResponse);
		intHighLow = 2;
		intActTypes = 2;
		intResp_NoResp = 2;
		matDetectLowHighActHetRNR = nan(intContrasts,intHighLow,intActTypes,intResp_NoResp);
		for intContrast=1:length(vecContrasts)
			vecResp = vecStimDecoded(cellSelectContrasts{intContrast} & indSelectResp) == 2;
			vecNoResp = vecStimDecoded(cellSelectContrasts{intContrast} & ~indSelectResp) == 2;
			
			[dblP,dblCI] = binofit(sum(vecResp),length(vecResp));
			matDetect(1,intContrast,1) = dblCI(2); %upper
			matDetect(2,intContrast,1) = dblP; %mean
			matDetect(3,intContrast,1) = dblCI(1); %lower
			
			
			[dblP,dblCI] = binofit(sum(vecNoResp),length(vecNoResp));
			matDetect(1,intContrast,2) = dblCI(2); %upper
			matDetect(2,intContrast,2) = dblP; %mean
			matDetect(3,intContrast,2) = dblCI(1); %lower
			
			%get act+het for these trials
			vecThisDecoded = vecStimDecoded(cellSelectContrasts{intContrast});
			vecAct = mean(matTrialResponse(:,cellSelectContrasts{intContrast}),1);
			vecHet = vecHeterogeneity(cellSelectContrasts{intContrast})';
			
			%act
			[dummy,vecSorted]=sort(vecAct);
			vecRespTrials = find(indSelectResp(cellSelectContrasts{intContrast}));
			vecNoRespTrials = find(~indSelectResp(cellSelectContrasts{intContrast}));
			vecLowAct = vecSorted(1:floor(length(vecSorted)/2));
			vecHighAct = vecSorted(ceil((length(vecSorted)/2)+0.1):end);
			vecRespLowAct = vecLowAct(ismember(vecLowAct,vecRespTrials));
			vecRespHighAct = vecHighAct(ismember(vecHighAct,vecRespTrials));
			vecNoRespLowAct = vecLowAct(ismember(vecLowAct,vecNoRespTrials));
			vecNoRespHighAct = vecHighAct(ismember(vecHighAct,vecNoRespTrials));
			
			dblP_LowActResp = binofit(sum(vecThisDecoded(vecRespLowAct) == 2),length(vecRespLowAct));
			dblP_HighActResp = binofit(sum(vecThisDecoded(vecRespHighAct) == 2),length(vecRespHighAct));
			dblP_LowActNoResp = binofit(sum(vecThisDecoded(vecNoRespLowAct) == 2),length(vecNoRespLowAct));
			dblP_HighActNoResp = binofit(sum(vecThisDecoded(vecNoRespHighAct) == 2),length(vecNoRespHighAct));
			
			
			matDetectLowHighActHetRNR(intContrast,1,1,1) = dblP_LowActResp; %low dF/F resp
			matDetectLowHighActHetRNR(intContrast,2,1,1) = dblP_HighActResp; %high dF/F resp
			matDetectLowHighActHetRNR(intContrast,1,1,2) = dblP_LowActNoResp; %low dF/F no-resp
			matDetectLowHighActHetRNR(intContrast,2,1,2) = dblP_HighActNoResp; %high dF/F no-resp
			
			%het
			[dummy,vecSorted]=sort(vecHet);
			vecLowHet = vecSorted(1:floor(length(vecSorted)/2));
			vecHighHet = vecSorted(ceil((length(vecSorted)/2)+0.1):end);
			vecRespLowHet = vecLowHet(ismember(vecLowHet,vecRespTrials));
			vecRespHighHet = vecHighHet(ismember(vecHighHet,vecRespTrials));
			vecNoRespLowHet = vecLowHet(ismember(vecLowHet,vecNoRespTrials));
			vecNoRespHighHet = vecHighHet(ismember(vecHighHet,vecNoRespTrials));
			
			dblP_LowHetResp = binofit(sum(vecThisDecoded(vecRespLowHet) == 2),length(vecRespLowHet));
			dblP_HighHetResp = binofit(sum(vecThisDecoded(vecRespHighHet) == 2),length(vecRespHighHet));
			dblP_LowHetNoResp = binofit(sum(vecThisDecoded(vecNoRespLowHet) == 2),length(vecNoRespLowHet));
			dblP_HighHetNoResp = binofit(sum(vecThisDecoded(vecNoRespHighHet) == 2),length(vecNoRespHighHet));
			
			matDetectLowHighActHetRNR(intContrast,1,2,1) = dblP_LowHetResp; %low het resp
			matDetectLowHighActHetRNR(intContrast,2,2,1) = dblP_HighHetResp; %high het resp
			matDetectLowHighActHetRNR(intContrast,1,2,2) = dblP_LowHetNoResp; %low het no-resp
			matDetectLowHighActHetRNR(intContrast,2,2,2) = dblP_HighHetNoResp; %high het no-resp
		end
		
		
		
		%plot
		hDecResp = figure;
		hold on;
		
		vecWindow = 1:length(vecContrasts);
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		
		for intDetect=[0 1]
			if intDetect == 1
				vecMeanTrace = matDetect(2,:,1);
				vecMinTrace = matDetect(3,:,1);
				vecMaxTrace = matDetect(1,:,1);
				vecColorFill = [0.7 1.0 0.7];
				vecColorLine = [0 1 0];
			else
				vecMeanTrace = matDetect(2,:,2);
				vecMinTrace = matDetect(3,:,2);
				vecMaxTrace = matDetect(1,:,2);
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
			end
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
		set(gca,'XScale','log','YScale','linear')
		set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
		title('Stimulus presence decoding with 95% CI')
		grid on
		xlabel('Contrast')
		ylabel('Decoded stimulus presence')
		xlim([min(vecContrasts(vecWindow)) max(vecContrasts(vecWindow))])
		ylim([0 1])
		legend({'95% CI','Miss','95% CI','Hit'},'Location','Best')
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_decodeStimPresence%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		matMeanDetect = squeeze(matDetect(2,:,:)); %[contrasts] x [hit/miss]
		cellSaveMatDetect{intPopulation} = matMeanDetect; %[contrasts] x [hit/miss]
		cellSaveMatDetectLowHighActHet{intPopulation} = matDetectLowHighActHetRNR;
		
		%% brain state analysis
		%take blocks of time of good/bad performance, then look at neural
		%correlates
		% => canned
		
		%% investigate multiplicative gain for 33% HCN from miss to hit
		% => canned
		%}
	end
	%{
	%overall figure
	matMeanAct = mean(matMeanAct,3);
	matMeanHet = mean(matMeanHet,3);
	
	h=figure;
	subplot(2,2,1)
	plot(vecCX,matMeanAct(:,1),'r');
	hold on;
	plot(vecCX,matMeanAct(:,2),'g');
	hold off;
	set(gca,'XTick',vecCX,'XTickLabel',vecContrasts, 'XScale','log')
	xlim([min(vecCX) max(vecCX)])
	ylabel('Pref-pop dF/F0')
	xlabel('Stimulus contrast (%)')
	
	subplot(2,2,2)
	plot(vecCX,matMeanHet(:,1),'r');
	hold on;
	plot(vecCX,matMeanHet(:,2),'g');
	hold off;
	set(gca,'XTick',vecCX,'XTickLabel',vecContrasts, 'XScale','log')
	xlim([min(vecCX) max(vecCX)])
	ylabel('Heterogeneity')
	xlabel('Stimulus contrast (%)')
	
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%sresp_over_contrasts_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	%}
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strRecs = num2str(vecRecordings);
		strFile = ['data_aggregate' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
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