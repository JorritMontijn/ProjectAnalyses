sParams.boolSaveData = true;
sParams.boolSavePlots = false;
for intMouse=[2 3]
	close all
	clearvars -except sParams intMouse strOldDir
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
		load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat']);
	end
	
	sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
	strOldDir = cd(sParams.strFigDir);
	

	for intPopulation=1:length(cellMultiSes)
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
		
		indStimResp = structStim.vecTrialResponse==1;
		
		
		%% calculate response to preferred stimulus over time
		%calculate response to preferred orientation over time
		matSelectC = repmat(cellSelectContrasts{end},[intNeurons 1]);
		vecLookupC = unique(cellMultiSes{intPopulation}.structStim.Contrast);
		vecLookupO = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		matPrefO = vecLookupO(repmat(vecNeuronPrefStim',[1 intTrials]));
		matStimO = repmat(cellMultiSes{intPopulation}.structStim.Orientation,[intNeurons 1]);
		matSelectO = matPrefO == matStimO;
		matSelect = matSelectC & matSelectO;
		%vecReps = sum(matSelect,2);
		%if any(vecReps ~= vecReps(1)),error([mfilename ':UnevenRepetitions'],'Number of repetitions unequal; pop%d',intPopulation);end
		%matPrefResp = reshape(matTrialResponse(matSelect),size(matSelect,1),vecReps(1));
		
		vecTrialOri = nan(1,intTrials);
		for intOri=1:length(vecLookupO)
			vecTrialOri(vecLookupO(intOri)==cellMultiSes{intPopulation}.structStim.Orientation) = intOri;
		end
		
		vecTrialContrast = nan(1,intTrials);
		for intC=1:length(vecLookupC)
			vecTrialContrast(vecLookupC(intC)==cellMultiSes{intPopulation}.structStim.Contrast) = intC;
		end
		
		%% get contrast responses
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		vecContrasts=[0 0.5 2 8 32 100];
		vecCX=[0.2 0.5 2 8 32 100];
		if intPopulation == 1
			cellHit = cell(length(cellSelectContrasts),1);%[c x resp/no-resp]
			cellMiss = cell(length(cellSelectContrasts),1); %[c x resp/no-resp]
		end
		for intTrial=1:intTrials
			intC = vecTrialContrast(intTrial);
			intO = vecTrialOri(intTrial);
			
			vecPrefPop = vecNeuronPrefStim==intO;
			dblMeanPrefPopAct = mean(matTrialResponse(vecPrefPop,intTrial));
			if indSelectRespTrials(intTrial) == 1
				cellHit{intC} = [cellHit{intC} dblMeanPrefPopAct];
			else
				cellMiss{intC} = [cellMiss{intC} dblMeanPrefPopAct];
			end
		end
		
		%% rise-time analysis
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
		intActTypes = 2; %hetero - dF/F
		intRespTypes = 3; %miss - slow - fast
		matSlopes = nan(length(vecC),intRespTypes,intActTypes);
		matRiseTimes = nan(length(vecC),intRespTypes,intActTypes);
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
					
					%get slopes/rise times
					vecBaseline = round(([-1 0]-vecWindowSecs(1))*cellMultiSes{intPopulation}.samplingFreq);
					dblBase = mean(vecMeanTrace(vecBaseline(1):vecBaseline(end)));
					dblMax = max(vecMeanTrace(vecBaseline(end):end));
					
					dblLow = (dblMax-dblBase)*0.1 + dblBase;
					dblHigh = (dblMax-dblBase)*0.9 + dblBase;
					
					indPoints = false(size(vecMeanTrace));
					indPoints(vecBaseline(end):end) = true;
					intPeak = find(vecMeanTrace > dblHigh & indPoints,1,'first');
					indPoints = false(size(vecMeanTrace));
					indPoints(1:intPeak) = true;
					intStart = find(vecMeanTrace < dblLow & indPoints,1,'last');
					
					if ~isempty(intPeak) && ~isempty(intStart)
						%rise time
						dblRiseTime = (intPeak-intStart)/cellMultiSes{intPopulation}.samplingFreq;
						
						%slope
						dblSlope = (dblHigh-dblLow)/dblRiseTime;
						
						%save data
						matSlopes(intC,intType,intAct+1) = dblSlope;
						matRiseTimes(intC,intType,intAct+1) = dblRiseTime;
					end
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
					ylim([-0.02 0.08])
				end
				grid on
				xlabel('Time after stim onset (s)')
				xlim(vecWindowSecs)
				drawnow;
			end
		end
		
		cellSaveNeuralRespSpeed{intPopulation,1} = matSlopes;
		cellSaveNeuralRespSpeed{intPopulation,2} = matRiseTimes;
		
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
		
	end
	%% plot overall sigmoid
	if sParams.boolSavePlots
		h=figure;
		errorfill(vecCX,cellfun(@mean,cellMiss)',cellfun(@std,cellMiss)' ./ sqrt(cellfun(@numel,cellMiss))',[1 0 0],[1 0.5 0.5]);
		hold on;
		errorfill(vecCX,cellfun(@mean,cellHit)',cellfun(@std,cellHit)' ./ sqrt(cellfun(@numel,cellHit))',[0 1 0],[0.5 1.0 0.5]);
		hold off;
		set(gca,'XTick',vecCX,'XTickLabel',vecContrasts, 'XScale','log')
		xlim([min(vecCX) max(vecCX)])
		ylabel('Pref-pop dF/F0')
		xlabel('Stimulus contrast (%)')
		
		%t-tests
		if any(cellfun(@isempty,cellMiss)),cellMiss{cellfun(@isempty,cellMiss)}=nan;end
		if any(cellfun(@isempty,cellHit)),cellHit{cellfun(@isempty,cellHit)}=nan;end
		[h,p]=cellfun(@ttest2,cellHit,cellMiss);
		[h,j,pCorr]=fdr_bh(p);
		
		[hOverall,pOverall]=ttest2(cell2mat(cellHit(2:end-1)'),cell2mat(cellMiss(2:end-1)'));
		
		title(sprintf('%s,T-tests; overall,p=%.3f;0,p=%.3f;0.5,p=%.3f;2,p=%.3f;8,p=%.3f;32,p=%.3f;100,p=%.3f',strSes,pOverall,pCorr))
		
		drawnow;
		strFig = sprintf('%sagg_resp_over_contrasts_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	%% save data
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['data_aggregateX' strSes '_' strrep(strDate,'    ','_')];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end
cd(strOldDir);