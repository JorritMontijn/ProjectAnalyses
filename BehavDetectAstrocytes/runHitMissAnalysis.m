clear all;
%close all;
for intMouse=5
	clearvars -except intMouse
	intUseNeuropilSubtraction = -1; %[-1=no prepro,0=none,1=pre,2=post]
	boolExcludeLocomotor = false;
	strAnalyzeType = 'neuron'; %'neuron' or 'astrocyte'
	boolUseDoPEP = false; %doPEP plots population mean across trials; if false, it plots mean+SEM across neurons
	
	%% load data
	loadSes;
	
	%% run populations
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		%msg
		fprintf('Now processing %s [pop %d], analyzing %ss [%s]\n',strSes,intPopulation,strAnalyzeType,getTime);
		
		%recalc dfof
		if intUseNeuropilSubtraction == 2
			dblNeuropilSubtractionFactor = [];
			cellMultiSes{intPopulation} = doRecalcdFoF(cellMultiSes{intPopulation},5,[],strAnalyzeType,[],[],dblNeuropilSubtractionFactor);
		elseif intUseNeuropilSubtraction == 0
			cellMultiSes{intPopulation} = doRecalcdFoF(cellMultiSes{intPopulation},3,[],strAnalyzeType);
		else
			%do nothing
		end
		
		%replace neurons with astrocytes
		sObject = cellMultiSes{intPopulation}.(strAnalyzeType);
		boolOnlyPresence = true;
		
		%% run header
		runHitMissHeader;
		intTotFrameNr = numel(sObject(1).dFoF);
		dblSampFreq = cellMultiSes{1}.structStim.FrameOff(end)/cellMultiSes{1}.structStim.SecsOff(end);
		
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
		
		%% run analysis
		intObjects = numel(sObject);
		if boolUseDoPEP
			%set PEP options
			sOptionsHit = struct;
			sOptionsHit.handleFig = 1;
			sOptionsHit.vecColor = [0 1 0];
			sOptionsHit.vecWindow = [-3 5]*dblSampFreq;
			sOptionsMiss = sOptionsHit;
			sOptionsMiss.vecColor = [1 0 0];
			figure;
			%get data
			vecDataAggregate = zeros(size(sObject(1).dFoF));
			vecTimestamps = 1:numel(vecDataAggregate);
			for intObject=1:intObjects
				vecData = sObject(intObject).dFoF-median(sObject(intObject).dFoF);
				vecDataAggregate = vecDataAggregate + vecData;
			end
			for intContrast=1:intContrasts
				
				%select trials
				indContrastTrials = cellSelectContrasts{intContrast};
				vecHits = find(indContrastTrials & indStimResp);
				vecMisses = find(indContrastTrials & ~indStimResp);
				
				%get starting frames
				vecHitStarts = structStim.FrameOn(vecHits);
				vecMissStarts = structStim.FrameOn(vecMisses);
				
				%plot
				subplot(2,3,intContrast);
				[vecHitMean,vecHitSEM] = doPEP(vecTimestamps,vecDataAggregate,vecHitStarts,sOptionsHit);
				[vecMissMean,vecMissSEM] = doPEP(vecTimestamps,vecDataAggregate,vecMissStarts,sOptionsMiss);
				ylim([-2 5])
			end
			
		else
			for intObject=1:intObjects
				vecData = sObject(intObject).dFoF-median(sObject(intObject).dFoF);
				
				
				%% do analysis
				%prep
				cellContsMiss = cell(1,intContrasts);
				cellContsHit = cell(1,intContrasts);
				vecPEP = round(-3*dblSampFreq):round(5*dblSampFreq);
				%go through contrasts
				%clf;
				for intContrast=1:intContrasts
					
					%select trials
					indContrastTrials = cellSelectContrasts{intContrast};
					vecHits = find(indContrastTrials & indStimResp);
					vecMisses = find(indContrastTrials & ~indStimResp);
					
					%prep data
					matHitPEP = nan(numel(vecHits),numel(vecPEP));
					matMissPEP = nan(numel(vecMisses),numel(vecPEP));
					
					%go through hits & misses
					for intHit=1:numel(vecHits)
						intTrial = vecHits(intHit);
						vecSelect = vecPEP+structStim.FrameOn(intTrial);
						indKeep = ~(vecSelect<1 | vecSelect>intTotFrameNr);
						matHitPEP(intHit,indKeep) = vecData(vecSelect(indKeep));
					end
					for intMiss=1:numel(vecMisses)
						intTrial = vecMisses(intMiss);
						vecSelect = vecPEP+structStim.FrameOn(intTrial);
						indKeep = ~(vecSelect<1 | vecSelect>intTotFrameNr);
						matMissPEP(intMiss,indKeep) = vecData(vecSelect(indKeep));
					end
					%save data
					cellSaveRespPEP{intPopulation}{intObject,intContrast,1} = matMissPEP;
					cellSaveRespPEP{intPopulation}{intObject,intContrast,2} = matHitPEP;
				end
			end
			%% plot
			figure;
			for intContrast=1:intContrasts
				cellMiss = cellSaveRespPEP{intPopulation}(:,intContrast,1);
				matMeanMiss = cell2mat(cellfun(@nanmean,cellMiss,cellfill(1,size(cellMiss)),'UniformOutput',false));
				cellSaveMeanMiss{intPopulation}(:,:,intContrast) = matMeanMiss;
				
				cellHit = cellSaveRespPEP{intPopulation}(:,intContrast,2);
				matMeanHit = cell2mat(cellfun(@nanmean,cellHit,cellfill(1,size(cellHit)),'UniformOutput',false));
				cellSaveMeanHit{intPopulation}(:,:,intContrast) = matMeanHit;
				
				%plot
				subplot(2,3,intContrast)
				errorfill(vecPEP/dblSampFreq,nanmean(matMeanHit,1),nanstd(matMeanHit,[],1)/sqrt(numel(vecHits)),[0 1 0],[0.7 1 0.7]);
				hold on
				errorfill(vecPEP/dblSampFreq,nanmean(matMeanMiss,1),nanstd(matMeanMiss,[],1)/sqrt(numel(vecMisses)),[1 0 0],[1 0.7 0.7]);
				hold off
				ylim([-0.015 0.035])
			end
			drawnow;
		end
	end
end
