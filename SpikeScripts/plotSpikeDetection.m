%load data
for intMouse=5
	close all
	drawnow;
	clearvars -except intMouse
	
	%use neuropil subtraction?
	boolUseNeuropilSubtraction = false;
	
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
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['D:\Data\Results\spikeAnalysis\dataPreProAssemblies' strSes '.mat']);
	end
	
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	%#ok<*AGROW>
	strSes = ['AA' strSes]; 
	
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
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
		
		%get response matrix
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		%get behavioral variables
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs -  cellMultiSes{intPopulation}.structStim.SecsOn;
		vecRTs(indMiss) = nan;
		indFast = vecRTs < nanmedian(vecRTs);
		indSlow = vecRTs >= nanmedian(vecRTs);
		dblCutOffRT = nanmedian(vecRTs);
		
		%get heterogeneity
		%[vecFrameHeterogeneity,vecFrameActivity] = calcSlidingHeteroGen(cellMultiSes{intPopulation});
		%[vecTrialNormHeterogeneity,vecTrialNormActivity] = calcMatRespHeteroGen(matRespNormPerContrast);
		
		%get assembly variables
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		%matRespAssemblies = sAssemblies{intPopulation}.matRespAssemblies; 
		matReorderedAssemblyCorrs=sAssemblies{intPopulation}.matAssemblyCorrelations(sAssemblies{intPopulation}.vecReorder,sAssemblies{intPopulation}.vecReorder);
		matM = getBlockMeans(matReorderedAssemblyCorrs,sAssemblies{intPopulation}.vecAssemblyIdentity);
		
		%get assembly activity during trials
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		dblStimSecs = 1;
		sStimTemp.FrameOff = sStimTemp.FrameOn + round(dblSamplingFreq*dblStimSecs);
		matRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp);%[assemblies x trials]
		
		%detect spikes
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
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
		
		
		%%
		intNeuron = 1
			vecFrames = 10000:13000;
			vecAct = cellMultiSes{intPopulation}.neuron(intNeuron).dFoF(vecFrames);
			%vecAct =vecAct+mean(vecAct);
			%vecFilt = normpdf(-2:1:2,0,1);
			%vecAct = conv(vecAct,vecFilt/sum(vecFilt),'same');
			%vecAct = (cellMultiSes{intPopulation}.neuron(intNeuron).dFoF(vecFrames-1) + cellMultiSes{intPopulation}.neuron(intNeuron).dFoF(vecFrames) + cellMultiSes{intPopulation}.neuron(intNeuron).dFoF(vecFrames+1))/3;
			[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(vecAct,dblSamplingFreq,0.5);
			
			clf
			plot(vecFrames/dblSamplingFreq,vecAct,'b')
			hold on
			plot(vecFrames/dblSamplingFreq,expFit,'g')
			xlim(vecFrames([1 end])/dblSamplingFreq)
		
		%% get population events
		fprintf('Retrieving population spiking data for assembly formation analysis. This might take a while...   [%s]\n',getTime);
		
		%get number of active neurons per frame
		vecNrActive=sum(matSpikeCountsHD>0,1);
		
		%select which frames are potentially interesting
		vecPotentialAssemblyFrames = vecNrActive > 0;
		
		%remove frames where single neurons fire single spikes in temporal isolation
		vecSingleSpikes = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))<2 & vecNrActive(2:(end-1))>0 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecSingleSpikes false]) = 0;
		
		%get frames where multiple neurons are firing in temporal isolation
		vecTransientAssemblies = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))>2 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecTransientAssemblies false]) = 1;
		
		%define start frames of assemblies
		vecAssemblyStarts = [false vecPotentialAssemblyFrames(1:(end-1))==0 & vecPotentialAssemblyFrames(2:(end))>0];
		vecAssemblyStops = [false vecPotentialAssemblyFrames(1:(end-1))>0 & vecPotentialAssemblyFrames(2:end)==0];
		
		%combine
		vecAssemblyStarts = find(vecAssemblyStarts);
		vecAssemblyStops = find(vecAssemblyStops);
		if vecAssemblyStarts(1) > vecAssemblyStops(1)
			vecAssemblyStarts = [1 vecAssemblyStarts];
		end
		[dummy,vecReorder] = sort(vecAssemblyStarts,'ascend');
		if length(vecAssemblyStops) < length(vecAssemblyStarts)
			vecAssemblyStops(end+1) = length(vecPotentialAssemblyFrames);
		end
		vecAssemblyStarts = vecAssemblyStarts(vecReorder);
		vecAssemblyStops = vecAssemblyStops(vecReorder);
		
		%get assemblies from spike matrix + assembly stop/starts
		matAssemblies = getAssemblies(matSpikeCountsHD,vecAssemblyStarts,vecAssemblyStops)>0;
		
		%remove assemblies with 2 or fewer neurons
		indRealAssemblies = sum(matAssemblies,1)>2 & sum(matAssemblies,1) < intNeurons;
		matAssemblies = matAssemblies(:,indRealAssemblies);
		vecAssemblyStarts = vecAssemblyStarts(indRealAssemblies);
		vecAssemblyStops = vecAssemblyStops(indRealAssemblies);
		
		%get subset
		vecFrames2 = 5100:5350;
		indAssemblySub = vecAssemblyStarts < vecFrames2(end) & vecAssemblyStops > vecFrames2(1);
		
		vecStarts = vecAssemblyStarts(indAssemblySub) - vecFrames2(1) + 1;
		vecStops = vecAssemblyStops(indAssemblySub) - vecFrames2(1) + 1;
		
		matAEs = matSpikeCountsHD(:,vecFrames2)>0;
		matPlotted = false(size(matAEs));
		
		%% plot
		figure
		hold on
		matColormap = jet(length(vecStarts));
		line(([0;300]+vecFrames2(1))/dblSamplingFreq,[intNeurons intNeurons]+5,'Color',[0.5 0.5 0.5],'LineWidth',5);
		for intPopulationEvent=1:length(vecStarts)
			intStart = vecStarts(intPopulationEvent);
			intStop = vecStops(intPopulationEvent);
			matSubSet = matAEs(:,intStart:intStop);
			matPlotted(:,intStart:intStop) = matSubSet;
			[vecNeurons,vecTimePoints] = find(matSubSet);
			line(([vecTimePoints';vecTimePoints']+intStart-1+vecFrames2(1))/dblSamplingFreq,[vecNeurons'-0.45;vecNeurons'+0.45],'Color',matColormap(intPopulationEvent,:),'LineWidth',3);
			
			line(([vecTimePoints(1)-0.5;vecTimePoints(end)+0.5]+intStart-1+vecFrames2(1))/dblSamplingFreq,[intNeurons intNeurons]+5,'Color',matColormap(intPopulationEvent,:),'LineWidth',5);
		end
		matSubSet = matAEs & ~matPlotted;
		[vecNeurons,vecTimePoints] = find(matSubSet);	
		line(([vecTimePoints';vecTimePoints']+vecFrames2(1))/dblSamplingFreq,[vecNeurons'-0.45;vecNeurons'+0.45],'Color',[0.5 0.5 0.5],'LineWidth',3);
		xlim([202 210])
		ylabel('Neuron')
		xlabel('Time (s)')
		hold off
		
		%% create miss/slow/fast selection vectors
		intFrames = size(matSpikeCountsHD,2);
		%miss
		vecMissOn = cellMultiSes{intPopulation}.structStim.FrameOn(indMiss);
		vecMissOff = vecMissOn + round(dblSamplingFreq*dblStimSecs);
		indMissHD = false(1,intFrames);
		for intTrial=1:length(vecMissOn)
			indMissHD(vecMissOn(intTrial):vecMissOff(intTrial)) = true;
		end
		vecMissHD = find(indMissHD);
		%slow
		vecSlowOn = cellMultiSes{intPopulation}.structStim.FrameOn(indSlow);
		vecSlowOff = vecSlowOn + round(dblSamplingFreq*dblStimSecs);
		indSlowHD = false(1,intFrames);
		for intTrial=1:length(vecSlowOn)
			indSlowHD(vecSlowOn(intTrial):vecSlowOff(intTrial)) = true;
		end
		vecSlowHD = find(indSlowHD);
		%fast
		vecFastOn = cellMultiSes{intPopulation}.structStim.FrameOn(indFast);
		vecFastOff = vecFastOn + round(dblSamplingFreq*dblStimSecs);
		indFastHD = false(1,intFrames);
		for intTrial=1:length(vecFastOn)
			indFastHD(vecFastOn(intTrial):vecFastOff(intTrial)) = true;
		end
		vecFastHD = find(indFastHD);
		
		%% plot temporal sequence
		%get subset
		vecOccFast = find(ismember(vecAssemblyStarts,vecFastHD));
		intOcc = vecOccFast(10);
		vecFrames3 = vecAssemblyStarts(intOcc):(vecAssemblyStops(intOcc)-1);
		matAEs = matSpikeCountsHD(:,vecFrames3);
		
		matPlotted = false(size(matAEs));
		matAEs = double(matAEs>0);
		dblCoM = centerOfMass(matAEs,2);
		vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblCoM)/dblSamplingFreq;
		
		matAEs = double(matAEs>0);
		matAEs(matAEs==0) = nan;
		matPlot = bsxfun(@times,matAEs,vecThisMuCC);
		
		vecRange = max(abs(matPlot(:)))*[-1 1];
		figure
		implot(matPlot,vecRange,colormap(redbluepurple(128)),3);caxis(vecRange*1000);colorbar;
		hold on
		plot([dblCoM dblCoM],get(gca,'ylim'),'k--')
		hold off
		set(gca,'xtick',1:length(vecFrames3),'xticklabel',vecFrames3/dblSamplingFreq)
		xlabel('Time (s)');
		ylabel('Neuron');
		title(sprintf('Activation sequence in population event (%d)',intOcc));
		
		
		%% 
		figure
		hold on
		vecAssemblies = find(indAssemblySub);
		matColormap = jet(length(vecAssemblies));
		intCounter = 0;
		for intAssembly=1:size(matAssemblies,2)
			vecNeurons = find(matAssemblies(:,intAssembly));
			if ismember(intAssembly,vecAssemblies)
				intCounter = intCounter + 1;
				vecColor = matColormap(intCounter,:);
				intW = 3;
			else
				vecColor = [0.5 0.5 0.5];
				intW = 2;
			end
			line([repmat(intAssembly,[1 length(vecNeurons)]);repmat(intAssembly,[1 length(vecNeurons)])],[vecNeurons'-0.45;vecNeurons'+0.45],'Color',vecColor,'LineWidth',intW);
		end
		ylim([0 intNeurons+0.5]);
		xlim([200 300]);
		%xlim([0 size(matAssemblies,2)+1]);
		ylabel('Neuron')
		xlabel('Population Event #')
		hold off
		
		
		%% plot location members assemblies
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
		
		% get assembly location data
		matLocationsXY = nan(intNeurons,2);
		for intNeuron = 1:intNeurons
			matLocationsXY(intNeuron,1) = cellMultiSes{intPopulation}.neuron(intNeuron).x;
			matLocationsXY(intNeuron,2) = cellMultiSes{intPopulation}.neuron(intNeuron).y;
		end
		dblPix2Micron = cellMultiSes{intPopulation}.xml.dblActualVoxelSizeX/1000;
		matLocationsXY_microns = matLocationsXY*dblPix2Micron; %transform pixels to microns
		matDiffX = abs(bsxfun(@minus,matLocationsXY_microns(:,1),matLocationsXY_microns(:,1)'));
		matDiffY = abs(bsxfun(@minus,matLocationsXY_microns(:,2),matLocationsXY_microns(:,2)'));
		matDist = sqrt(matDiffX.^2+matDiffY.^2);
		
		%compare mean distance of core assembly members to random group of same size
		vecPotentialMembers = find(any(matAssemblyCoreMembers,1));
		vecMembersAssembly2 = find(matAssemblyCoreMembers(2,:));
		
		%get real dist
		matSelect = tril(true(length(vecMembersAssembly2)),-1);
		matRealDist = matDist(vecMembersAssembly2,vecMembersAssembly2);
		dblMeanDist = mean(matRealDist(matSelect));

		figure
		scatter(matLocationsXY(vecPotentialMembers,1),matLocationsXY(vecPotentialMembers,2),'bo');
		hold on
		scatter(matLocationsXY(vecMembersAssembly2,1),matLocationsXY(vecMembersAssembly2,2),'rx');
		plot([1 1],[1 512],'k');
		plot([1 512],[1 1],'k');
		hold off
		xlim([0 512]);
		ylim([0 512]);
		title(sprintf('Mean inter-member distance: %.3f micron',dblMeanDist))
		return
	end
end