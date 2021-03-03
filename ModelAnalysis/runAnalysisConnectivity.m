strPath = 'F:\Data\Results\SimAggregates\';
%strSearchFile = '*LargeRetOriTuning*';
strSearchFile = '*LargeRetNoise00_2021-01-12*';
sFiles=dir(fullfile(strPath,strSearchFile));
intUseFile = ~contains({sFiles(:).name},'prepro');
strSourceFile = sFiles(intUseFile).name;
strSourcePath = sFiles(intUseFile).folder;

sLoad=load(fullfile(strSourcePath,strSourceFile));

cellSplit = strsplit(strSourceFile,'.');
strSourceNoExt = strjoin(cellSplit(1:(end-1)),'.');
cellSource = strsplit(strSourceNoExt,'_');
strFigBase = strjoin(cellSource((end-1):end),'_');

%%
dblStepSize = (sLoad.vecOverallT(end)-sLoad.vecOverallT(1))/(numel(sLoad.vecOverallT)-1);
cellSpikeTimes = doSpikeIntToDouble(sLoad.cellSpikeTimesCortex,dblStepSize,sLoad.vecOverallT(1));

vecNrSpikes = cellfun(@numel,sLoad.cellSpikeTimesCortex);
vecStimOn = sLoad.vecStimStartSecs;
vecStimOff = sLoad.vecStimStopSecs;
vecOrientation = sLoad.vecTrialOris;
vecMeanR = mean(matSpikeRate,2);
dblTrialDur = mean(vecStimOff-vecStimOn);

%%
matSpikeRate = getSpikeCounts(cellSpikeTimes,vecStimOn,vecStimOff)./(vecStimOff-vecStimOn);
[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matSpikeRate,vecOrientation);
matMeanNS = mean(matRespNSR,3);
dblThreshold = 8;%9.5
vecMaxResp = max(matMeanNS,[],2); %in any orientation
vecAQ=vecMaxResp>dblThreshold;
vecEI = sLoad.vecCellTypes;

%% prep
vecPrefOri = sLoad.vecPrefOri*2;
vecUniques = unique(roundi(diff(sort(vecPrefOri)),6));
dblStep = vecUniques(2);
vecBinEdges = (-pi-(dblStep/2)):dblStep:(pi+(dblStep/2));
vecBinCenters = vecBinEdges(2:end)-(dblStep/2);
matSynFromTo = sLoad.matSynFromTo;

%% plot
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	for intTargetAQ=[1 0]
		if intTargetAQ==1
			strTargetAQ = 'active';
		else
			strTargetAQ = 'silent';
		end
		
		
			for intSourceEI=1:2
				if intSourceEI==1
					strSourceEI = 'E';
				else
					strSourceEI = 'I';
				end
				for intSourceAQ=[1 0]
			if intSourceAQ==1
				strSourceAQ = 'active';
			else
				strSourceAQ = 'silent';
			end
				vecSourceN = find(intSourceAQ(:) == vecAQ & vecEI(:)==intSourceEI);
				vecTargetN = find(intTargetAQ(:) == vecAQ & vecEI(:)==intTargetEI);
				
				vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
				vecSynW = sLoad.vecSynWeight(vecPlotSyn);
				vecSynOriDiff = nan(size(vecSynW));
				for intSynIdx=1:numel(vecSynW)
					vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
					vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
				end
				
				intSubplot = intSubplot + 1;
				subplot(4,4,intSubplot);
				scatter(rad2deg(vecSynOriDiff)/2,vecSynW,20,'.');
				title(sprintf('%s %s to %s %s',strSourceAQ,strSourceEI,strTargetAQ,strTargetEI));
				xlabel('Orientation difference (degs)');
				ylabel('Synaptic weight');
				ylim([0 20]);
				xlim([-90 90]);
				fixfig;grid off;
			end
		end
	end
end
maxfig(gcf,0.85);
drawnow;
strFigFile1 = ['SynWeightByPrefOriDiff_' strFigBase];
export_fig(fullfile(strSourcePath,[strFigFile1 '.tif']));
export_fig(fullfile(strSourcePath,[strFigFile1 '.jpg']));
export_fig(fullfile(strSourcePath,[strFigFile1 '.pdf']));

%% plot 2
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	for intTargetAQ=[1 0]
		if intTargetAQ==1
			strTargetAQ = 'active';
		else
			strTargetAQ = 'silent';
		end
		
		
			for intSourceEI=1:2
				if intSourceEI==1
					strSourceEI = 'E';
				else
					strSourceEI = 'I';
				end
				for intSourceAQ=[1 0]
			if intSourceAQ==1
				strSourceAQ = 'active';
			else
				strSourceAQ = 'silent';
			end
				vecSourceN = find(intSourceAQ(:) == vecAQ & vecEI(:)==intSourceEI);
				vecTargetN = find(intTargetAQ(:) == vecAQ & vecEI(:)==intTargetEI);
				
				vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
				vecSynW = sLoad.vecSynWeight(vecPlotSyn);
				vecSynOriDiff = nan(size(vecSynW));
				for intSynIdx=1:numel(vecSynW)
					vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
					vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
				end
				
				vecCounts = histcounts(vecSynOriDiff,vecBinEdges);
				intSubplot = intSubplot + 1;
				
				subplot(4,4,intSubplot);
				
				plot(rad2deg(vecBinCenters)/2,vecCounts);
				title(sprintf('%s %s to %s %s',strSourceAQ,strSourceEI,strTargetAQ,strTargetEI));
				xlabel('Orientation difference (degs)');
				ylabel('# of synapses');
				%ylim([0 20]);
				xlim([-90 90]);
				fixfig;grid off;
			end
		end
	end
end
maxfig(gcf,0.85);
drawnow;
strFigFile2 = ['SynNumberByPrefOriDiff_' strFigBase];
export_fig(fullfile(strSourcePath,[strFigFile2 '.tif']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.jpg']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.pdf']));
