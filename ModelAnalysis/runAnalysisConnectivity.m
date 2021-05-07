strPath = 'F:\Code\Simulations\SimulationsEVS\Connectivity\';
%strSearchFile = '*LargeRetOriTuning*';
strSearchFile = '*2021-01-08*';
%strSearchFile = '*Conn*2020-10-29*';
%strSearchFile = '*2021-03-22*';
%strSearchFile = '*2018-07-19*';
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
sConn = sLoad.sConnectivity;
vecEI = sConn.vecCellTypes;
%% prep
vecPrefOri = sConn.vecPrefOri*2;
vecUniques = unique(roundi(diff(sort(vecPrefOri)),6));
dblStep = vecUniques(2);
vecBinEdges = (-pi-(dblStep/2)):dblStep:(pi+(dblStep/2));
vecBinCenters = vecBinEdges(2:end)-(dblStep/2);
matSynFromTo = sConn.matSynFromTo;

%% plot 0
if ~isfield(sConn,'matTuningSimilarity')
	matPrefGabors = sConn.matPrefGabors;
	intCellsV1 = size(matPrefGabors,3);
	matTuningSimilarity=nan(intCellsV1,intCellsV1);
	parfor i1=1:intCellsV1
		a = abs(matPrefGabors(:,:,i1));
		a = a - (sum(a(:),'double') / numel(a));
		for i2=1:intCellsV1
			b = abs(matPrefGabors(:,:,i2)); %#ok<PFBNS>
			b = b - (sum(b(:),'double') / numel(b));
			r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
			matTuningSimilarity(i1,i2) = r;
		end
	end
	matTuningSimilarity(diag(diag(true(intCellsV1)))) = 0;
else
	matTuningSimilarity = sConn.matTuningSimilarity;
end
subplot(2,3,2);
imagesc(matTuningSimilarity);colorbar
xlabel('Neuron #');
ylabel('Neuron #');
%title('Tuning similarity 2021-03-22');
title('Tuning similarity 2021-01-08');

fixfig;grid off


strFigFile1 = ['SynWeightByPrefOriDiff_' strFigBase];
return
%% plot
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	for intSourceEI=1:2
		if intSourceEI==1
			strSourceEI = 'E';
		else
			strSourceEI = 'I';
		end
		
		vecSourceN = find(vecEI(:)==intSourceEI);
		vecTargetN = find(vecEI(:)==intTargetEI);
		
		vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
		vecSynW = sConn.vecSynWeight(vecPlotSyn);
		vecSynOriDiff = nan(size(vecSynW));
		for intSynIdx=1:numel(vecSynW)
			vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
			vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
		end
		
		intSubplot = intSubplot + 1;
		subplot(2,2,intSubplot);
		scatter(rad2deg(vecSynOriDiff)/2,vecSynW,20,'.');
		title(sprintf('%s %s to %s %s',strSourceAQ,strSourceEI,strTargetAQ,strTargetEI));
		xlabel('Orientation difference (degs)');
		ylabel('Synaptic weight');
		ylim([0 20]);
		xlim([-90 90]);
		fixfig;grid off;
	end
end

maxfig(gcf,0.85);
drawnow;
strFigFile1 = ['SynWeightByPrefOriDiff_' strFigBase];
%export_fig(fullfile(strSourcePath,[strFigFile1 '.tif']));
%export_fig(fullfile(strSourcePath,[strFigFile1 '.jpg']));
%export_fig(fullfile(strSourcePath,[strFigFile1 '.pdf']));

%% plot 2
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	
	for intSourceEI=1:2
		if intSourceEI==1
			strSourceEI = 'E';
		else
			strSourceEI = 'I';
		end
		
		vecSourceN = find(vecEI(:)==intSourceEI);
		vecTargetN = find(vecEI(:)==intTargetEI);
		
		vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
		vecSynW = sConn.vecSynWeight(vecPlotSyn);
		vecSynOriDiff = nan(size(vecSynW));
		for intSynIdx=1:numel(vecSynW)
			vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
			vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
		end
		
		vecCounts = histcounts(vecSynOriDiff,vecBinEdges);
		intSubplot = intSubplot + 1;
		
		subplot(2,2,intSubplot);
		
		plot(rad2deg(vecBinCenters)/2,vecCounts);
		title(sprintf('%s %s to %s %s',strSourceAQ,strSourceEI,strTargetAQ,strTargetEI));
		xlabel('Orientation difference (degs)');
		ylabel('# of synapses');
		%ylim([0 20]);
		xlim([-90 90]);
		fixfig;grid off;
	end
end
maxfig(gcf,0.85);
drawnow;
strFigFile2 = ['SynNumberByPrefOriDiff_' strFigBase];
export_fig(fullfile(strSourcePath,[strFigFile2 '.tif']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.jpg']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.pdf']));
