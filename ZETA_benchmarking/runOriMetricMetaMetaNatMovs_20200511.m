clear all;
%close all;
strPath = 'D:\Data\Processed\ZETA\NatMovs\';
strFigPath = 'D:\Data\Results\ZETA\NatMovs\';
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	'Poisson',...Area 3
	'Retina',...Area 4
	'CaNM',...Area 5
	'CaDG',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual'};%,...Area 14
%{
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area',...Area 18
	'Subiculum',...Area 19
	'Field CA1',...Area 20
	'Field CA2',...Area 21
	'Field CA3',...Area 22
	'Dentate gyrus',...Area 23
	'Retrosplenial'...Area 24
	};
%}
cellRunStim = {...
	'',...Stim 1 
	'RunDriftingGratings',...Stim 2 
	'RunNaturalMovie'...Stim 3
	};
cellRunRand = {...
	'',...Rand 1 
	'-Rand',...Rand 2 
	};
cellRepStr = {...
	'RunDriftingGratings','';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
};

%% prep
cellDatasetNames = {};
matZeta =[];
matNumCells = [];
matSignifZ = [];
matSignifHz = [];
matBinsAnovaFracP = [];
matBinsAnovaFracP_corr = [];
intIdx = 0;
for intArea=7:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if intArea < 7
		vecRunStims = 1;
	else
		vecRunStims = [3];
	end
	for intStimType=vecRunStims
		intIdx = intIdx + 1;
		strStim = cellRunStim{intStimType};
		strName = replace([strArea strStim],cellRepStr(:,1),cellRepStr(:,2));
		cellDatasetNames{intIdx} = strName;
	for intRandType=1:2
		%set var
		strRand = cellRunRand{intRandType};
		
		
		%% load data
		strRunType = [strArea strRand strStim];
		sDir=dir([strPath 'ZetaDataBins*' strRunType '*']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
			sLoad=load([strPath strFile]);
			if isfield(sLoad,'vecZeta')
				vecZP=sLoad.vecP;
				matNumCells(intIdx,intRandType) = numel(vecZP);
				matSignifZ(intIdx,intRandType) = sum(vecZP<0.05);
				matSignifHz(intIdx,intRandType) = sum(sLoad.vecHzP<0.05);
			else
				cellBinsAnovaP{intIdx,intRandType} = sLoad.matBinAnova;
				matBinsAnovaFracP(intIdx,intRandType,:) = sum(sLoad.matBinAnova<0.05,2);
				matBinsAnovaFracP_corr(intIdx,intRandType,:) = sum(sLoad.matBinAnova<(0.05/size(sLoad.matBinAnova,1)),2);
			end
		end
	end
	end
end
%%
indRem = any(matNumCells < 20,2);
matSignifZ(indRem,:) = [];
matNumCells(indRem,:) = [];
matSignifHz(indRem,:) = [];
matBinsAnovaFracP(indRem,:,:) = [];
matBinsAnovaFracP_corr(indRem,:,:) = [];
cellDatasetNames(indRem) = [];

[vecP_ZI,matCI_ZI] = binofit(matSignifZ(:,1),matNumCells(:,1),0.25);
[vecP_ZFA,matCI_ZFA] = binofit(matSignifZ(:,2),matNumCells(:,2),0.25);

[vecP_HzI,matCI_HzI] = binofit(matSignifHz(:,1),matNumCells(:,1),0.25);
[vecP_HzFA,matCI_HzFA] = binofit(matSignifHz(:,2),matNumCells(:,2),0.25);
figure
for intDataset=1:size(matSignifZ,1)
	matTable = zeros(2,2);
	matTable(1,1) = matSignifZ(intDataset,1);
	matTable(1,2) = matNumCells(intDataset,1)-matSignifZ(intDataset,1);
	matTable(2,1) = matSignifHz(intDataset,1);
	matTable(2,2) = matNumCells(intDataset,1)-matSignifHz(intDataset,1);
	[hI,pI] = fishertest(matTable);
	
	matTable = zeros(2,2);
	matTable(1,1) = matSignifZ(intDataset,2);
	matTable(1,2) = matNumCells(intDataset,2)-matSignifZ(intDataset,2);
	matTable(2,1) = matSignifHz(intDataset,2);
	matTable(2,2) = matNumCells(intDataset,2)-matSignifHz(intDataset,2);
	[hF,pF] = fishertest(matTable);
	
	subplot(4,6,intDataset)
	errorbar(2,vecP_ZI(intDataset),matCI_ZI(intDataset,1)-vecP_ZI(intDataset),matCI_ZI(intDataset,2)-vecP_ZI(intDataset),'bx')
	hold on
	errorbar(4,vecP_ZFA(intDataset),matCI_ZFA(intDataset,1)-vecP_ZFA(intDataset),matCI_ZFA(intDataset,2)-vecP_ZFA(intDataset),'bx')
	errorbar(1,vecP_HzI(intDataset),matCI_HzI(intDataset,1)-vecP_HzI(intDataset),matCI_HzI(intDataset,2)-vecP_HzI(intDataset),'kx')
	errorbar(3,vecP_HzFA(intDataset),matCI_HzFA(intDataset,1)-vecP_HzFA(intDataset),matCI_HzFA(intDataset,2)-vecP_HzFA(intDataset),'kx')
	hold off
	set(gca,'xtick',1:4,'xticklabel',{getGreek('Mu'),getGreek('Zeta'),['FA ' getGreek('Mu')],['FA ' getGreek('Zeta')]})
	%title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
	title(sprintf('%s; N=%d',cellDatasetNames{intDataset},matNumCells(intDataset,1)));
	ylabel('Fraction of significant cells')
	xlim([0.5 4.5])
	ylim([0 1])
	fixfig;
end

[h,pI]=ttest(vecP_ZI,vecP_HzI);
[h,pF]=ttest(vecP_ZFA,vecP_HzFA);
subplot(4,6,24)
plot(repmat([1 2],[intDataset 1])',[vecP_HzI vecP_ZI]','g')
hold on
plot(repmat([3 4],[intDataset 1])',[vecP_HzFA vecP_ZFA]','r')
title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
hold off	
xlim([0.5 4.5])
ylim([0 1])
set(gca,'xtick',1:4,'xticklabel',{getGreek('Mu'),getGreek('Zeta'),['FA ' getGreek('Mu')],['FA' getGreek('Zeta')]})
ylabel('Fraction of significant cells')
fixfig;
maxfig()
return
drawnow;
export_fig(sprintf('%sMetaSummaryFig.tif',strFigPath));
export_fig(sprintf('%sMetaSummaryFigEF.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sMetaSummaryFig.pdf',strFigPath));

%% make plot
vecBinDurs = sort([(2.^(-9:9))*(1/60)]);

matBinclude = bsxfun(@rdivide,matBinsAnovaFracP,matNumCells);
matBinclude_corr = bsxfun(@rdivide,matBinsAnovaFracP_corr,matNumCells);

matBinTP = squeeze(matBinclude(:,1,:));
matBinFP = squeeze(matBinclude(:,2,:));
matBinInclusion = matBinTP;% ./ (matBinTP + matBinFP); 

matBinInclusion_corr = squeeze(matBinclude_corr(:,1,:));

vecZTP = vecP_ZI;
vecZFP = vecP_ZFA;
vecZInclusion = vecZTP;% ./ (vecZTP + vecZFP); 

vecBinMean = mean(matBinInclusion);
vecBinSEM = std(matBinInclusion)./sqrt(size(matBinInclusion,1));


[h,vecP]=ttest(matBinInclusion,repmat(vecZInclusion,[1 size(matBinInclusion,2)]));
[h,h2,vecP] = fdr_bh(vecP);

[h,vecP_corr]=ttest(matBinInclusion_corr,repmat(vecZInclusion,[1 size(matBinInclusion_corr,2)]));
[h,h2,vecP_corr] = fdr_bh(vecP_corr);

vecBinMean_corr = mean(matBinInclusion_corr);
vecBinSEM_corr = std(matBinInclusion_corr)./sqrt(size(matBinInclusion,1));

figure;
subplot(2,2,1)
scatter(1/60,0.6,'og')
hold on
errorbar(vecBinDurs,vecBinMean_corr,vecBinSEM_corr,'kx')
errorbar([vecBinDurs(1) vecBinDurs(end)],mean(vecZInclusion)*[1 1],(std(vecZInclusion)/sqrt(size(matBinInclusion,1)))*[1 1],'bx-')
hold off
title('Bonferroni-corrected; blue=ZETA, black=binned ANOVA')
set(gca,'xscale','log')
xlabel('Bin size (s)')
ylim([0 0.7]);
ylabel('Inclusion rate')
fixfig;

subplot(2,2,2)
scatter(1/60,0.6,'og')
hold on
errorbar(vecBinDurs,vecBinMean,vecBinSEM,'kx')
errorbar([vecBinDurs(1) vecBinDurs(end)],mean(vecZInclusion)*[1 1],(std(vecZInclusion)/sqrt(size(matBinInclusion,1)))*[1 1],'bx-')
hold off
title('Uncorrected; blue=ZETA, black=binned ANOVA')
set(gca,'xscale','log')
xlabel('Bin size (s)')
ylim([0 0.7]);
ylabel('Inclusion rate')
fixfig;

subplot(2,2,3)
hold on
for intBin=1:numel(vecP_corr)
	text(vecBinDurs(intBin),0.7,sprintf('%.3f',vecP_corr(intBin)),'rotation',45);
text(vecBinDurs(intBin),0.3,sprintf('%.3e',vecP_corr(intBin)),'rotation',45);
end
hold off
xlim([min(vecBinDurs) max(vecBinDurs)]);
set(gca,'xscale','log')

subplot(2,2,4)
hold on
for intBin=1:numel(vecP_corr)
	text(vecBinDurs(intBin),0.7,sprintf('%.3f',vecP(intBin)),'rotation',45);
text(vecBinDurs(intBin),0.3,sprintf('%.3e',vecP(intBin)),'rotation',45);
end
hold off
xlim([min(vecBinDurs) max(vecBinDurs)]);
set(gca,'xscale','log')

%%
export_fig([strFigPath 'NatMovSummary.tif']);
export_fig([strFigPath 'NatMovSummary.pdf']);

%% relative
