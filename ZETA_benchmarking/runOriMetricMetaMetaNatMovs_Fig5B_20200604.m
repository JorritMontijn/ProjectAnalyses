clear all;
%close all;
strDisk = 'F:';
strPath = [strDisk '\Data\Processed\ZETA\NatMovs\'];
strFigPath = [strDisk '\Data\Results\ZETA\NatMovs\'];
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
		sDir=dir([strPath 'ZetaData*' strRunType '*']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
			sLoad=load([strPath strFile]);
			if isfield(sLoad,'vecZeta')
				vecZP=sLoad.vecP;
				cellAreasZ{intIdx} = sLoad.cellArea{1};
				matNumCells(intIdx,intRandType) = numel(vecZP);
				cellZetaP{intIdx,intRandType} = vecZP;
				matSignifZ(intIdx,intRandType) = sum(vecZP<0.05);
				matSignifHz(intIdx,intRandType) = sum(sLoad.vecHzP<0.05);
			elseif isfield(sLoad,'vecZetaP')
				vecZP=sLoad.vecZetaP;
				matNumCells(intIdx,intRandType) = numel(vecZP);
				matSignifZ(intIdx,intRandType) = sum(vecZP<0.05);
				cellBinsAnovaP{intIdx,intRandType} = sLoad.matBinAnova;
				matBinsAnovaFracP(intIdx,intRandType,:) = sum(sLoad.matBinAnova<0.05,2);
				matBinsAnovaFracP_corr(intIdx,intRandType,:) = sum(sLoad.matBinAnova<(0.05/size(sLoad.matBinAnova,1)),2);
			elseif isfield(sLoad,'matBinAnova')
				%cellAreasA{intIdx} = sLoad.cellArea{1};
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
cellBinsAnovaP(indRem,:) = [];
cellZetaP(indRem,:) = [];
matBinsAnovaFracP(indRem,:,:) = [];
matBinsAnovaFracP_corr(indRem,:,:) = [];
cellDatasetNames(indRem) = [];

[vecP_ZI,matCI_ZI] = binofit(matSignifZ(:,1),matNumCells(:,1),0.25);
[vecP_ZFA,matCI_ZFA] = binofit(matSignifZ(:,2),matNumCells(:,2),0.25);

%% plot ROC
vecTotCells = sum(matNumCells,1);
vecAlpha=0.0001:0.0001:1;
vecAlpha=[10.^[-10:-5] 1.01.^[-1000:0]];
vecTP_Z = nan(1,numel(vecAlpha));
vecTP_A = nan(1,numel(vecAlpha));
vecFP_Z = nan(1,numel(vecAlpha));
vecFP_A = nan(1,numel(vecAlpha));
for intAlpha=1:numel(vecAlpha)
	dblAlpha=vecAlpha(intAlpha);
	
	%bins signifcalt
intBins = size(cellBinsAnovaP{1},1);
vecIncl = any(matBinsAnovaFracP_corr,3);
y=cellfun(@(x) any(x<dblAlpha),cellBinsAnovaP,'uniformoutput',false);
matInclA=cellfun(@(x) sum(x)/numel(x),y);
matInclA(cellfun(@isempty,cellBinsAnovaP)) = nan;

%zeta signicfalt
matInclZ = cellfun(@(x) sum(x<dblAlpha)/numel(x),cellZetaP);
matInclZ(cellfun(@isempty,cellZetaP)) = nan;

%fisher
dblZ_TP = sum(matInclZ(:,1) .* matNumCells(:,1));
dblZ_FP = sum(matInclZ(:,2) .* matNumCells(:,2));

dblA_TP = sum(matInclA(:,1) .* matNumCells(:,1));
dblA_FP = sum(matInclA(:,2) .* matNumCells(:,2));

mat2x2 = [dblZ_TP dblA_TP; dblZ_FP dblA_FP];
%[p,chi2stat] = chi2test(mat2x2)

% [h, p, stats] = fishertest(mat2x2)

 vecTP = [dblZ_TP vecTotCells(1) dblA_TP vecTotCells(2)];
pTP=bino2test(vecTP);

 vecFP = [dblZ_FP vecTotCells(1) dblA_FP vecTotCells(2)];
pFP=bino2test(vecFP);

vecTP_Z(intAlpha) = dblZ_TP/vecTotCells(1);
vecTP_A(intAlpha) = dblA_TP/vecTotCells(1);
vecFP_Z(intAlpha) = dblZ_FP/vecTotCells(2);
vecFP_A(intAlpha) = dblA_FP/vecTotCells(2);
end
% plot
 plot(vecFP_Z,vecTP_Z)
 hold on
  plot(vecFP_A,vecTP_A,'r')
 hold off
xlim([0 1]);ylim([0 1])
xlabel('False positive rate');
ylabel('Inclusion rate');
fixfig;

dblStepAUC = 0.01;
vecAUC_edges = dblStepAUC:dblStepAUC:1;
vecAUC_centers = vecAUC_edges(2:end) - dblStepAUC/2;
[vecCountsZ,vecMeansZ] = makeBins(vecFP_Z,vecTP_Z,vecAUC_edges);
vecAUC_Z = interp1(vecAUC_centers(~isnan(vecMeansZ)),vecMeansZ(~isnan(vecMeansZ)),vecAUC_centers);
dblAUC_Z = sum(vecAUC_Z(~isnan(vecAUC_Z)))/sum(~isnan(vecAUC_Z));
[vecCountsA,vecMeansA] = makeBins(vecFP_A,vecTP_A,vecAUC_edges);
vecAUC_A = interp1(vecAUC_centers(~isnan(vecMeansA)),vecMeansA(~isnan(vecMeansA)),vecAUC_centers);
dblAUC_A = sum(vecAUC_A(~isnan(vecAUC_A)))/sum(~isnan(vecAUC_A));


%% statistical test AUC
x=cellfun(@(x) min(x,[],1),cellBinsAnovaP,'uniformoutput',false);
vecP_TP_A=cell2vec(x(:,1));
vecP_FP_A=cell2vec(x(:,2));

vecP_TP_Z=cell2vec(cellZetaP(:,1));
vecP_FP_Z=cell2vec(cellZetaP(:,2));

vecAllA = cat(1,vecP_TP_A,vecP_FP_A);
vecAllZ = cat(1,vecP_TP_Z,vecP_FP_Z);

vecClass = cat(1,0*vecP_TP_A,0*vecP_FP_A+1);
[AUC_A,AUC_A_ci,AUC_A_se] = auc([vecClass vecAllA],0.05,'mann-whitney');
[AUC_Z,AUC_Z_ci,AUC_Z_se] = auc([vecClass vecAllZ],0.05,'mann-whitney');

% Observed data
m0 = AUC_A - AUC_Z;
s0 = (AUC_A_se + AUC_Z_se)/2;
z = m0/s0;
AUC_p = 1 - abs(normcdf(z)-normcdf(-z));

title(sprintf('Z=%.3f/%.3f(B);A=%.3f/sd=%.3f (R); p=%.3f',AUC_Z,AUC_Z_se,AUC_A,AUC_A_se,AUC_p));
drawnow;
%%
export_fig([strFigPath 'Fig5C_ROC.tif']);
export_fig([strFigPath 'Fig5C_ROC.pdf']);

%% make plot
vecBinDurs = sort([(2.^(-9:9))*(1/60)]);
%vecBinDurs = sort([(2.^(0:9))*(1/60)]);

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


[h,vecP1]=ttest(matBinInclusion,repmat(vecZInclusion,[1 size(matBinInclusion,2)]));
[h,h2,vecP] = fdr_bh(vecP1);

[h,vecP1_corr]=ttest(matBinInclusion_corr,repmat(vecZInclusion,[1 size(matBinInclusion_corr,2)]));
[h,h2,vecP_corr] = fdr_bh(vecP1_corr);

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
