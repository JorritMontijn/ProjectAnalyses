
%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
vecResamps = 100;%[100 200 500 1000 2000];
intResamps= numel(vecResamps);
%vecResamps = [5000 10000];
boolSave = true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
strFileSearch = ['Zeta2DataAnovaPoissonPeakResamp250.mat'];
sDir = dir(fullpath(strDataSourcePath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matAggTtest2 = sLoad.matTtest2;
matAggZeta2 = sLoad.matTsZeta2;
strRec = sLoad.strRec;

%% plot
intDiff = 1;
intResampIdx = 1;
matZeta2 = matAggZeta2(:,:,intResampIdx);
matTtest2 = matAggTtest2(:,:,intResampIdx);

figure;maxfig;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta2(:,1) > 2 & matTtest2(:,1) < 2) + 2*(matZeta2(:,1) < 2 & matTtest2(:,1) > 2) + 3*(matZeta2(:,1) > 2 & matTtest2(:,1) > 2);
scatter(matTtest2(:,1),matZeta2(:,1),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('T-test (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('Diff %d; TP at %s=0.05: %s=%.3f, %s=%.3f',intDiff,getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta2(:,1)>2)/numel(matZeta2(:,1)),'T2',sum(matTtest2(:,1)>2)/numel(matTtest2(:,1))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h2=subplot(2,3,2)
vecColor2 = 1 + 1*(matZeta2(:,2) < 2 & matTtest2(:,2) > 2) + 2*(matZeta2(:,2) > 2 & matTtest2(:,2) < 2) + 3*(matZeta2(:,2) > 2 & matTtest2(:,2) > 2);
scatter(matTtest2(:,2),matZeta2(:,2),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('T-test (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta2(:,2)>2)/numel(matZeta2(:,2)),'T2',sum(matTtest2(:,2)>2)/numel(matTtest2(:,2))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {lines(1),'k'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
h2=subplot(2,3,3)
maxfig;
hold on;

for intTest=1:2
	if intTest == 2
		matDataZ = matTtest2;
	else
		matDataZ = matZeta2;
	end
	matData = normcdf(matDataZ,'upper')*2;
	intCells = size(matData,1);
	vecBothData = cat(1,matData(:,1),matData(:,2));
	vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
	vecThresholds = sort(vecBothData);
	vecRealP = matData(:,1);
	vecShuffP = matData(:,2);
	
	vecTP = sum(vecRealP<vecThresholds',1)/intCells;
	vecFP = sum(vecShuffP<vecThresholds',1)/intCells;
	
	plot(vecFP,vecTP,'Color',cellColor{intTest});
	
	[dblAUC,Aci] = getAuc(vecTP,vecFP);
	vecAUC(intTest) = dblAUC;
end
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
legend({sprintf('ZETA-test, AUC=%.4f',vecAUC(1)),sprintf('T-test, AUC=%.4f',vecAUC(2))},'location','best','interpreter','none');
fixfig;


%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));
