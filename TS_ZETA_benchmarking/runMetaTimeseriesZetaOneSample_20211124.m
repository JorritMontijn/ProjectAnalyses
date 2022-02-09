clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
%strFileSearch = ['SimOneSampleTsZetaQ0.mat'];
strFileSearch = ['SimOneSampleTsZetaQ0.mat'];
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matAggTtest = -norminv(sLoad.matTtest/2);
matAggZeta = -norminv(sLoad.matTsZeta/2);
matAggAnova = -norminv(sLoad.matAnova/2);
strRec = sLoad.strRec;

%% plot
intDiff = 1;
intResampIdx = 1;
matZeta = matAggZeta(:,:);
matTtest = matAggTtest(:,:);
matAnova = matAggAnova(:,:);

figure;maxfig;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta(:,1) > 2 & matTtest(:,1) < 2) + 2*(matZeta(:,1) < 2 & matTtest(:,1) > 2) + 3*(matZeta(:,1) > 2 & matTtest(:,1) > 2);
scatter(matTtest(:,1),matZeta(:,1),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('T-test (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('Diff %d; TP at %s=0.05: %s=%.3f, %s=%.3f',intDiff,getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta(:,1)>2)/numel(matZeta(:,1)),'T2',sum(matTtest(:,1)>2)/numel(matTtest(:,1))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h2=subplot(2,3,2);
vecColor2 = 1 + 1*(matZeta(:,2) < 2 & matTtest(:,2) > 2) + 2*(matZeta(:,2) > 2 & matTtest(:,2) < 2) + 3*(matZeta(:,2) > 2 & matTtest(:,2) > 2);
scatter(matTtest(:,2),matZeta(:,2),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('T-test (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta(:,2)>2)/numel(matZeta(:,2)),'T2',sum(matTtest(:,2)>2)/numel(matTtest(:,2))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h4 =subplot(2,3,4);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta(:,1) > 2 & matAnova(:,1) < 2) + 2*(matZeta(:,1) < 2 & matAnova(:,1) > 2) + 3*(matZeta(:,1) > 2 & matAnova(:,1) > 2);
scatter(matAnova(:,1),matZeta(:,1),100,vecColor1,'.');
colormap(h4,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('ANOVA (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('Diff %d; TP at %s=0.05: %s=%.3f, %s=%.3f',intDiff,getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta(:,1)>2)/numel(matZeta(:,1)),'A',sum(matAnova(:,1)>2)/numel(matAnova(:,1))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h5=subplot(2,3,5);
vecColor2 = 1 + 1*(matZeta(:,2) < 2 & matAnova(:,2) > 2) + 2*(matZeta(:,2) > 2 & matAnova(:,2) < 2) + 3*(matZeta(:,2) > 2 & matAnova(:,2) > 2);
scatter(matAnova(:,2),matZeta(:,2),100,vecColor1,'.');
colormap(h5,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('ANOVA (z)')
ylabel('ZETA (\delta\zeta_c)')
title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta(:,2)>2)/numel(matZeta(:,2)),'A',sum(matAnova(:,2)>2)/numel(matAnova(:,2))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {'k','r',lines(1)};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
h2=subplot(2,3,3)
maxfig;
hold on;

for intTest=1:3
	if intTest == 1
		matDataZ = matTtest;
	elseif intTest == 2
		matDataZ = matAnova;
	elseif intTest == 3
		matDataZ = matZeta;
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
	
	[dblAUC,Aci] = getAuc(vecShuffP,vecRealP);
	vecAUC(intTest) = dblAUC;
end
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
legend({sprintf('T-test, AUC=%.4f',vecAUC(1)),sprintf('ANOVA, AUC=%.4f',vecAUC(2)),sprintf('TS-ZETA-test, AUC=%.4f',vecAUC(3))},'location','best','interpreter','none');
fixfig;


%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));
