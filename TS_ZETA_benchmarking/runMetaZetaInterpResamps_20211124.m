clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
intResamps = 10000;
strFileSearch = ['Zeta2SimZeta2*Resamp' num2str(intResamps) '.mat'];
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matZeta2 = sLoad.matZeta2;
matZeta = sLoad.matZeta;
matZetaTime = sLoad.matZeta;
matZeta2Time = sLoad.matZeta2Time;
strRec = sLoad.strRec;

hMegaFig = figure;maxfig;

%% plot
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta2(:,1) > 2 & matZeta(:,1) < 2) + 2*(matZeta2(:,1) < 2 & matZeta(:,1) > 2) + 3*(matZeta2(:,1) > 2 & matZeta(:,1) > 2);
scatter(matZeta2(:,1),matZeta(:,1),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('ZETA2 (no interp) (\zeta_c)')
ylabel('ZETA (no interp) (\zeta_c)')
title(sprintf('TP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta2(:,1)>2)/numel(matZeta2(:,1)),getGreek('zeta'),sum(matZeta(:,1)>2)/numel(matZeta(:,1))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h2=subplot(2,3,2)
vecColor2 = 1 + 1*(matZeta2(:,2) < 2 & matZeta(:,2) > 2) + 2*(matZeta2(:,2) > 2 & matZeta(:,2) < 2) + 3*(matZeta2(:,2) > 2 & matZeta(:,2) > 2);
scatter(matZeta2(:,2),matZeta(:,2),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('ZETA2 (no interp) (\zeta_c)')
ylabel('ZETA (no interp) (\zeta_c)')
title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '2'],sum(matZeta2(:,2)>2)/numel(matZeta2(:,2)),getGreek('zeta'),sum(matZeta(:,2)>2)/numel(matZeta(:,2))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {lines(1),'k'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
h2=subplot(2,3,3)
maxfig;
hold on;

for intTest=1:2
	if intTest == 1
		matDataZ = matZeta;
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
	
	[dblAUC,Aci] = auc(cat(2,vecBothLabels,vecBothData));
	vecAUC(intTest) = dblAUC;
end
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
legend({sprintf('ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA2-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
fixfig;

subplot(2,3,4)
dblStep = 0.5;
vecBinEdges=0:dblStep:10;
vecBinCenters = vecBinEdges(2:end) - dblStep/2;
matZetaTime = sLoad.matZeta;
matZeta2Time = sLoad.matZeta2Time;
vecZetaTime = histcounts(matZetaTime(:,1),vecBinEdges);
vecZeta2Time = histcounts(matZeta2Time(:,1),vecBinEdges);

plot(vecBinCenters,vecZetaTime,'color',cellColor{1});
hold on
plot(vecBinCenters,vecZeta2Time,'color',cellColor{2});
hold off
xlabel('Computation time (s)');
ylabel('Number of cells (count)');
legend({'ZETA-test','ZETA2-test'},'location','best','interpreter','none');
title(sprintf('Median time, ZETA=%.3fs, ZETA2 (no interp)=%.3fs',median(matZetaTime(:,1)),median(matZeta2Time(:,1))));
fixfig;


%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));
