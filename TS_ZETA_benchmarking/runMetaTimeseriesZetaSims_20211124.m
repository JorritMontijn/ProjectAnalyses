clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
strFileSearch = '*SimSampFreqs*.mat';
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matTsZeta = sLoad.matTsZeta;
matZeta = sLoad.matZeta;
strRec = sLoad.strRec;
vecSampFreqs = sLoad.vecSampFreqs;

hMegaFig = figure;maxfig;
vecR = nan(1,numel(vecSampFreqs));
matR_ci = nan(2,numel(vecSampFreqs));
for intSampIdx=1:numel(vecSampFreqs)
	dblSampFreq=vecSampFreqs(intSampIdx);
	
	%% plot
	h1 =subplot(3,4,intSampIdx);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matTsZeta(:,1,intSampIdx) > 2 & matZeta(:,1,intSampIdx) < 2) + 2*(matTsZeta(:,1,intSampIdx) < 2 & matZeta(:,1,intSampIdx) > 2) + 3*(matTsZeta(:,1,intSampIdx) > 2 & matZeta(:,1,intSampIdx) > 2);
	scatter(matTsZeta(:,1,intSampIdx),matZeta(:,1,intSampIdx),100,vecColor1,'.');
	colormap(h1,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('TS-ZETA (\zeta_c)')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('Incl. at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),['TS-' getGreek('zeta')],sum(matTsZeta(:,1,intSampIdx)>2)/numel(matTsZeta(:,1,intSampIdx)),getGreek('zeta'),sum(matZeta(:,1,intSampIdx)>2)/numel(matZeta(:,1,intSampIdx))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;
	
	h2=subplot(3,4,4+intSampIdx)
	vecColor2 = 1 + 1*(matTsZeta(:,2,intSampIdx) < 2 & matZeta(:,2,intSampIdx) > 2) + 2*(matTsZeta(:,2,intSampIdx) > 2 & matZeta(:,2,intSampIdx) < 2) + 3*(matTsZeta(:,2,intSampIdx) > 2 & matZeta(:,2,intSampIdx) > 2);
	scatter(matTsZeta(:,2,intSampIdx),matZeta(:,2,intSampIdx),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('TS-ZETA (\zeta_c)')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),['TS-' getGreek('zeta')],sum(matTsZeta(:,2,intSampIdx)>2)/numel(matTsZeta(:,2,intSampIdx)),getGreek('zeta'),sum(matZeta(:,2,intSampIdx)>2)/numel(matZeta(:,2,intSampIdx))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	
	%% plot ROC
	cellColor = {lines(1),'k'};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	h2=subplot(3,4,8+intSampIdx)
	maxfig;
	hold on;
	
	for intTest=1:2
		if intTest == 1
			matDataZ = matTsZeta;
		else
			matDataZ = matZeta;
		end
		matData = normcdf(matDataZ,'upper')*2;
		intCells = size(matData,1);
		vecBothData = cat(1,matData(:,1,intSampIdx),matData(:,2,intSampIdx));
		vecBothLabels = cat(1,zeros(size(matData(:,1,intSampIdx))),ones(size(matData(:,1,intSampIdx))));
		vecThresholds = sort(vecBothData);
		vecRealP = matData(:,1,intSampIdx);
		vecShuffP = matData(:,2,intSampIdx);
		
		vecTP = sum(vecRealP<vecThresholds',1)/intCells;
		vecFP = sum(vecShuffP<vecThresholds',1)/intCells;
		
		plot(vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci] = getAuc(vecTP,vecFP);
		vecAUC(intTest) = dblAUC;
	end
	hold off;
	xlabel('False positive fraction');
	ylabel('Inclusion fraction');
	title(sprintf('%d Hz; ROC analysis',dblSampFreq));
	legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
	fixfig;
	
	%% save correlation
	vecRandZeta = matZeta(:,2,intSampIdx);
	vecRandTS_Zeta = matTsZeta(:,2,intSampIdx);
	[R,P,RL,RU] = corrcoef(vecRandZeta,vecRandTS_Zeta);
	vecR(intSampIdx) = R(1,2);
	matR_ci(:,intSampIdx) = [RL(1,2) RU(1,2)];
end

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));

%% plot corr
figure
errorbar(vecSampFreqs,vecR,vecR-matR_ci(1,:),vecR-matR_ci(2,:))
ylim([0 1]);
set(gca,'xscale','log');
fixfig;
xlabel('Sampling frequency (Hz)');
ylabel('Pearson correlation, r(ZETA,TS-ZETA)');
grid off;
set(gca,'box','off');
xlim([0.9 101]);

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.pdf']));
