clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
intResamps = 100;
strFileSearch = ['TsZetaSimTau*Resamp' num2str(intResamps) '.mat'];
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matTsZeta = sLoad.matTsZeta;
matZeta = sLoad.matZeta;
dblTau0 = 63/1000;%s
vecTau = sLoad.vecTau(2:end)*dblTau0;
strRec = sLoad.strRec;
intTauNum = numel(vecTau);

%% plot
hMegaFig = figure;maxfig;
vecR = nan(1,numel(vecTau));
matR_ci = nan(2,numel(vecTau));
vecAUC_Tau = nan(1,numel(vecTau));
matAUC_TauCi = nan(2,numel(vecTau));
vecRunTau = 1:numel(vecTau);
matTP = zeros(numel(vecTau),intResamps*2);
matFP = zeros(numel(vecTau),intResamps*2);
		
intRunTauNum = numel(vecRunTau);

for intTauIdx=1:numel(vecRunTau)
	intTau = vecRunTau(intTauIdx);
	dblTau=vecTau(intTau);
	
	%% plot
	h1 =subplot(3,intRunTauNum,intTauIdx);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matTsZeta(:,1,intTau) > 2 & matZeta(:,1) < 2) + 2*(matTsZeta(:,1,intTau) < 2 & matZeta(:,1) > 2) + 3*(matTsZeta(:,1,intTau) > 2 & matZeta(:,1) > 2);
	scatter(matTsZeta(:,1,intTau),matZeta(:,1),100,vecColor1,'.');
	colormap(h1,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	if intTauIdx==1
		xlabel('TS-ZETA (\zeta_c)')
		ylabel('ZETA (\zeta_c)')
		%legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
		fixfig;
	else
		fixfig;
		set(gca,'yticklabel',[]);
	end
	title(sprintf('TP=%.2f',sum(matTsZeta(:,1,intTau)>2)/numel(matTsZeta(:,1,intTau))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;
	xlim([0 8]);
	ylim([0 8]);
	
	h2=subplot(3,intRunTauNum,intRunTauNum+intTauIdx)
	vecColor2 = 1 + 1*(matTsZeta(:,2,intTau) < 2 & matZeta(:,2) > 2) + 2*(matTsZeta(:,2,intTau) > 2 & matZeta(:,2) < 2) + 3*(matTsZeta(:,2,intTau) > 2 & matZeta(:,2) > 2);
	scatter(matTsZeta(:,2,intTau),matZeta(:,2),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	if intTauIdx==1
		xlabel('TS-ZETA (\zeta_c)')
		ylabel('ZETA (\zeta_c)')
		%legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
		fixfig;
	else
		fixfig;
		set(gca,'yticklabel',[]);
	end
	title(sprintf('FP=%.2f',sum(matTsZeta(:,2,intTau)>2)/numel(matTsZeta(:,2,intTau))));
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	xlim([0 3]);
	ylim([0 3]);
	
	%% plot ROC
	cellColor = {'k',lines(1)};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	h2=subplot(3,intRunTauNum,intRunTauNum*2+intTauIdx)
	maxfig;
	hold on;
	
	for intTest=1:2
		if intTest == 1
			matDataZ = matTsZeta(:,:,intTau);
		else
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
		
		[dblAUC,Aci] = auc(cat(2,vecBothLabels,vecBothData));
		vecAUC(intTest) = dblAUC;
		
		if intTest == 1
			vecAUC_Tau(intTau) = dblAUC;
			matAUC_TauCi(:,intTau) = Aci;
			matTP(intTau,:) = vecTP;
			matFP(intTau,:) = vecFP;
		
		else
			dblAUC_zeta = dblAUC;
			vecAUC_zeta = Aci;
			vecTP_zeta = vecTP;
			vecFP_zeta = vecFP;
		end
	end
	hold off;
	title(sprintf('tau=%.3fs',dblTau));
	if intTauIdx==1
		xlabel('False positive fraction');
		ylabel('Inclusion fraction');
		legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
		fixfig;
	else
		fixfig;
		set(gca,'yticklabel',[]);
	end
	%% save correlation
	vecRandZeta = matZeta(:,1);
	vecRandTS_Zeta = matTsZeta(:,1,intTau);
	[R,P,RL,RU] = corrcoef(vecRandZeta,vecRandTS_Zeta);
	vecR(intTau) = R(1,2);
	matR_ci(:,intTau) = [RL(1,2) RU(1,2)];
end

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));

%% plot corr & auc
figure;maxfig;
subplot(2,3,1)
h=plot(matFP',matTP');
cMap = copper(size(matFP,1));
cMap = flipud(cMap);
colororder(gca,cMap);
colormap(cMap);
h=colorbar
h.Label.String='tau (s)';
h.Label.HorizontalAlignment='left';
h.Label.Position = [1 0.5 0];
set(h,'ytick',[0 1],'yticklabel',{sprintf('%.1e',vecTau(1)),sprintf('%.1e',vecTau(end))});
xlabel('False positive fraction');
ylabel('Inclusion fraction');
hold on
plot(vecFP_zeta,vecTP_zeta,'color',cellColor{2});
hold off
fixfig;

subplot(2,3,2)
hold on
errorbar(vecTau,vecR,vecR-matR_ci(1,:),vecR-matR_ci(2,:),'color',cellColor{1})
hold off
%ylim([0 1]);
set(gca,'xscale','log');
fixfig;
xlabel('Ca-indicator mean lifetime tau (s)');
ylabel('Pearson correlation, r(ZETA,TS-ZETA)');
grid off;
set(gca,'box','off');

subplot(2,3,3)
hold on
errorbar(vecTau,vecAUC_Tau,vecAUC_Tau-matAUC_TauCi(1,:),vecAUC_Tau-matAUC_TauCi(2,:),'color',cellColor{1})
plot(vecTau,dblAUC_zeta*ones(size(vecTau)),'color',cellColor{2})
hold off
set(gca,'xscale','log');
fixfig;
xlabel('Ca-indicator timescale (tau)');
ylabel('AUC TS-ZETA');
grid off;
set(gca,'box','off');

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.pdf']));
