clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

%% prep
intResamps = 250;
intNeurons = 100;
strFileSearch = ['TsZetaSim*Tau*Resamp' num2str(intResamps) '.mat'];
sDir = dir(fullpath(strDataPath,strFileSearch));
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
matTP = zeros(numel(vecTau),intNeurons*2);
matFP = zeros(numel(vecTau),intNeurons*2);
		
intRunTauNum = numel(vecRunTau);
figure;maxfig;
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
	else
		set(gca,'yticklabel',[]);
	end
	title(h1,sprintf('TP=%.2f',sum(matTsZeta(:,1,intTau)>2)/numel(matTsZeta(:,1,intTau))))
	%set(gca,'xscale','log','yscale','log');
	xlim(h1,[0 8]);
	ylim(h1,[0 8]);
	
	h2=subplot(3,intRunTauNum,intRunTauNum+intTauIdx)
	vecColor2 = 1 + 1*(matTsZeta(:,2,intTau) < 2 & matZeta(:,2) > 2) + 2*(matTsZeta(:,2,intTau) > 2 & matZeta(:,2) < 2) + 3*(matTsZeta(:,2,intTau) > 2 & matZeta(:,2) > 2);
	scatter(matTsZeta(:,2,intTau),matZeta(:,2),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	if intTauIdx==1
		xlabel('TS-ZETA (\zeta_c)')
		ylabel('ZETA (\zeta_c)')
		%legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
	else
		set(gca,'yticklabel',[]);
	end
	title(h2,sprintf('FP=%.2f',sum(matTsZeta(:,2,intTau)>2)/numel(matTsZeta(:,2,intTau))));
	%set(gca,'xscale','log','yscale','log');
	xlim(h2,[0 3]);
	ylim(h2,[0 3]);
	
	%% plot ROC
	cellColor = {'k',lines(1)};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	h3=subplot(3,intRunTauNum,intRunTauNum*2+intTauIdx)
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
		
		
		plot(h3,vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci] = getAuc(vecTP,vecFP);
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
	title(h3,sprintf('tau=%.3fs',dblTau));
	if intTauIdx==1
		xlabel(h3,'False positive fraction');
		ylabel(h3,'Inclusion fraction');
		legend(h3,{sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
	else
		set(h3,'yticklabel',[]);
	end
	%% save correlation
	vecRandZeta = matZeta(:,1);
	vecRandTS_Zeta = matTsZeta(:,1,intTau);
	[R,P,RL,RU] = corrcoef(vecRandZeta,vecRandTS_Zeta);
	
	R_reflective = sum(vecRandZeta.*vecRandTS_Zeta)/sqrt(sum(vecRandZeta.^2)*sum(vecRandTS_Zeta.^2));
	R_reflective_se = sqrt((1-R_reflective^2)/(intNeurons-2));
	vecR(intTau) = R_reflective;
	matR_ci(:,intTau) = [R_reflective-R_reflective_se R_reflective+R_reflective_se];
end
fixfig;

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

subplot(2,3,2)
hold on
errorbar(vecTau,vecR,vecR-matR_ci(1,:),vecR-matR_ci(2,:),'color',cellColor{1})
hold off
%ylim([0 1]);
set(gca,'xscale','log');
xlabel('Ca-indicator mean lifetime tau (s)');
ylabel('Reflective correlation, r(ZETA,TS-ZETA)');
grid off;
set(gca,'box','off');

subplot(2,3,3)
hold on
errorbar(vecTau,vecAUC_Tau,vecAUC_Tau-matAUC_TauCi(1,:),vecAUC_Tau-matAUC_TauCi(2,:),'color',cellColor{1})
plot(vecTau,dblAUC_zeta*ones(size(vecTau)),'color',cellColor{2})
hold off
set(gca,'xscale','log');
xlabel('Ca-indicator mean lifetime tau (s)');
ylabel('AUC TS-ZETA');
grid off;
set(gca,'box','off');
fixfig;

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Summary.pdf']));
