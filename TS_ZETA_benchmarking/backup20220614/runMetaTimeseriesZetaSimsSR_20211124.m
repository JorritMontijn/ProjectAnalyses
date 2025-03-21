clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
intResamps = 250;
intSampFreq = 10;
strFileSearch = ['*SimNoSuperResolution*Resamp' num2str(intResamps) 'SampFreq' num2str(intSampFreq) '.mat'];
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matTtest = sLoad.matTtest;
matAnova = sLoad.matAnova;
matSrInTsZeta = sLoad.matSrInTsZeta;
matSrTsZeta = sLoad.matSrTsZeta;
matTsZeta = sLoad.matTsZeta;
matZeta = sLoad.matZeta;
dblTau0 = 63/1000;%s
vecTau = sLoad.vecTau*dblTau0;
strRec = sLoad.strRec;
intTauNum = numel(vecTau);

%% plot
hMegaFig = figure;maxfig;
vecR = nan(1,numel(vecTau));
matR_ci = nan(2,numel(vecTau));
vecAUC_Tau = nan(1,numel(vecTau));
matAUC_TauCi = nan(2,numel(vecTau));
vecRunTau = 1:numel(vecTau);
matTP = zeros(numel(vecTau),200);
matFP = zeros(numel(vecTau),200);

intRunTauNum = numel(vecRunTau);

for intTauIdx=1:numel(vecRunTau)
	intTau = vecRunTau(intTauIdx);
	dblTau=vecTau(intTau);
	intPlotTypes=3;
	intGraphs = 4;
	for intPlotType=1:intPlotTypes
		%% get data
		if intPlotType == 1
			matDat1 = matTsZeta;
			matDat2 = matSrTsZeta;
			strDat1 = 'TS-ZETA';
			strDat2 = 'SR-TS-ZETA';
		elseif intPlotType == 2
			matDat1 = matTsZeta;
			matDat2 = matSrInTsZeta;
			strDat1 = 'TS-ZETA';
			strDat2 = 'SR-In-TS-ZETA';
			
		elseif intPlotType == 3
			matDat1 = matSrTsZeta;
			matDat2 = matSrInTsZeta;
			strDat1 = 'SR-TS-ZETA';
			strDat2 = 'SR-In-TS-ZETA';
			
		end
		%% plot
		
		h1 =subplot(intPlotTypes,intGraphs,intPlotType);
		matC = [0.5 0.5 0.5;...
			0 0.8 0;...
			0.8 0 0;...
			0 0 0.8];
		vecColor1 = 1 + (matDat1(:,1,intTau) > 2 & matDat2(:,1) < 2) + 2*(matDat1(:,1,intTau) < 2 & matDat2(:,1) > 2) + 3*(matDat1(:,1,intTau) > 2 & matDat2(:,1) > 2);
		vecColsUnique = flat(unique(vecColor1))';
		colormap(h1,matC(1:max(vecColor1),:));
		hold on
		for intCol=vecColsUnique
			scatter(matDat1(vecColor1==intCol,1,intTau),matDat2(vecColor1==intCol,1),100,matC(intCol,:),'.');
		end
		hold off
		%xlim([0 1]);ylim([0 1]);
		if intTauIdx==1
			xlabel(strDat1)
			ylabel(strDat2)
			%legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
			fixfig;
		else
			fixfig;
			set(gca,'yticklabel',[]);
		end
		title(sprintf('TP=%.2f',sum(matDat1(:,1,intTau)>2)/numel(matDat1(:,1,intTau))))
		%set(gca,'xscale','log','yscale','log');
		fixfig;
		xlim([0 8]);
		ylim([0 8]);
		
		h2=subplot(intPlotTypes,intGraphs,intGraphs+intPlotType)
		vecColor2 = 1 + 1*(matDat1(:,2,intTau) < 2 & matDat2(:,2) > 2) + 2*(matDat1(:,2,intTau) > 2 & matDat2(:,2) < 2) + 3*(matDat1(:,2,intTau) > 2 & matDat2(:,2) > 2);
		hold on
		for intCol=vecColsUnique
			scatter(matDat1(vecColor1==intCol,2,intTau),matDat2(vecColor1==intCol,2),100,matC(intCol,:),'.');
		end
		hold off
		colormap(h2,matC(1:max(vecColor1),:));
		%xlim([0 1]);ylim([0 1]);
		if intTauIdx==1
			xlabel(strDat1)
			ylabel(strDat2)
			%legend({sprintf('TS-ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('ZETA-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
			fixfig;
		else
			fixfig;
			set(gca,'yticklabel',[]);
		end
		title(sprintf('%s,FP=%.3f',strDat1,sum(matDat1(:,2,intTau)>2)/numel(matDat1(:,2,intTau))));
		%set(gca,'xscale','log','yscale','log');
		fixfig;maxfig;
		xlim([0 3]);
		ylim([0 3]);
		
		%% plot ROC
		cellColor = {'k',lines(1)};
		%vecH(intResampNpx) = subplot(4,3,intResampNpx);
		h3=subplot(intPlotTypes,intGraphs,intGraphs*2+intPlotType)
		maxfig;
		hold on;
		
		for intTest=1:2
			if intTest == 1
				matDataZ = matDat1(:,:,intTau);
			else
				matDataZ = matDat2;
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
			
			[dblAUC,Aci] = getAuc(matDataZ(:,1),matDataZ(:,2));
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
			legend({sprintf('%s, AUC=%.2f',strDat1,vecAUC(1)),sprintf('%s, AUC=%.2f',strDat2,vecAUC(2))},'location','best','interpreter','none');
			fixfig;
		else
			fixfig;
			set(gca,'yticklabel',[]);
		end
		
	end
	
	%% calc AUC t-test & anova
	cellColor = {'k','r',[0 0 1],[0.5 0 0.5],lines(1)};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	h4=subplot(intPlotTypes,intGraphs,intPlotTypes*intGraphs);
	maxfig;
	hold on;
	for intTest=1:5
		if intTest == 1
			matDataZ = -norminv(matTtest/2);
			strLeg = 'T-test';
		elseif intTest == 2
			matDataZ = -norminv(matAnova/2);
			strLeg = 'Anova';
		elseif intTest == 3
			matDataZ = matTsZeta;
			strLeg = 'TS-Z';
		elseif intTest == 4
			matDataZ = matSrTsZeta;
			strLeg = 'SR-TS-Z';
		elseif intTest == 5
			matDataZ = matSrInTsZeta;
			strLeg = 'In-TS-Z';
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
		
		[dblAUC,Aci] = getAuc(matDataZ(:,1),matDataZ(:,2));
		vecAUC(intTest) = dblAUC;
		cellLegend{intTest} = sprintf('%s, AUC=%.2f',strLeg,dblAUC);
	end
	xlabel('False positive fraction');
	ylabel('Inclusion fraction');
	legend(cellLegend,'location','best','interpreter','none');
	fixfig;
	
	
	%% plot difference
	matDiff_TS_SRTS = matSrTsZeta-matTsZeta;
	matDiff_TS_SRINTS = matSrInTsZeta-matTsZeta;
	matDiff_SRTS_SRINTS = matSrInTsZeta-matSrTsZeta;
	
	h4=subplot(intPlotTypes,intGraphs,intGraphs);
	hold on
	bplot(matDiff_TS_SRTS(:,1),1);
	bplot(matDiff_TS_SRINTS(:,1),2);
	bplot(matDiff_SRTS_SRINTS(:,1),3);
	hold off
	set(gca,'xtick',1:3,'xticklabel',{'SR v base','SRI v base','SRI v SR'});
	ylabel('Z-score change over base');
	title('TP');
	fixfig;%grid off;
	vecLimY = [-1 1]*max(abs(get(gca,'ylim')));
	ylim(vecLimY);
	xtickangle(20);
	
	h4=subplot(intPlotTypes,intGraphs,intGraphs*2);
	hold on
	bplot(matDiff_TS_SRTS(:,2),1);
	bplot(matDiff_TS_SRINTS(:,2),2);
	bplot(matDiff_SRTS_SRINTS(:,2),3);
	hold off
	set(gca,'xtick',1:3,'xticklabel',{'SR v base','SRI v base','SRI v SR'});
	ylabel('Z-score change over base');
	title('FP');
	fixfig;%grid off;
	ylim(vecLimY);
	xtickangle(20);
	
	
	maxfig;
	
	
	%% save correlation
	vecRandZeta = matDat2(:,1);
	vecRandTS_Zeta = matDat1(:,1,intTau);
	[R,P,RL,RU] = corrcoef(vecRandZeta,vecRandTS_Zeta);
	vecR(intTau) = R(1,2);
	matR_ci(:,intTau) = [RL(1,2) RU(1,2)];
end

%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));
return
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
