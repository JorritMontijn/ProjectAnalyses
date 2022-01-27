clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
intType=2;
if intType==0
%	strFileSearch = ['ZetaTimeSimResamp*100.mat'];
elseif intType==1
	strFileSearch = ['ZetaTimeSimResamp*250.mat'];
elseif intType==2
	strFileSearch = ['ZetaTimeSimBurstsResamp*250.mat'];
elseif intType==3
	strFileSearch = ['ZetaTimeSimBurstsOnlyResamp*250.mat'];
elseif intType==4
	strFileSearch = 'ZetaTimeV1*NpxResamp100.mat';
end
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));

matAnova = sLoad.matAnova;
matTtest = sLoad.matTtest;
matZeta = sLoad.matZeta;
matZetaTime = sLoad.matZetaTime;
matCompTimeAnova = sLoad.matCompTimeAnova;
matCompTimeZeta = sLoad.matCompTimeZeta;
matCompTimeZetaTime = sLoad.matCompTimeZetaTime;
intResampNum = sLoad.intResampNum;
strRec = sLoad.strRec;

%% plot
	matA_Z = -norminv(matAnova./2);
	matT_Z = -norminv(matTtest./2);
	matZeta_Z = -norminv(matZeta./2);
	matZetaT_Z = -norminv(matZetaTime(:,:,end)./2);
	matZetaT_Z(isinf(matZetaT_Z))=8;
	%matZetaT_Z(matZetaT_Z<0)=0;
	matTimeZeta = log(matCompTimeZeta);
	matTimeZetaT = log(matCompTimeZetaTime);
	matTimeAnova = log(matCompTimeAnova);
	
	
	figure;maxfig;
	h1 =subplot(2,3,1);
	vecColor1 = 1 + (matZetaT_Z(:,1) > 2 & matZeta_Z(:,1) < 2) + 2*(matZetaT_Z(:,1) < 2 & matZeta_Z(:,1) > 2) + 3*(matZetaT_Z(:,1) > 2 & matZeta_Z(:,1) > 2);
	matC = [0.5 0.5 0.5;...
			0 0.8 0;...
			0.8 0 0;...
			0 0 0.8];
	vecColsUnique = flat(unique(vecColor1))';
	colormap(h1,matC(1:max(vecColor1),:));
	hold on
	for intCol=vecColsUnique
		scatter(matZetaT_Z(vecColor1==intCol,1),matZeta_Z(vecColor1==intCol,1),100,matC(intCol,:),'.');
	end
	hold off
	colormap(h1,matC(1:max(vecColor1),:));
	vecLim = [0 max([get(gca,'xlim') get(gca,'ylim')])];
	xlim(vecLim);ylim(vecLim);
	xlabel('ZETA+ (IFR-perm)')
	ylabel('ZETA')
	title(sprintf('Resamp %d; TP at %s=0.05: %s=%.3f, %s=%.3f',intResampNum,getGreek('alpha'),[getGreek('zeta') '+'],sum(matZetaT_Z(:,1)>2)/numel(matZetaT_Z(:,1)),getGreek('zeta'),sum(matZeta_Z(:,1)>2)/numel(matZeta_Z(:,1))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;
	
	h2=subplot(2,3,2)
	vecColor2 = 1 + 1*(matZetaT_Z(:,2) < 2 & matZeta_Z(:,2) > 2) + 2*(matZetaT_Z(:,2) > 2 & matZeta_Z(:,2) < 2) + 3*(matZetaT_Z(:,2) > 2 & matZeta_Z(:,2) > 2);
	scatter(matZetaT_Z(:,2),matZeta_Z(:,2),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	vecLim = [0 max([get(gca,'xlim') get(gca,'ylim')])];
	xlim(vecLim);ylim(vecLim);
	xlabel('ZETA+ (IFR-perm)')
	ylabel('ZETA')
	title(sprintf('FP at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),[getGreek('zeta') '+'],sum(matZetaT_Z(:,2)>2)/numel(matZetaT_Z(:,2)),getGreek('zeta'),sum(matZeta_Z(:,2)>2)/numel(matZeta_Z(:,2))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	
	%% plot ROC
	cellColor = {'k','r',lines(1),[0.4 0 0.4]};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	h2=subplot(2,3,3)
	maxfig;
	hold on;
	
	for intTest=1:4
		if intTest == 1
			matDataZ = matT_Z;
		elseif intTest == 2
			matDataZ = matA_Z;
		elseif intTest == 3
			matDataZ = matZeta_Z;
		elseif intTest == 4
			matDataZ = matZetaT_Z;
		end
		%remove nans
		indRem=any(isnan(matDataZ),2);
		matDataZ(indRem,:) = [];
		matData = normcdf(matDataZ,'upper')*2;
		intCells = size(matData,1);
		vecBothData = cat(1,matData(:,1),matData(:,2));
		vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
		vecThresholds = sort(vecBothData);
		vecRealP = matData(:,1);
		vecShuffP = matData(:,2);
		
		vecTP = sum(vecRealP<=vecThresholds',1)/intCells;
		vecFP = sum(vecShuffP<=vecThresholds',1)/intCells;
		
		plot(vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci] = getAuc(vecShuffP,vecRealP);
		vecAUC(intTest) = dblAUC;
	end
	hold off;
	xlabel('False positive fraction');
	ylabel('Inclusion fraction');
	legend({sprintf('t-test, AUC=%.4f',vecAUC(1)),...
		sprintf('ANOVA, AUC=%.4f',vecAUC(2)),...
		sprintf('ZETA-test, AUC=%.4f',vecAUC(3)),...
		sprintf('ZETA+-test, AUC=%.4f',vecAUC(4))},'location','best','interpreter','none');
	fixfig;
	
	subplot(2,3,4)
	dblStep = 0.1;
	vecBinEdges=-4:dblStep:2;
	vecBinCenters = vecBinEdges(2:end) - dblStep/2;
	vecZetaTime = histcounts(matTimeZeta(:,1),vecBinEdges);
	vecZeta2Time = histcounts(matTimeZetaT(:,1),vecBinEdges);
	
	plot(vecBinCenters,vecZetaTime,'color',cellColor{1});
	hold on
	plot(vecBinCenters,vecZeta2Time,'color',cellColor{2});
	hold off
	xlabel('Computation time (log(s))');
	ylabel('Number of cells (count)');
	legend({'ZETA-test','ZETA+-test'},'location','best','interpreter','none');
	title(sprintf('Median log time, ZETA=%.3f, ZETA+ =%.3f',median(matTimeZeta(:,1)),median(matTimeZetaT(:,1))));
	fixfig;
	
	%% save
drawnow;
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.tif']));
export_fig(fullpath(strFigPath,[strrep(strFile,'.mat',''),'Overview.pdf']));


