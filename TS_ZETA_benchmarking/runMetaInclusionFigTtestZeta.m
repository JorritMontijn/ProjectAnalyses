%% file locs
clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

strIndicator = 'GCaMP';

%% load
dblAlpha = 0.05;
cellTests = {'Spike times','Time-series','Two-sample spike t','Two-sample t-s'};
cellTests = {'Neuropixels','GCaMP6f','Two-sample Npx','Two-sample GCaMP6f'};
figure;hAx0=axes;hold on;plot([dblAlpha dblAlpha],[0 1],'--','color',[0.5 0.5 0.5]);
hFig=figure;maxfig;
vecBinoP = nan(1,4);
for intCase=1:4
	if intCase==1
		%zeta
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visualRunDriftingGratingsResamp500.mat');
		sLoad = load(strFile);
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visual-RandRunDriftingGratingsResamp500.mat');
		sLoadRand = load(strFile);
		matMeanP = cat(1,sLoad.vecTtestP,sLoadRand.vecTtestP);
		matZetaP = cat(1,sLoad.vecZetaP_UniStitch,sLoadRand.vecZetaP_UniStitch); %with replacement
		matAnovaP = cat(1,sLoad.vecAnovaP_optimal,sLoadRand.vecAnovaP_optimal);
	elseif intCase==2
		%ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta(strIndicator);
	elseif intCase==3
		%2-sample zeta
		%strFile = fullpath(strDataPath, 'Zeta2DataAnovaV1RunDriftingGratingsResamp500.mat');
		strFile = fullpath(strDataPath, 'Zeta2DataStimDiffV1RunDriftingGratingsResamp500Q0.mat');
		%strFile = fullpath(strDataPath, 'Zeta2DataShiftRespV1RunDriftingGratingsResamp500Q0.mat');
		
		sLoad = load(strFile);
		matMeanP = sLoad.matTtest2';
		matZetaP = sLoad.matZeta2';
		matAnovaP = sLoad.matAnova2_optimal';
	elseif intCase==4
		%2-sample ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta2('');
	end
	
	%remove nans
	indRem = any(isnan(matMeanP) | isnan(matZetaP) | isnan(matAnovaP),1);
	matMeanP(:,indRem)=[];
	matZetaP(:,indRem)=[];
	matAnovaP(:,indRem)=[];
	
	intN = size(matZetaP,2);
	[dblTPRZ,vecTPRZci] = binofit(sum(matZetaP(1,:)<dblAlpha),intN);
	[dblTPRT,vecTPRTci] = binofit(sum(matMeanP(1,:)<dblAlpha),intN);
	[dblTPRA,vecTPRAci] = binofit(sum(matAnovaP(1,:)<dblAlpha),intN);
	
	[dblFPRZ,vecFPRZci] = binofit(sum(matZetaP(2,:)<dblAlpha),intN);
	[dblFPRT,vecFPRTci] = binofit(sum(matMeanP(2,:)<dblAlpha),intN);
	[dblFPRA,vecFPRAci] = binofit(sum(matAnovaP(2,:)<dblAlpha),intN);
	
	dblSensitivityZ = dblTPRZ/dblFPRZ;
	dblSensitivityT = dblTPRT/dblFPRT;
	dblSensitivityA = dblTPRA/dblFPRA;
	
	vecBinoP(intCase) = bino2test(sum(matMeanP(1,:)<dblAlpha),intN,sum(matZetaP(1,:)<dblAlpha),intN);

	%% plot
	figure(hFig)
	vecColZ = lines(1);
	vecColT = [0 0 0];
	vecColA = [0.8 0 0];
	subplot(2,4,intCase);
	hold on
	bar(1,100*dblTPRT,0.5,'Facecolor',vecColT,'Edgecolor','none');
	bar(2,100*dblTPRZ,0.5,'Facecolor',vecColZ,'Edgecolor','none');
	bar(3,100*dblTPRA,0.5,'Facecolor',vecColA,'Edgecolor','none');
	errorbar(1,100*dblTPRT,...
		100*[dblTPRT-vecTPRTci(1)],...
		100*[dblTPRT-vecTPRTci(2)],'linestyle','none','color',vecColT,'capsize',20);
	errorbar(2,100*[dblTPRZ],...
		100*[dblTPRZ-vecTPRZci(1)],...
		100*[dblTPRZ-vecTPRZci(2)],'linestyle','none','color',vecColZ,'capsize',20);
	errorbar(3,100*[dblTPRA],...
		100*[dblTPRA-vecTPRAci(1)],...
		100*[dblTPRA-vecTPRAci(2)],'linestyle','none','color',vecColA,'capsize',20);
	xlim([0 4]);
	ylim([0 100]);
	title(sprintf('%s, +%d%%, +%d%%',cellTests{intCase},round((dblTPRZ/dblTPRT - 1)*100),round((dblTPRZ/dblTPRA - 1)*100)));
	ylabel('Cell inclusion rate (%)');%,'Color',[0 0.7 0]);
	ax = gca;
	set(ax,'xtick',[1 2 3]);
	ax.XTickLabel{1} = '\color{black}T-test';
	ax.XTickLabel{2} = '\color[rgb]{0,0.4470,0.7410}ZETA';
	ax.XTickLabel{3} = '\color[rgb]{0.8,0,0}ANOVA';
	
	%%{
	subplot(2,4,4+intCase);
	hold on
	bar(1,100*dblFPRT,0.5,'Facecolor',vecColT,'Edgecolor','none');
	bar(2,100*dblFPRZ,0.5,'Facecolor',vecColZ,'Edgecolor','none');
	bar(3,100*dblFPRA,0.5,'Facecolor',vecColA,'Edgecolor','none');
	errorbar(1,100*dblFPRT,...
		100*[dblFPRT-vecFPRTci(1)],...
		100*[dblFPRT-vecFPRTci(2)],'linestyle','none','color',vecColT,'capsize',20);
	errorbar(2,100*[dblFPRZ],...
		100*[dblFPRZ-vecFPRZci(1)],...
		100*[dblFPRZ-vecFPRZci(2)],'linestyle','none','color',vecColZ,'capsize',20);
	errorbar(3,100*[dblFPRA],...
		100*[dblFPRA-vecFPRAci(1)],...
		100*[dblFPRA-vecFPRAci(2)],'linestyle','none','color',vecColA,'capsize',20);
	xlim([0 4]);
	ylim([0 70]);
	title(sprintf('%s, +%d%%, +%d%%',cellTests{intCase},round((dblFPRZ/dblFPRT - 1)*100),round((dblFPRZ/dblFPRA - 1)*100)));
	ylabel('False positive rate (%)');%,'Color',[0 0.7 0]);
	ax = gca;
	set(ax,'xtick',[1 2 3]);
	ax.XTickLabel{1} = '\color{black}T-test';
	ax.XTickLabel{2} = '\color[rgb]{0,0.4470,0.7410}ZETA';
	ax.XTickLabel{3} = '\color[rgb]{0.8,0,0}ANOVA';
	%%}
	
	
	axes(hAx0);
	errorbar(dblFPRT,dblTPRT,... %errorbar(x,y,yneg,ypos,xneg,xpos)
		[dblTPRT-vecTPRTci(1)],...
		[dblTPRT-vecTPRTci(2)],...
		[dblFPRT-vecFPRTci(1)],...
		[dblFPRT-vecFPRTci(2)],...
		'linestyle','none','color',vecColT,'capsize',5);
	errorbar(dblFPRZ,dblTPRZ,... %errorbar(x,y,yneg,ypos,xneg,xpos)
		[dblTPRZ-vecTPRZci(1)],...
		[dblTPRZ-vecTPRZci(2)],...
		[dblFPRZ-vecFPRZci(1)],...
		[dblFPRZ-vecFPRZci(2)],...
		'linestyle','none','color',vecColZ,'capsize',5);
	errorbar(dblFPRA,dblTPRA,... %errorbar(x,y,yneg,ypos,xneg,xpos)
		[dblTPRA-vecTPRAci(1)],...
		[dblTPRA-vecTPRAci(2)],...
		[dblFPRA-vecFPRAci(1)],...
		[dblFPRA-vecFPRAci(2)],...
		'linestyle','none','color',vecColA,'capsize',5);
	
end
xlabel('False positive fraction');
ylabel('Cell inclusion fraction');
title(sprintf('Alpha=%.3f',dblAlpha));
legend({'Norm','T-test','ZETA','ANOVA'},'location','best');
ylim([0 1]);
xlim([0 1]);
fixfig([],[],2,16);

saveas(gcf,fullpath(strFigPath,'SummaryTPFP.tif'));
saveas(gcf,fullpath(strFigPath,'SummaryTPFP.pdf'));


%%
figure(hFig);
fixfig([],[],2,16);


saveas(gcf,fullpath(strFigPath,'SummaryTtestZeta.tif'));
saveas(gcf,fullpath(strFigPath,'SummaryTtestZeta.pdf'));
