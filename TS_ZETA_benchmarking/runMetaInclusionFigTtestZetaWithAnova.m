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

%% load
dblAlpha = 0.05;
cellTests = {'Spike times','Time-series','Two-sample spike t','Two-sample t-s'};
cellTests = {'Neuropixels','GCaMP6f','Two-sample Npx','Two-sample GCaMP6f'};
figure;hAx0=axes;hold on;%plot([dblAlpha dblAlpha],[0 1],'--','color',[0.5 0.5 0.5]);
hFig=figure;maxfig;
vecBinoP = nan(1,4);
for intCase=1:8
	if intCase==1
		%zeta
		strFile=fullpath(strDataPath, 'ZetaDataV1RunNaturalMovieResamp500.mat');
		sLoad = load(strFile);
		matMeanP = sLoad.matTtestP';
		matZetaP = sLoad.matZetaP_Stitch';
		matAnovaP = sLoad.matAnovaP_optimal';
	elseif intCase==2
		%ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta('NM');
	elseif intCase==3
		%2-sample zeta
		strFile = 'Zeta2DataV1RunNaturalMovieResamp500.mat';
		[matMeanP,matZetaP,matAnovaP] = loadZeta2(strFile);
		matMeanP = matMeanP';
		matZetaP = matZetaP';
		matAnovaP = matAnovaP';
	elseif intCase==4
		%2-sample ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta2('TsZeta2NM');
	elseif intCase==5
		%zeta
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visualRunDriftingGratingsResamp500.mat');
		%strFile=fullpath(strDataPath, 'ZetaDataDevPrimary visualRunDriftingGratingsResamp250.mat');
		sLoad = load(strFile);
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visual-RandRunDriftingGratingsResamp500.mat');
		%strFile=fullpath(strDataPath, 'ZetaDataDevPrimary visual-RandRunDriftingGratingsResamp250.mat');
		sLoadRand = load(strFile);
		matMeanP = cat(1,sLoad.vecTtestP,sLoadRand.vecTtestP);
		matZetaP = cat(1,sLoad.vecZetaP_UniStitch,sLoadRand.vecZetaP_UniStitch); %with replacement
		matAnovaP = cat(1,sLoad.vecAnovaP_optimal,sLoadRand.vecAnovaP_optimal);
	elseif intCase==6
		%ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta('GCaMP');
	elseif intCase==7
		%2-sample zeta
		%strFile = 'Zeta2DataAnovaV1RunDriftingGratingsResamp500.mat';
		strFile = 'Zeta2DataStimDiffV1RunDriftingGratingsResamp500.mat';
		%strFile = 'Zeta2DataShiftRespV1RunDriftingGratingsResamp500Q0.mat';
		[matMeanP,matZetaP,matAnovaP] = loadZeta2(strFile);
		matMeanP = matMeanP';
		matZetaP = matZetaP';
		matAnovaP = matAnovaP';
	elseif intCase==8
		%2-sample ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta2('TsZeta3_DiffStims');
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
	
	%dblWLOR_Z = getWeightedLogOddsRatio(matZetaP(2,:)); %fractional logarithmic divergence integral
	dblWLOR_Z = getLiberalness(matZetaP(2,:)); %fractional logarithmic divergence integral
	
	vecFPRZci = [0.9 1.1]*dblWLOR_Z;
	dblFPRZ = dblWLOR_Z;
	%dblWLOR_T = getWeightedLogOddsRatio(matMeanP(2,:));
	dblWLOR_T = getLiberalness(matMeanP(2,:));
	vecFPRTci = [0.9 1.1]*dblWLOR_T;
	dblFPRT = dblWLOR_T;
	%dblWLOR_A = getWeightedLogOddsRatio(matAnovaP(2,:));
	dblWLOR_A = getLiberalness(matAnovaP(2,:));
	vecFPRAci = [0.9 1.1]*dblWLOR_A;
	dblFPRA = dblWLOR_A;

	vecBinoP(intCase) = bino2test(sum(matMeanP(1,:)<dblAlpha),intN,sum(matZetaP(1,:)<dblAlpha),intN);

	%% plot
	intPlot = modx(intCase,4);
	figure(hFig)
	vecColZ = lines(1);
	vecColT = [0 0 0];
	vecColA = [0.8 0 0];
	hAx = subplot(2,4,intPlot);
	hold on
	%bar(1+0.5*(intCase>4),100*dblTPRT,0.5,'Facecolor',vecColT,'Edgecolor','none');
	%bar(2+0.5*(intCase>4),100*dblTPRZ,0.5,'Facecolor',vecColZ,'Edgecolor','none');
	%bar(3+0.5*(intCase>4),100*dblTPRA,0.5,'Facecolor',vecColA,'Edgecolor','none');
	vecX = 1:3;
	vecY = [dblTPRT dblTPRZ dblTPRA];
	vecLow = [dblTPRT-vecTPRTci(1) dblTPRZ-vecTPRZci(1) dblTPRA-vecTPRAci(1)];
	vecHigh = [dblTPRT-vecTPRTci(2) dblTPRZ-vecTPRZci(2) dblTPRA-vecTPRAci(2)];
	plot(hAx,vecX, 100*vecY, 'color', [0.5 0.5 0.5]);
	errorbar(1,100*dblTPRT,...
		100*[dblTPRT-vecTPRTci(1)],...
		100*[dblTPRT-vecTPRTci(2)],'linestyle','none','color',vecColT,'capsize',20);
	errorbar(2,100*[dblTPRZ],...
		100*[dblTPRZ-vecTPRZci(1)],...
		100*[dblTPRZ-vecTPRZci(2)],'linestyle','none','color',vecColZ,'capsize',20);
	errorbar(3,100*[dblTPRA],...
		100*[dblTPRA-vecTPRAci(1)],...
		100*[dblTPRA-vecTPRAci(2)],'linestyle','none','color',vecColA,'capsize',20);
	%xlim([0 4]);
	vecLimY = [100 100 100 100];
	ylim([0 vecLimY(intPlot)]);
	title(sprintf('%s, n=%d, +%d%%, +%d%%',cellTests{intPlot},size(matMeanP(1,:),2),round((dblTPRZ/dblTPRT - 1)*100),round((dblTPRZ/dblTPRA - 1)*100)));
	ylabel('Cell inclusion fraction (%)');%,'Color',[0 0.7 0]);
	ax = gca;
	set(ax,'xtick',[1 2 3]);
	ax.XTickLabel{1} = '\color{black}T-test';
	ax.XTickLabel{2} = '\color[rgb]{0,0.4470,0.7410}ZETA';
	ax.XTickLabel{3} = '\color[rgb]{0.8,0,0}ANOVA';
	
	%%{
	subplot(2,4,4+intPlot);
	hold on
	bar(1+0.5*(intCase>4),dblFPRT,0.5,'Facecolor',vecColT,'Edgecolor','none');
	bar(2+0.5*(intCase>4),dblFPRZ,0.5,'Facecolor',vecColZ,'Edgecolor','none');
	bar(3+0.5*(intCase>4),dblFPRA,0.5,'Facecolor',vecColA,'Edgecolor','none');
	xlim([0 4]);
	%ylim([0 70]);
	title(sprintf('%s, +%d%%, +%d%%',cellTests{intPlot},round((dblFPRZ/dblFPRT - 1)*100),round((dblFPRZ/dblFPRA - 1)*100)));
	ylabel('FPR deviation from norm (WLOR)');%,'Color',[0 0.7 0]);
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
		'linestyle','none','color',vecColT,'capsize',5);
	errorbar(dblFPRZ,dblTPRZ,... %errorbar(x,y,yneg,ypos,xneg,xpos)
		[dblTPRZ-vecTPRZci(1)],...
		[dblTPRZ-vecTPRZci(2)],...
		'linestyle','none','color',vecColZ,'capsize',5);
	errorbar(dblFPRA,dblTPRA,... %errorbar(x,y,yneg,ypos,xneg,xpos)
		[dblTPRA-vecTPRAci(1)],...
		[dblTPRA-vecTPRAci(2)],...
		'linestyle','none','color',vecColA,'capsize',5);
end
%set(gca,'xscale','log');
%xlabel('FPR deviation from norm (WLORt)');
xlabel('Inclusion error (excess FPR)');
ylabel('Cell inclusion fraction');
title(sprintf('Alpha=%.3f',dblAlpha));
legend({'T-test','ZETA','ANOVA'},'location','best');
ylim([0 1]);
%xlim([0 1]);
fixfig([],[],2,16);
xlim([0 2.5]);drawnow;
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_0_2.5.jpg'));
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_0_2.5.pdf'));

%%
xlim([2.5 80]);drawnow;
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_2.5_80.jpg'));
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_2.5_80.pdf'));

%%
xlim([0 80]);drawnow;
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_0_80.jpg'));
saveas(gcf,fullpath(strFigPath,'SummaryExcessError_0_80.pdf'));

%%
figure(hFig);
fixfig([],[],2,16);


saveas(gcf,fullpath(strFigPath,'SummaryTtestZeta.jpg'));
saveas(gcf,fullpath(strFigPath,'SummaryTtestZeta.pdf'));
