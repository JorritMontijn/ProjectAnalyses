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
cellTests = {'Spike times','Time-series','Two-sample spike t','Two-sample t-s'};
cellTests = {'Neuropixels','GCaMP6f','Two-sample Npx','Two-sample GCaMP6f'};
figure;maxfig;
for intCase=1:4
	if intCase==1
		%zeta
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visualRunDriftingGratingsResamp250.mat');
		sLoad = load(strFile);
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visual-RandRunDriftingGratingsResamp250.mat');
		sLoadRand = load(strFile);
		matMeanP = cat(1,sLoad.vecTtestP,sLoadRand.vecTtestP);
		matZetaP = cat(1,sLoad.vecZetaP_UniStitch,sLoadRand.vecZetaP_UniStitch); %with replacement
	elseif intCase==2
		%ts-zeta
		[matMeanP,matZetaP] = loadTsZeta(strIndicator);
	elseif intCase==3
		%2-sample zeta
		strFile = fullpath(strDataPath, 'Zeta2DataStimDiffV1RunDriftingGratingsResamp500Q0.mat');
		sLoad = load(strFile);
		matMeanP = sLoad.matTtest2';
		matZetaP = sLoad.matZeta2_old'; %with replacement
	elseif intCase==4
		%2-sample ts-zeta
		[matMeanP,matZetaP] = loadTsZeta2(strIndicator);
	end
	
	intN = size(matZetaP,2);
	[dblTPRZ,vecTPRZci] = binofit(sum(matZetaP(1,:)<0.05),intN);
	[dblTPRT,vecTPRTci] = binofit(sum(matMeanP(1,:)<0.05),intN);
	[dblFPRZ,vecFPRZci] = binofit(sum(matZetaP(2,:)<0.05),intN);
	[dblFPRT,vecFPRTci] = binofit(sum(matMeanP(2,:)<0.05),intN);
	
	%% plot
	vecColZ = lines(1);
	vecColT = [0 0 0];
	subplot(2,4,intCase);
	hold on
	bar(1,100*dblTPRT,0.5,'Facecolor',vecColT,'Edgecolor','none');
	bar(2,100*dblTPRZ,0.5,'Facecolor',vecColZ,'Edgecolor','none');
	errorbar(1,100*dblTPRT,...
		100*[dblTPRT-vecTPRTci(1)],...
		100*[dblTPRT-vecTPRTci(2)],'linestyle','none','color',vecColT,'capsize',20);
	errorbar(2,100*[dblTPRZ],...
		100*[dblTPRZ-vecTPRZci(1)],...
		100*[dblTPRZ-vecTPRZci(2)],'linestyle','none','color',vecColZ,'capsize',20);
	xlim([0 3]);
	ylim([0 70]);
	title(sprintf('%s, +%d%%',cellTests{intCase},round((dblTPRZ/dblTPRT - 1)*100)));
	ylabel('Cell inclusion rate (%)');%,'Color',[0 0.7 0]);
	ax = gca;
	set(ax,'xtick',[1 2]);
	ax.XTickLabel{1} = '\color{black}T-test';
	ax.XTickLabel{2} = '\color[rgb]{0,0.4470,0.7410}ZETA';
	
	%{
	subplot(2,4,4+intCase);
	errorbar([1 2],100*[dblFPRT dblFPRZ],...
		100*[dblFPRT-vecFPRTci(1) dblFPRZ-vecFPRZci(1)],...
		100*[dblFPRT-vecFPRTci(2) dblFPRZ-vecFPRZci(2)],'color',[0.8 0 0]);
	xlim([0 3]);
	ylim([0 70]);
	ylabel('False positive rate (%)','Color',[0.7 0 0]);
	ax = gca;
	set(ax,'xtick',[1 2]);
	ax.XTickLabel{1} = '\color{black}T-test';
	ax.XTickLabel{2} = '\color[rgb]{0,0.4470,0.7410}ZETA';
	%}
end
fixfig([],[],2,20)