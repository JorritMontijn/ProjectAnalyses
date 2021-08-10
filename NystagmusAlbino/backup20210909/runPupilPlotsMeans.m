
%get data
strDataSourcePath = 'F:\Data\Results\EyeTracking\';
%strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
sFiles = dir([strDataSourcePath '*.mat']);
cellFiles = {sFiles(:).name}';
close all;

%% pre-allocate
sLoad=load([strDataSourcePath sFiles(1).name]);
vecT = sLoad.sPupilResults.vecT;
intBins = numel(vecT);
vecAllRecIdx = [];
vecAllOrientation = [];
matAllX = zeros(0,intBins);
matAllY = zeros(0,intBins);
matAllR = zeros(0,intBins);
matAllXderiv = zeros(0,intBins);
matAllYderiv = zeros(0,intBins);
matAllRderiv = zeros(0,intBins);
matAllL = zeros(0,intBins);
vecAllInclude = false(0);
vecAllAlbino = [];
boolZscore = true;
%% go through files
for intFile=1:numel(cellFiles)
	%load
	sLoad=load([strDataSourcePath sFiles(intFile).name]);
	fprintf('Processing %s\n',sFiles(intFile).name);
	
	%get data
	matThisX = sLoad.sPupilResults.matX;
	matThisY = sLoad.sPupilResults.matY;
	matThisR = sLoad.sPupilResults.matR;
	matThisXderiv = sLoad.sPupilResults.matXderiv;
	matThisYderiv = sLoad.sPupilResults.matYderiv;
	matThisRderiv = sLoad.sPupilResults.matRderiv;
	matThisL = sLoad.sPupilResults.matL;
	
	%z-score?
	if boolZscore
		matThisX = (matThisX - nanmean(matThisX(:)))/nanstd(matThisX(:));
		matThisY = (matThisY - nanmean(matThisY(:)))/nanstd(matThisY(:));
		matThisR = (matThisR - nanmean(matThisR(:)))/nanstd(matThisR(:));
		matThisXderiv = (matThisXderiv - nanmean(matThisXderiv(:)))/nanstd(matThisXderiv(:));
		matThisYderiv = (matThisYderiv - nanmean(matThisYderiv(:)))/nanstd(matThisYderiv(:));
		matThisRderiv = (matThisRderiv - nanmean(matThisRderiv(:)))/nanstd(matThisRderiv(:));
		matThisL = ((matThisL - nanmean(matThisL(:)))/nanstd(matThisL(:)))*0.2;
	end
	
	%assign
	intTrials = numel(sLoad.sPupilResults.Orientation);
	vecAllRecIdx = cat(1,vecAllRecIdx,ones(intTrials,1)*intFile);
	vecAllOrientation = cat(1,vecAllOrientation,sLoad.sPupilResults.Orientation);
	matAllX = cat(1,matAllX,matThisX);
	matAllY = cat(1,matAllY,matThisY);
	matAllR = cat(1,matAllR,matThisR);
	matAllXderiv = cat(1,matAllXderiv,matThisXderiv);
	matAllYderiv = cat(1,matAllYderiv,matThisYderiv);
	matAllRderiv = cat(1,matAllRderiv,matThisRderiv);
	matAllL = cat(1,matAllL,matThisL);
	vecAllInclude = cat(1,vecAllInclude,repmat(~any(isnan(mean(sLoad.sPupilResults.matL,1))),[intTrials 1]));
	vecAllAlbino = cat(1,vecAllAlbino,repmat(contains(sLoad.sPupilResults.strRec,'MA'),[intTrials 1]));
end

%% remove bad recordings
vecAllRecIdx = vecAllRecIdx(vecAllInclude);
vecAllOrientation = vecAllOrientation(vecAllInclude);
matAllX = matAllX(vecAllInclude,:);
matAllY = matAllY(vecAllInclude,:);
matAllR = matAllR(vecAllInclude,:);
matAllXderiv = matAllXderiv(vecAllInclude,:);
matAllYderiv = matAllYderiv(vecAllInclude,:);
matAllRderiv = matAllRderiv(vecAllInclude,:);
matAllL = matAllL(vecAllInclude,:);
vecAllAlbino = vecAllAlbino(vecAllInclude);
vecLimX = [vecT(1) vecT(end)];

%% plot grand mean
for boolDeriv = [false true]
	for boolAlbino = [false true]
		for boolCentered = [false true]
			%% select data
			indUseData = vecAllAlbino==boolAlbino;
			if boolAlbino
				strAlbino = 'DBA';
			else
				strAlbino = 'BL6';
			end
			
			if boolDeriv
				strDeriv = 'Speed';
				matX = matAllXderiv(indUseData,:);
				matY = matAllYderiv(indUseData,:);
				matR = matAllRderiv(indUseData,:);
				matL = matAllL(indUseData,:);
			else
				strDeriv = 'Loc';
				matX = matAllX(indUseData,:);
				matY = matAllY(indUseData,:);
				matR = matAllR(indUseData,:);
				matL = matAllL(indUseData,:);
			end
			%plot
			if boolCentered
				strCentered = 'Normal';
			else
				matX = matX - nanmean(matX,1);
				matY = matY - nanmean(matY,1);
				strCentered = 'Mean-subtracted';
			end
			
			h=figure;
			maxfig;
			%hAx=gca;
			%plot(hAx,vecRefT,matSyncLumPerTrial');
			%hAx.ColorOrder = redbluepurple;
			subplot(3,3,5);
			hold on;
			plot(vecT,nanmean(matL,1),'k');
			errorfill(vecT, nanmean(matR,1), nanstd(matR,[],1)/sqrt(sum(indUseData)),'color',lines(1))
			hold off
			title(sprintf('%s %s %s; radius+sync',strCentered,strAlbino,strDeriv));
			fixfig;
			xlim(vecLimX);
			
			intPlotOris = 8;
			vecPlotOrder = [6 3 2 1 4 7 8 9];
			vecHandles = nan(size(vecPlotOrder));
			dblOriStep = (pi/(intPlotOris/2));
			vecOriRad = deg2rad(vecAllOrientation(indUseData));
			for intPlotOri=1:8
				%get which oris to include
				dblCenterOri = (intPlotOri-1)*dblOriStep;
				indUseOriTrials = abs(circ_dist(dblCenterOri,vecOriRad)) < (dblOriStep*0.49);
				
				%plot
				subplot(3,3,vecPlotOrder(intPlotOri));
				hold on;
				errorfill(vecT, nanmean(matX(indUseOriTrials,:),1), nanstd(matX(indUseOriTrials,:),[],1)/sqrt(sum(indUseOriTrials)),[1 0 0]);
				errorfill(vecT, nanmean(matY(indUseOriTrials,:),1), nanstd(matY(indUseOriTrials,:),[],1)/sqrt(sum(indUseOriTrials)),[0 0 1]);
				hold off
				fixfig;
				title(sprintf('X=r,Y=b,Ori=%.1f deg',rad2deg(dblCenterOri)));
				xlim(vecLimX);
				
			end
			
			%save figure
			strFigFile = sprintf('AveragePupilPlots%s%s%s',strCentered,strAlbino,strDeriv);
			export_fig([strDataSourcePath strFigFile '.tif']);
			export_fig([strDataSourcePath strFigFile '.pdf']);
		end
	end
end