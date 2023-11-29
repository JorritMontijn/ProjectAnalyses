%%
%https://github.com/meagmohit/EEG-Datasets
%https://purl.stanford.edu/xd109qh3109
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004660

%% locs
area_lbls={...
	'Temporal pole', ... %1
	'Parahippocampal gyrus',... %2, parahippocampal part of the medial occipito-temporal gyrus
	'Inferior temporal gyrus',... %3
	'Middle temporal gyrus',... %4
	'fusiform gyrus',... %5 Lateral occipito-temporal gyrus,
	'Lingual gyrus',... %6, lingual part of  the medial occipito-temporal gyrus
	'Inferior occipital gyrus',... %7
	'Cuneus',... %8
	'Post-ventral cingulate gyrus',... %9 Posterior-ventral part of the
	'Middle Occipital gyrus',... %10
	'occipital pole',... %11
	'precuneus',... %12
	'Superior occipital gyrus',... %13
	'Post-dorsal cingulate gyrus',... %14 Posterior-dorsal part of the cingulate gyrus
	' ',...%15
	' ',...%16
	' ',...%17
	' ',...%18
	' ',...%19
	'Non-included area',... %20
	};

%%
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

strSourceDataPath = 'D:\Data\Processed\EEG\fhpred\data';
strMasterPath = 'D:\Data\Processed\EEG\fhpred';
sFiles = dir(strSourceDataPath);
cellFiles = {sFiles.name};
cellSubjects = cellFiles(~contains(cellFiles,'.'));
%cellSubjects = cellfun(@(x) x(1:2),cellFiles,'uniformoutput',false);
strOldPath = cd(strMasterPath);
boolMakeAllPlots = false;
if boolMakeAllPlots,intPlotZ = 4;else,intPlotZ=0;end

try
	load(fullpath(strDataPath,'ZetaEEG'),'sZetaEEG');
	
	cellFields = fieldnames(sZetaEEG);
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		eval([strField '=sZetaEEG.' strField ';']);
	end
	
	%load first subject
	strSubject = cellSubjects{1};
	sLocs = load(fullpath([strMasterPath filesep 'locs'], [strSubject '_xslocs']));
	%discard non-included areas
	matCoords = sLocs.locs;
	cellElectrodeLocs = area_lbls(sLocs.elcode)';
	cellElectrodeLocs(sLocs.elcode==20)=[];
	matCoords(sLocs.elcode==20,:)=[];
	
catch
	%prep vars
	matAllCoords = [];
	matAllCoordsT = [];
	vecAllZetaZ_faces = [];
	vecAllZetaZ_houses = [];
	vecAllZetaZ_diff = [];
	vecAllTtestZ_faces = [];
	vecAllTtestZ_houses = [];
	vecAllTtestZ_diff = [];
	vecAllClustZ_diff = [];
	vecAllR2_faces = [];
	vecAllR2_houses = [];
	vecAllSubject = [];
	
	%%
	%pre-run all subjects
	for intSubject=1:numel(cellSubjects)
		strSubject = cellSubjects{intSubject};
		strTargetFile = [strMasterPath filesep 'data' filesep strSubject filesep strSubject '_erp_cross_folds.mat'];
		%check if files exist
		if ~exist(strTargetFile,'file')
			fhpred_master;
		end
		close all;
		
		%% get data
		[vecR2_faces,vecR2_houses,dblR2_cutoff] = getfhpredr2(strSubject,strSourceDataPath);
		sData = load(fullpath([strMasterPath filesep 'data' filesep strSubject], [strSubject '_faceshouses'])); %note that bad channels have already been rejected
		sLocs = load(fullpath([strMasterPath filesep 'locs'], [strSubject '_xslocs']));
		%discard non-included areas
		matCoords = sLocs.locs;
		matData = sData.data(:,sLocs.elcode~=20);
		cellElectrodeLocs = area_lbls(sLocs.elcode)';
		cellElectrodeLocs(sLocs.elcode==20)=[];
		vecR2_faces(sLocs.elcode==20)=[];
		vecR2_houses(sLocs.elcode==20)=[];
		matCoords(sLocs.elcode==20,:)=[];
		
		vecStim = sData.stim;
		
		% transform coords
		strNifti = [strMasterPath filesep 'brains' filesep strSubject filesep strSubject '_mri.nii'];
		info = niftiinfo(strNifti);
		tform = affine3d(info.Transform.T);
		%transform coords
		matCoords4 = matCoords;
		matCoords4(:,4)=1;
		matCoords4T = matCoords4 * tform.T;
		matCoordsTransformed = matCoords4T(:,1:3);
		
		%define stims
		intITI = 101;
		intPrePost = 0;
		vecHouses = 1:50;
		vecFaces = 51:100;
		
		%% re-reference / regress out 1st mode
		matData=car(matData);
		dblFs = sData.srate;
		vecTime = (1:size(matData,1))/dblFs;
		
		%% get events
		vecHouseOn = find(~ismember(sData.stim(1:(end-1)),vecHouses) & ismember(sData.stim(2:end),vecHouses))/dblFs;
		vecHouseOff = find(ismember(sData.stim(1:(end-1)),vecHouses) & ~ismember(sData.stim(2:end),vecHouses))/dblFs;
		vecFaceOn = find(~ismember(sData.stim(1:(end-1)),vecFaces) & ismember(sData.stim(2:end),vecFaces))/dblFs;
		vecFaceOff = find(ismember(sData.stim(1:(end-1)),vecFaces) & ~ismember(sData.stim(2:end),vecFaces))/dblFs;
		dblUseMaxDur = 0.6;
		matCol = lines(2);
		%go through electrodes
		
		for intElectrode = 1:size(matData,2)
			fprintf('Running subject %s (%d/%d), electrode %d/%d [%s]\n',...
				strSubject,intSubject,numel(cellSubjects),intElectrode,size(matData,2),getTime);
			
			%calc zeta
			vecData = matData(:,intElectrode);
			[dblP,sZETA2]=zetatstest2(vecTime,vecData,cat(2,vecHouseOn,vecHouseOff),vecTime,vecData,cat(2,vecFaceOn,vecFaceOff),dblUseMaxDur,250,intPlotZ);
			[dblP_Houses,sZH]=zetatstest(vecTime,vecData,cat(2,vecHouseOn,vecHouseOff),dblUseMaxDur,250,0,0,[],true,1);
			[dblP_Faces,sZF]=zetatstest(vecTime,vecData,cat(2,vecFaceOn,vecFaceOff),dblUseMaxDur,250,0,0,[],true,1);
			
			%save data
			if intElectrode==1
				sMultiZ2 = sZETA2;
				sMultiZH = sZH;
				sMultiZF = sZF;
				vecClustZ_diff = nan(size(matData,2),1);
			else
				sMultiZ2(intElectrode) = sZETA2;
				sMultiZH(intElectrode) = sZH;
				sMultiZF(intElectrode) = sZF;
			end
			
			%cluster analysis
			%intResampNum = 250;
			%dblJitterSize = 2;
			%boolStitch=true;
			%dblSuperResFactor=100;
			%sClustHouses=clustertest(vecTime,vecData,cat(2,vecHouseOn,vecHouseOff),dblUseMaxDur,intResampNum,dblJitterSize,boolStitch,dblSuperResFactor);
			
			%cluster analysis
			matCond1 = sZETA2.matDataPerTrial1';
			matCond2 = sZETA2.matDataPerTrial2';
			clustersPos = ez_clusterstat_time(matCond1,matCond2,1000);
			clustersNeg = ez_clusterstat_time(matCond2,matCond1,1000);
			vecPosP = cell2vec({clustersPos.p});
			vecNegP = cell2vec({clustersNeg.p});
			dblClustP=min(bonf_holm([min(vecPosP) min(vecNegP)]));
			vecClustZ_diff(intElectrode) = -norminv(dblClustP/2);
			
			%plot
			if boolMakeAllPlots
				hF=gcf;
				h1 = hF.Children(2);
				h2 = hF.Children(4);
				title(h1,sprintf('%s, Electrode %d, %s; Houses',strSubject,intElectrode,cellElectrodeLocs{intElectrode}));
				title(h2,sprintf('%s, Electrode %d, %s; Faces',strSubject,intElectrode,cellElectrodeLocs{intElectrode}));
				
				subplot(2,3,2);cla;
				hold all
				intPlotIters = 50;
				for intIter=1:intPlotIters
					plot(sZH.cellRandT{intIter},sZH.cellRandDiff{intIter},'Color',[0.5 0.5 0.5]);
				end
				plot(sZH.vecTimeSR,sZH.vecD,'Color',matCol(1,:));
				hold off
				xlabel('Time after event (s)');
				ylabel('Data deviation');
				title(sprintf('Houses T-ZETA p=%.3f, t-test p=%.3f, R^2=%.3f',sZH.dblZetaP,sZH.dblMeanP,vecR2_houses(intElectrode)));
				
				subplot(2,3,5);cla;
				hold all
				intPlotIters = 50;
				for intIter=1:intPlotIters
					plot(sZF.cellRandT{intIter},sZF.cellRandDiff{intIter},'Color',[0.5 0.5 0.5]);
				end
				plot(sZF.vecTimeSR,sZF.vecD,'Color',matCol(1,:));
				hold off
				xlabel('Time after event (s)');
				ylabel('Data deviation');
				title(sprintf('Faces T-ZETA p=%.3f, t-test p=%.3f, R^2=%.3f',sZF.dblZetaP,sZF.dblMeanP,vecR2_faces(intElectrode)));
				
				subplot(2,3,6);
				intReps1 = numel(vecHouseOn);
				intReps2 = numel(vecFaceOn);
				vecT = sZETA2.vecRefT;
				vecMu1 = mean(sZETA2.matDataPerTrial1);
				vecSd1 = std(sZETA2.matDataPerTrial1)./sqrt(intReps1);
				vecMu2 = mean(sZETA2.matDataPerTrial2);
				vecSd2 = std(sZETA2.matDataPerTrial2)./sqrt(intReps2);
				hold on
				h1=plot(vecT,vecMu1);
				h2=plot(vecT,vecMu2);
				plot(vecT,vecMu1-vecSd1,'--','color',h1.Color);
				plot(vecT,vecMu1+vecSd1,'--','color',h1.Color);
				plot(vecT,vecMu2-vecSd2,'--','color',h2.Color);
				plot(vecT,vecMu2+vecSd2,'--','color',h2.Color);
				hold off
				title(sprintf('%s, Electrode %d, %s',strSubject,intElectrode,cellElectrodeLocs{intElectrode}));
				legend({'House','Face'},'location','best');
				xlabel('Time after event (s)');
				ylabel('Normalized signal (a.u.)');
				xlim([0 dblUseMaxDur]);
				
				%test single-condition responsiveness
				fixfig;
				export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('%s_E%02d',strSubject,intElectrode) '.png']));
				export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('%s_E%02d',strSubject,intElectrode) '.pdf']));
				
				%% clustering figure
				figure;
				subplot(2,3,1)
				matDiff = matCond1 - matCond2;
				vecMuD = mean(matDiff,2);
				vecSeD = std(matDiff,[],2)/sqrt(size(matDiff,2));
				plot(vecT,vecMuD);
				hold on
				plot(vecT,vecMuD-vecSeD,'--','color',lines(1));
				plot(vecT,vecMuD+vecSeD,'--','color',lines(1));
				hold off
				xlabel('Time (s)');
				ylabel('Signal difference');
				title(sprintf('%s E%d',strSubject,intElectrode));
				
				subplot(2,3,2)
				hold on
				plot(vecT,clustersPos.tmap);
				plot(vecT,clustersNeg.tmap);
				plot(vecT,clustersPos.map);
				plot(vecT,-clustersNeg.map);
				plot(vecT([1 end]),1.963*[1 1]);
				plot(vecT([1 end]),-1.963*[1 1]);
				hold off
				title(sprintf('crit-sum=%.3f; +Sum=%.3f,p=%.3f; -Sum=%.3f,p=%.3f',clustersNeg.cluscrit,clustersPos.clustsum,clustersPos.p,clustersNeg.clustsum,clustersNeg.p));
				xlim([0 dblUseMaxDur]);
				xlabel('Time (s)');
				ylabel('t-tstatistic per 1 ms bin');
				maxfig;fixfig;
				export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('Clust%s_E%02d',strSubject,intElectrode) '.png']));
				export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('Clust%s_E%02d',strSubject,intElectrode) '.pdf']));
				
			end
		end
		
		%plot electrode responsiveness in 3d
		vecR2_faces = vecR2_faces';
		vecR2_houses = vecR2_houses';
		vecZetaZ_faces = cell2vec({sMultiZF.dblZETA});
		vecZetaZ_houses = cell2vec({sMultiZH.dblZETA});
		vecZetaZ_diff = cell2vec({sMultiZ2.dblZETA});
		vecTtestZ_faces = cell2vec({sMultiZF.dblMeanZ});
		vecTtestZ_houses = cell2vec({sMultiZH.dblMeanZ});
		vecTtestZ_diff = cell2vec({sMultiZ2.dblMeanZ});
		
		matAllCoords = cat(1,matAllCoords,matCoords);
		matAllCoordsT = cat(1,matAllCoordsT,matCoordsTransformed);
		vecAllZetaZ_faces = cat(1,vecAllZetaZ_faces,vecZetaZ_faces);
		vecAllZetaZ_houses = cat(1,vecAllZetaZ_houses,vecZetaZ_houses);
		vecAllZetaZ_diff = cat(1,vecAllZetaZ_diff,vecZetaZ_diff);
		vecAllTtestZ_faces = cat(1,vecAllTtestZ_faces,vecTtestZ_faces);
		vecAllTtestZ_houses = cat(1,vecAllTtestZ_houses,vecTtestZ_houses);
		vecAllTtestZ_diff = cat(1,vecAllTtestZ_diff,vecTtestZ_diff);
		vecAllClustZ_diff = cat(1,vecAllClustZ_diff,vecClustZ_diff);
		vecAllR2_faces = cat(1,vecAllR2_faces,vecR2_faces);
		vecAllR2_houses = cat(1,vecAllR2_houses,vecR2_houses);
		vecAllSubject = cat(1,vecAllSubject,intSubject*ones(size(vecZetaZ_diff)));
	end
	% save
	sZetaEEG=struct;
	sZetaEEG.matAllCoords = matAllCoords;
	sZetaEEG.matAllCoordsT = matAllCoordsT;
	sZetaEEG.vecAllZetaZ_faces =vecAllZetaZ_faces;
	sZetaEEG.vecAllZetaZ_houses = vecAllZetaZ_houses;
	sZetaEEG.vecAllZetaZ_diff = vecAllZetaZ_diff;
	sZetaEEG.vecAllTtestZ_faces = vecAllTtestZ_faces;
	sZetaEEG.vecAllTtestZ_houses = vecAllTtestZ_houses;
	sZetaEEG.vecAllTtestZ_diff = vecAllTtestZ_diff;
	sZetaEEG.vecAllClustZ_diff = vecAllClustZ_diff;
	sZetaEEG.vecAllR2_faces = vecAllR2_faces;
	sZetaEEG.vecAllR2_houses = vecAllR2_houses;
	sZetaEEG.vecAllSubject = vecAllSubject;
	save(fullpath(strDataPath,'ZetaEEG'),'sZetaEEG');
end

%% plot
%prep summary fig
hFigAll = figure;maxfig;
hAx1=subplot(2,3,1);hold on;
dblSize = 20;
scatter(hAx1,vecAllR2_houses,vecAllZetaZ_houses,dblSize,lines(1),'o','filled','markerfacealpha',0.8);
scatter(hAx1,vecAllR2_houses,vecAllTtestZ_houses,dblSize,[0.2 0.2 0.2],'o','filled','markerfacealpha',0.8);
title(sprintf('Houses-signif, Z: %d%%, T: %d%%, R: %d%%',...
	round(sum(vecAllZetaZ_houses>1.96)/numel(vecAllZetaZ_houses)*100),...
	round(sum(vecAllTtestZ_houses>1.96)/numel(vecAllTtestZ_houses)*100),...
	round(sum(vecAllR2_houses>0.05)/numel(vecAllR2_houses)*100)...
	))
legend({'ZETA','T-test'},'location','best');
ylabel('ZETA/t-test significance (\sigma)');
xlabel('Signal predictability (R^2)');

hAx2=subplot(2,3,2);hold on;
scatter(hAx2,vecAllR2_faces,vecAllZetaZ_faces,dblSize,lines(1),'o','filled','markerfacealpha',0.8);
scatter(hAx2,vecAllR2_faces,vecAllTtestZ_faces,dblSize,[0.2 0.2 0.2],'o','filled','markerfacealpha',0.8);
title(sprintf('Faces-signif, Z: %d%%, T: %d%%, R: %d%%',...
	round(sum(vecAllZetaZ_faces>1.96)/numel(vecAllZetaZ_faces)*100),...
	round(sum(vecAllTtestZ_faces>1.96)/numel(vecAllTtestZ_faces)*100),...
	round(sum(vecAllR2_faces>0.05)/numel(vecAllR2_faces)*100)...
	))
legend({'ZETA','T-test'},'location','best');
ylabel('ZETA/t-test Significance (\sigma)');
xlabel('Signal predictability (R^2)');

vecAllClustZ_diff(isinf(vecAllClustZ_diff))=0;
vecAllClustZ_diff = abs(vecAllClustZ_diff);
hAx3=subplot(2,3,3);cla(hAx3);hold on;
[pBino2,z]=bino2test(sum(vecAllZetaZ_diff>1.96),numel(vecAllZetaZ_diff),sum(vecAllClustZ_diff>1.96),numel(vecAllClustZ_diff));
dblZTS = 1.96;
ind00 = vecAllZetaZ_diff<dblZTS & vecAllClustZ_diff<dblZTS;
ind01 = vecAllZetaZ_diff<dblZTS & vecAllClustZ_diff>=dblZTS;
ind10 = vecAllZetaZ_diff>=dblZTS & vecAllClustZ_diff<dblZTS;
ind11 = vecAllZetaZ_diff>=dblZTS & vecAllClustZ_diff>=dblZTS;
scatter(hAx3,vecAllClustZ_diff(ind00),vecAllZetaZ_diff(ind00),dblSize,[0.5 0.5 0.5],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecAllClustZ_diff(ind01),vecAllZetaZ_diff(ind01),dblSize,[1 0 0],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecAllClustZ_diff(ind10),vecAllZetaZ_diff(ind10),dblSize,[0 1 0],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecAllClustZ_diff(ind11),vecAllZetaZ_diff(ind11),dblSize,[0 0 1],'o','filled','markerfacealpha',0.8);
title(sprintf('Diff-signif, Z: %d%%, C: %d%%; bino2-p=%.3f',...
	round(sum(vecAllZetaZ_diff>1.96)/numel(vecAllZetaZ_diff)*100),...
	round(sum(vecAllClustZ_diff>1.96)/numel(vecAllClustZ_diff)*100),pBino2...
	))
xlabel('Clustering test significance (\sigma)');
ylabel('ZETA significance (\sigma)');

% finish figs
hAx4=subplot(2,3,4);hold on;
hAx5=subplot(2,3,5);hold on;
hAx6=subplot(2,3,6);hold on;

matCol = cat(2,vecAllZetaZ_faces,zeros(size(vecAllZetaZ_houses)),vecAllZetaZ_houses);
matCol = matCol./max(matCol(:));
vecOpacity = (vecAllZetaZ_houses+vecAllZetaZ_faces);
vecOpacity = vecOpacity/max(vecOpacity(:));
vecSize = 60;
%{
axes(hAx3);
scatter3(matAllCoords(:,1),matAllCoords(:,2),-matAllCoords(:,3),vecOpacity*vecSize,matCol,'filled');
colormap(hAx3,matCol);
title('Redness=face resp; blueness=house resp');
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
xlim([-75 75]);
ylim([-100 50]);
%}
axes(hAx4);
scatter3(matAllCoords(:,1),matAllCoords(:,2),matAllCoords(:,3),(vecAllZetaZ_faces/max(vecAllZetaZ_faces))*vecSize,vecAllZetaZ_faces,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('Face responsiveness');
xlim([-100 100]);
ylim([-100 100]);

%show locations
cellUniqueNames = unique(cellElectrodeLocs);
for intArea=1:numel(cellUniqueNames)
	intEl = find(ismember(cellElectrodeLocs,cellUniqueNames{intArea}),1);
	text(matCoords(intEl,1),matCoords(intEl,2),matCoords(intEl,3), cellUniqueNames{intArea});
end

axes(hAx5);
scatter(matAllCoords(:,1),matAllCoords(:,2),(vecAllZetaZ_houses/max(vecAllZetaZ_houses))*vecSize,vecAllZetaZ_houses,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('House responsiveness');
xlim([-100 100]);
ylim([-100 100]);

axes(hAx6);
dblMaxZ = max(max(vecAllZetaZ_diff),max(vecAllClustZ_diff));
scatter(matAllCoords(:,1),matAllCoords(:,2),(vecAllZetaZ_diff/dblMaxZ)*vecSize,vecAllZetaZ_diff,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('Differential responsiveness');
xlim([-100 100]);
ylim([-100 100]);
fixfig;

%%
export_fig(fullpath(strFigPath,'ZetaEEG.png'));
export_fig(fullpath(strFigPath,'ZetaEEG.pdf'));

%%
figure;maxfig
subplot(2,3,1)
scatter(matAllCoords(:,1),matAllCoords(:,2),(vecAllR2_houses/max(vecAllR2_houses))*vecSize,vecAllR2_houses,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('R2 houses');
xlim([-100 100]);
ylim([-100 100]);

subplot(2,3,2)
scatter(matAllCoords(:,1),matAllCoords(:,2),(vecAllR2_faces/max(vecAllR2_faces))*vecSize,vecAllR2_faces,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('R2 faces');
xlim([-100 100]);
ylim([-100 100]);

hax=subplot(2,3,3);
matBrain = niftiread(niftiinfo('D:\Data\Processed\EEG\fhpred\brains\ca\ca_mri.nii'));
matIm = flipud(squeeze(matBrain(:,round(128-mean(matAllCoords(:,3))),:))');
him=imagesc(-128:127,-128:127,matIm);
axis on
hax.Colormap = gray;
xlim([-100 100]);
ylim([-100 100]);

vecAllClustZ_diff = abs(vecAllClustZ_diff);
subplot(2,3,6)
scatter(matAllCoords(:,1),matAllCoords(:,2),0.1+(vecAllClustZ_diff/dblMaxZ)*vecSize,vecAllClustZ_diff,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('Differential responsiveness');
caxis([0 dblMaxZ])
xlim([-100 100]);
ylim([-100 100]);
fixfig;

export_fig(fullpath(strFigPath,'ZetaEEG2.png'));
export_fig(fullpath(strFigPath,'ZetaEEG2.pdf'));
