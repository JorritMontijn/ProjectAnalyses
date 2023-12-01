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
dblUseMaxDur = 0.6;
			
try
	load(fullpath(strDataPath,'ZetaEEG_0.6s.mat'),'sZetaEEG');
	
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
	matAllZetaZ_faces = [];
	matAllZetaZ_houses = [];
	matAllZetaZ_diff = [];
	matAllTtestZ_faces = [];
	matAllTtestZ_houses = [];
	matAllTtestZ_diff = [];
	matAllClustZ_diff = [];
	matAllR2_faces = [];
	matAllR2_houses = [];
	matAllSubject = [];
	
	%%
	%pre-run all subjects
	for intSubject=1:numel(cellSubjects)
		strSubject = cellSubjects{intSubject};
		strTargetFile = [strMasterPath filesep 'data' filesep strSubject filesep strSubject '_erp_cross_folds.mat'];
		%check if files exist
		if ~exist(strTargetFile,'file')
			fhpred_master;
		end
		cellRand = {'normal; ','rand; '};
		for intRand=1:2
			%% rand?
			boolRand =intRand==2;
			%% get data
			[vecR2_faces,vecR2_houses,dblR2_cutoff] = getfhpredr2(strSubject,strSourceDataPath,boolRand);
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
			
			%% re-reference / regress out 1st mode
			matData=car(matData);
			dblFs = sData.srate;
			vecTime = (1:size(matData,1))/dblFs;
			
			%% get events
			matCol = lines(2);
			%go through electrodes
			
			for intElectrode = 1:size(matData,2)
				fprintf('Running %ssubject %s (%d/%d), electrode %d/%d [%s]\n',...
					cellRand{intRand},strSubject,intSubject,numel(cellSubjects),intElectrode,size(matData,2),getTime);
				
				%% define stims
				intITI = 101;
				intPrePost = 0;
				vecHouses = 1:50;
				vecFaces = 51:100;
				if boolRand
					intHouses = numel(vecHouses);
					intFaces = numel(vecFaces);
					dblHouseRatio = intHouses/(intHouses+intFaces);
					
					intHouse2House = round(intHouses*dblHouseRatio);
					intHouse2Face = intHouses-intHouse2House;
					intFace2Face = round(intFaces*(1-dblHouseRatio));
					intFace2House = intFaces-intFace2Face;
					
					vecRandHouses = vecHouses(randperm(intHouses));
					vecRandFaces = vecFaces(randperm(intFaces));
					
					vecHouses = sort(cat(2,vecRandHouses(1:intHouse2House),vecRandFaces(1:intFace2House)));
					vecFaces = sort(cat(2,vecRandHouses((intHouse2House+1):end),vecRandFaces((intFace2House+1):end)));
				end
				%stim times
				vecHouseOn = find(~ismember(sData.stim(1:(end-1)),vecHouses) & ismember(sData.stim(2:end),vecHouses))/dblFs;
				vecHouseOff = find(ismember(sData.stim(1:(end-1)),vecHouses) & ~ismember(sData.stim(2:end),vecHouses))/dblFs;
				vecFaceOn = find(~ismember(sData.stim(1:(end-1)),vecFaces) & ismember(sData.stim(2:end),vecFaces))/dblFs;
				vecFaceOff = find(ismember(sData.stim(1:(end-1)),vecFaces) & ~ismember(sData.stim(2:end),vecFaces))/dblFs;
				if boolRand
					vecHouseOn = sort(vecHouseOn + (rand(size(vecHouseOn))-0.5)*dblUseMaxDur*4);
					vecHouseOff = sort(vecHouseOff + (rand(size(vecHouseOff))-0.5)*dblUseMaxDur*4);
					vecFaceOn = sort(vecFaceOn + (rand(size(vecFaceOn))-0.5)*dblUseMaxDur*4);
					vecFaceOff = sort(vecFaceOff + (rand(size(vecFaceOff))-0.5)*dblUseMaxDur*4);
				end
				
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
				matCond1 = sZETA2.matDataPerTrial1;
				matCond2 = sZETA2.matDataPerTrial2;
				[dblClustP,sClustPos,sClustNeg] = clustertest(matCond1,matCond2,1000);
				vecClustZ_diff(intElectrode) = -norminv(dblClustP/2);
				
				%% plot
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
					plot(vecT,sClustPos.tmap);
					plot(vecT,sClustNeg.tmap);
					plot(vecT,sClustPos.map);
					plot(vecT,-sClustNeg.map);
					plot(vecT([1 end]),1.963*[1 1]);
					plot(vecT([1 end]),-1.963*[1 1]);
					hold off
					title(sprintf('crit-sum=%.3f; +Sum=%.3f,p=%.3f; -Sum=%.3f,p=%.3f',sClustNeg.cluscrit,sClustPos.clustsum,sClustPos.p,sClustNeg.clustsum,sClustNeg.p));
					xlim([0 dblUseMaxDur]);
					xlabel('Time (s)');
					ylabel('t-tstatistic per 1 ms bin');
					maxfig;fixfig;
					export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('Clust%s_E%02d',strSubject,intElectrode) '.png']));
					export_fig(fullpath(strFigPath,['EEG_examples\ExampleEEG_' sprintf('Clust%s_E%02d',strSubject,intElectrode) '.pdf']));
					
				end
			end
			
			%plot electrode responsiveness in 3d
			if ~boolRand
				matZetaZ_faces = cell2vec({sMultiZF.dblZETA});
				matZetaZ_houses = cell2vec({sMultiZH.dblZETA});
				matZetaZ_diff = cell2vec({sMultiZ2.dblZETA});
				matTtestZ_faces = cell2vec({sMultiZF.dblMeanZ});
				matTtestZ_houses = cell2vec({sMultiZH.dblMeanZ});
				matTtestZ_diff = cell2vec({sMultiZ2.dblMeanZ});
				matClustZ_diff = vecClustZ_diff;
				matR2_faces = vecR2_faces';
				matR2_houses = vecR2_houses';
				
			else
				matZetaZ_faces(:,2) = cell2vec({sMultiZF.dblZETA});
				matZetaZ_houses(:,2) = cell2vec({sMultiZH.dblZETA});
				matZetaZ_diff(:,2) = cell2vec({sMultiZ2.dblZETA});
				matTtestZ_faces(:,2) = cell2vec({sMultiZF.dblMeanZ});
				matTtestZ_houses(:,2) = cell2vec({sMultiZH.dblMeanZ});
				matTtestZ_diff(:,2) = cell2vec({sMultiZ2.dblMeanZ});
				matClustZ_diff(:,2) = vecClustZ_diff;
				matR2_faces(:,2) = vecR2_faces';
				matR2_houses(:,2) = vecR2_houses';
			end
		end
		matAllCoords = cat(1,matAllCoords,matCoords);
		matAllCoordsT = cat(1,matAllCoordsT,matCoordsTransformed);
		matAllZetaZ_faces = cat(1,matAllZetaZ_faces,matZetaZ_faces);
		matAllZetaZ_houses = cat(1,matAllZetaZ_houses,matZetaZ_houses);
		matAllZetaZ_diff = cat(1,matAllZetaZ_diff,matZetaZ_diff);
		matAllTtestZ_faces = cat(1,matAllTtestZ_faces,matTtestZ_faces);
		matAllTtestZ_houses = cat(1,matAllTtestZ_houses,matTtestZ_houses);
		matAllTtestZ_diff = cat(1,matAllTtestZ_diff,matTtestZ_diff);
		matAllClustZ_diff = cat(1,matAllClustZ_diff,matClustZ_diff);
		matAllR2_faces = cat(1,matAllR2_faces,matR2_faces);
		matAllR2_houses = cat(1,matAllR2_houses,matR2_houses);
		matAllSubject = cat(1,matAllSubject,intSubject*ones(size(matZetaZ_diff)));
	end
	% save
	sZetaEEG=struct;
	sZetaEEG.matAllCoords = matAllCoords;
	sZetaEEG.matAllCoordsT = matAllCoordsT;
	sZetaEEG.matAllZetaZ_faces =matAllZetaZ_faces;
	sZetaEEG.matAllZetaZ_houses = matAllZetaZ_houses;
	sZetaEEG.matAllZetaZ_diff = matAllZetaZ_diff;
	sZetaEEG.matAllTtestZ_faces = matAllTtestZ_faces;
	sZetaEEG.matAllTtestZ_houses = matAllTtestZ_houses;
	sZetaEEG.matAllTtestZ_diff = matAllTtestZ_diff;
	sZetaEEG.matAllClustZ_diff = matAllClustZ_diff;
	sZetaEEG.matAllR2_faces = matAllR2_faces;
	sZetaEEG.matAllR2_houses = matAllR2_houses;
	sZetaEEG.matAllSubject = matAllSubject;
	save(fullpath(strDataPath,'ZetaEEG'),'sZetaEEG');
end

%% plot
%prep summary fig
hFigAll = figure;maxfig;
hAx1=subplot(2,3,1);hold on;
dblSize = 20;
scatter(hAx1,matAllR2_houses(:,1),matAllZetaZ_houses(:,1),dblSize,lines(1),'o','filled','markerfacealpha',0.8);
scatter(hAx1,matAllR2_houses(:,1),matAllTtestZ_houses(:,1),dblSize,[0.2 0.2 0.2],'o','filled','markerfacealpha',0.8);
title(sprintf('Houses-signif, Z: %d%%/%d%%, T: %d%%/%d%%, R: %d%%/%d%%',...
	round(sum(matAllZetaZ_houses(:,1)>1.96)/numel(matAllZetaZ_houses(:,1))*100),...
	round(sum(matAllZetaZ_houses(:,2)>1.96)/numel(matAllZetaZ_houses(:,2))*100),...
	round(sum(matAllTtestZ_houses(:,1)>1.96)/numel(matAllTtestZ_houses(:,1))*100),...
	round(sum(matAllTtestZ_houses(:,2)>1.96)/numel(matAllTtestZ_houses(:,2))*100),...
	round(sum(matAllR2_houses(:,1)>0.05)/numel(matAllR2_houses(:,1))*100),...
	round(sum(matAllR2_houses(:,2)>0.05)/numel(matAllR2_houses(:,2))*100)...
	))
legend({'ZETA','T-test'},'location','best');
ylabel('ZETA/t-test significance (\sigma)');
xlabel('Signal predictability (R^2)');

hAx2=subplot(2,3,2);hold on;
scatter(hAx2,matAllR2_faces(:,1),matAllZetaZ_faces(:,1),dblSize,lines(1),'o','filled','markerfacealpha',0.8);
scatter(hAx2,matAllR2_faces(:,1),matAllTtestZ_faces(:,1),dblSize,[0.2 0.2 0.2],'o','filled','markerfacealpha',0.8);
title(sprintf('Faces-signif, Z: %d%%/%d%%, T: %d%%/%d%%, R: %d%%/%d%%',...
	round(sum(matAllZetaZ_faces(:,1)>1.96)/numel(matAllZetaZ_faces(:,1))*100),...
	round(sum(matAllZetaZ_faces(:,2)>1.96)/numel(matAllZetaZ_faces(:,2))*100),...
	round(sum(matAllTtestZ_faces(:,1)>1.96)/numel(matAllTtestZ_faces(:,1))*100),...
	round(sum(matAllTtestZ_faces(:,2)>1.96)/numel(matAllTtestZ_faces(:,2))*100),...
	round(sum(matAllR2_faces(:,1)>0.05)/numel(matAllR2_faces(:,1))*100),...
	round(sum(matAllR2_faces(:,2)>0.05)/numel(matAllR2_faces(:,2))*100)...
	))
legend({'ZETA','T-test'},'location','best');
ylabel('ZETA/t-test Significance (\sigma)');
xlabel('Signal predictability (R^2)');

hAx3=subplot(2,3,3);cla(hAx3);hold on;
vecRealC = matAllClustZ_diff(:,1);
vecRealZ = matAllZetaZ_diff(:,1);
vecRealC(isinf(vecRealC))=0;
vecRealC = abs(vecRealC);
[pBino2,z]=bino2test(sum(vecRealZ>1.96),numel(vecRealZ),sum(vecRealC>1.96),numel(vecRealC));
dblZTS = 1.96;
ind00 = vecRealZ<dblZTS & vecRealC<dblZTS;
ind01 = vecRealZ<dblZTS & vecRealC>=dblZTS;
ind10 = vecRealZ>=dblZTS & vecRealC<dblZTS;
ind11 = vecRealZ>=dblZTS & vecRealC>=dblZTS;
scatter(hAx3,vecRealC(ind00),vecRealZ(ind00),dblSize,[0.5 0.5 0.5],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecRealC(ind01),vecRealZ(ind01),dblSize,[1 0 0],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecRealC(ind10),vecRealZ(ind10),dblSize,[0 1 0],'o','filled','markerfacealpha',0.8);
scatter(hAx3,vecRealC(ind11),vecRealZ(ind11),dblSize,[0 0 1],'o','filled','markerfacealpha',0.8);
title(sprintf('Diff-signif, Z: %d%%/%d%%, C: %d%%/%d%%; bino2-p=%.3f',...
	round(sum(vecRealZ>1.96)/numel(vecRealZ)*100),...
	round(sum(matAllZetaZ_diff(:,2)>1.96)/numel(matAllZetaZ_diff(:,2))*100),...
	round(sum(vecRealC>1.96)/numel(vecRealC)*100),...
	round(sum(matAllClustZ_diff(:,2)>1.96)/numel(matAllClustZ_diff(:,2))*100),pBino2...
	))
xlabel('Clustering test significance (\sigma)');
ylabel('ZETA significance (\sigma)');


% finish figs
hAx4=subplot(2,3,4);hold on;
hAx5=subplot(2,3,5);hold on;
hAx6=subplot(2,3,6);hold on;

matCol = cat(2,matAllZetaZ_faces(:,1),zeros(size(matAllZetaZ_houses(:,1))),matAllZetaZ_houses(:,1));
matCol = matCol./max(matCol(:));
vecOpacity = (matAllZetaZ_houses(:,1)+matAllZetaZ_faces(:,1));
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
scatter3(matAllCoords(:,1),matAllCoords(:,2),matAllCoords(:,3),(matAllZetaZ_faces(:,1)/max(matAllZetaZ_faces(:,1)))*vecSize,matAllZetaZ_faces(:,1),'filled');
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
scatter(matAllCoords(:,1),matAllCoords(:,2),(matAllZetaZ_houses(:,1)/max(matAllZetaZ_houses(:,1)))*vecSize,matAllZetaZ_houses(:,1),'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('House responsiveness');
xlim([-100 100]);
ylim([-100 100]);

axes(hAx6);
dblMaxZ = max(max(matAllZetaZ_diff(:,1)),max(matAllClustZ_diff(:,1)));
scatter(matAllCoords(:,1),matAllCoords(:,2),(matAllZetaZ_diff(:,1)/dblMaxZ)*vecSize,matAllZetaZ_diff(:,1),'filled');
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
scatter(matAllCoords(:,1),matAllCoords(:,2),(matAllR2_houses/max(matAllR2_houses))*vecSize,matAllR2_houses,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('R2 houses');
xlim([-100 100]);
ylim([-100 100]);

subplot(2,3,2)
scatter(matAllCoords(:,1),matAllCoords(:,2),(matAllR2_faces/max(matAllR2_faces))*vecSize,matAllR2_faces,'filled');
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

subplot(2,3,4)
vecRespDiff = matAllZetaZ_diff-max(matAllZetaZ_faces,matAllZetaZ_houses);
vecRespDiff(vecRespDiff<0)=0;
vecPlotSize = (0.1+imnorm(vecRespDiff))*vecSize;
scatter(matAllCoords(:,1),matAllCoords(:,2),vecPlotSize,vecRespDiff,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('Z2-max(Zf,Zh)');
xlim([-100 100]);
ylim([-100 100]);

subplot(2,3,5)
vecPlotSize = 0.1+(abs(matAllZetaZ_faces-matAllZetaZ_houses)/max(abs(matAllZetaZ_faces-matAllZetaZ_houses)))*vecSize;
scatter(matAllCoords(:,1),matAllCoords(:,2),vecPlotSize,matAllZetaZ_faces-matAllZetaZ_houses,'filled');
colorbar
xlabel('ML? coords')
ylabel('AP? coords')
zlabel('DV? coords')
title('Zf-Zh');
xlim([-100 100]);
ylim([-100 100]);

matAllClustZ_diff = abs(matAllClustZ_diff);
subplot(2,3,6)
scatter(matAllCoords(:,1),matAllCoords(:,2),0.1+(matAllClustZ_diff/dblMaxZ)*vecSize,matAllClustZ_diff,'filled');
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

%% inclusion summaries
figure;maxfig;
%houses
subplot(2,3,1);hold on;
[pZR,z]=bino2test(sum(matAllZetaZ_houses>1.96),numel(matAllZetaZ_houses),sum(matAllR2_houses>0.05),numel(matAllR2_houses));
[pZT,z]=bino2test(sum(matAllZetaZ_houses>1.96),numel(matAllZetaZ_houses),sum(matAllTtestZ_houses>1.96),numel(matAllTtestZ_houses));
[pTR,z]=bino2test(sum(matAllTtestZ_houses>1.96),numel(matAllTtestZ_houses),sum(matAllR2_houses>0.05),numel(matAllR2_houses));
[dblZ,vecZ_ci]=binofit(sum(matAllZetaZ_houses>1.96),numel(matAllZetaZ_houses));
[dblT,vecT_ci]=binofit(sum(matAllTtestZ_houses>1.96),numel(matAllTtestZ_houses));
[dblR,vecR_ci]=binofit(sum(matAllR2_houses>0.05),numel(matAllR2_houses));

errorbar(1,dblZ,dblZ-vecZ_ci(1),dblZ-vecZ_ci(2),'x','color',lines(1));
errorbar(2,dblT,dblT-vecT_ci(1),dblT-vecT_ci(2),'x','color','k');
errorbar(3,dblR,dblR-vecR_ci(1),dblR-vecR_ci(2),'x','color','r');
hold off;
xlim([0 4]);
set(gca,'xtick',[1 2 3],'xticklabel',{'T-ZETA2','T-test','R^2'});
ylim([0 1]);
ylabel('# of house resp. sites');
title(sprintf('houses, Bino2-p; Z-R=%.1e;Z-T=%.1e;T-R=%.1e',pZR,pZT,pTR));

%faces
subplot(2,3,2);hold on;
[pZR,z]=bino2test(sum(matAllZetaZ_faces>1.96),numel(matAllZetaZ_faces),sum(matAllR2_faces>0.05),numel(matAllR2_faces));
[pZT,z]=bino2test(sum(matAllZetaZ_faces>1.96),numel(matAllZetaZ_faces),sum(matAllTtestZ_faces>1.96),numel(matAllTtestZ_faces));
[pTR,z]=bino2test(sum(matAllTtestZ_faces>1.96),numel(matAllTtestZ_faces),sum(matAllR2_faces>0.05),numel(matAllR2_faces));
[dblZ,vecZ_ci]=binofit(sum(matAllZetaZ_faces>1.96),numel(matAllZetaZ_faces));
[dblT,vecT_ci]=binofit(sum(matAllTtestZ_faces>1.96),numel(matAllTtestZ_faces));
[dblR,vecR_ci]=binofit(sum(matAllR2_faces>0.05),numel(matAllR2_faces));

errorbar(1,dblZ,dblZ-vecZ_ci(1),dblZ-vecZ_ci(2),'x','color',lines(1));
errorbar(2,dblT,dblT-vecT_ci(1),dblT-vecT_ci(2),'x','color','k');
errorbar(3,dblR,dblR-vecR_ci(1),dblR-vecR_ci(2),'x','color','r');
hold off;
xlim([0 4]);
set(gca,'xtick',[1 2 3],'xticklabel',{'T-ZETA2','T-test','R^2'});
ylim([0 1]);
ylabel('# of face resp. sites');
title(sprintf('faces, Bino2; Z-R=%.1e;Z-T=%.1e;T-R=%.1e',pZR,pZT,pTR));

%diff
subplot(2,3,3);hold on;
matAllClustZ_diff(isinf(matAllClustZ_diff))=0;
matAllClustZ_diff = abs(matAllClustZ_diff);
[pZC,z]=bino2test(sum(matAllZetaZ_diff>1.96),numel(matAllZetaZ_diff),sum(matAllClustZ_diff>1.96),numel(matAllClustZ_diff));
[pZT,z]=bino2test(sum(matAllZetaZ_diff>1.96),numel(matAllZetaZ_diff),sum(matAllTtestZ_diff>1.96),numel(matAllTtestZ_diff));
[pTC,z]=bino2test(sum(matAllTtestZ_diff>1.96),numel(matAllTtestZ_diff),sum(matAllClustZ_diff>1.96),numel(matAllClustZ_diff));
[dblZ,vecZ_ci]=binofit(sum(matAllZetaZ_diff>1.96),numel(matAllZetaZ_diff));
[dblT,vecT_ci]=binofit(sum(matAllTtestZ_diff>1.96),numel(matAllTtestZ_diff));
[dblC,vecC_ci]=binofit(sum(matAllClustZ_diff>1.96),numel(matAllClustZ_diff));

errorbar(1,dblZ,dblZ-vecZ_ci(1),dblZ-vecZ_ci(2),'x','color',lines(1));
errorbar(2,dblT,dblT-vecT_ci(1),dblT-vecT_ci(2),'x','color','k');
errorbar(3,dblC,dblC-vecC_ci(1),dblC-vecC_ci(2),'x','color','r');
hold off;
xlim([0 4]);
set(gca,'xtick',[1 2 3],'xticklabel',{'T-ZETA2','T-test','Clustering'});
ylim([0 1]);
ylabel('# of diff. resp. sites');
title(sprintf('Diff, Bino2; Z-C=%.1e;Z-T=%.1e;T-C=%.1e',pZC,pZT,pTC));
fixfig;


export_fig(fullpath(strFigPath,'ZetaEEG3.png'));
export_fig(fullpath(strFigPath,'ZetaEEG3.pdf'));
