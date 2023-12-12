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
vecCutOffs = 2.^(-13:0);
intReps = 10000;
%prep vars
matAllClustZ_diff = [];
matAllSubject = [];
matAllCoords = [];

%%

try
	load(fullpath(strDataPath,'ClustEEG.mat'),'sClustEEG');
	
	cellFields = fieldnames(sClustEEG);
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		eval([strField '=sClustEEG.' strField ';']);
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
	%pre-run all subjects
	for intSubject=1:numel(cellSubjects)
		strSubject = cellSubjects{intSubject};
		
		cellRand = {'normal; ','rand; '};
		for intRand=1:2
			%% rand?
			boolRand =intRand==2;
			%% get data
			sData = load(fullpath([strMasterPath filesep 'data' filesep strSubject], [strSubject '_faceshouses'])); %note that bad channels have already been rejected
			sLocs = load(fullpath([strMasterPath filesep 'locs'], [strSubject '_xslocs']));
			%discard non-included areas
			matCoords = sLocs.locs;
			matData = sData.data(:,sLocs.elcode~=20);
			cellElectrodeLocs = area_lbls(sLocs.elcode)';
			cellElectrodeLocs(sLocs.elcode==20)=[];
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
				
				%save data
				if intElectrode==1
					matClustZ_diff_temp = nan(size(matData,2),numel(vecCutOffs));
				end
				
				%get matrices
				%% build reference time and matrices
				%time
				vecRefT1 = getTsRefT(vecTime,vecHouseOn,dblUseMaxDur);
				vecRefT2 = getTsRefT(vecTime,vecFaceOn,dblUseMaxDur);
				
				%set tol
				dblSuperResFactor = 1;
				dblSampInterval = (median(diff(vecRefT1)) + median(diff(vecRefT2)))/2;
				dblTol = dblSampInterval/dblSuperResFactor;
				vecRefC = cat(1,vecRefT1(:),vecRefT2(:));
				vecRefT = uniquetol(vecRefC,(dblTol*0.99)/max(abs(vecRefC)));
				intT = numel(vecRefT);
				
				%matrices
				vecData = matData(:,intElectrode);
				matCond1 = getInterpolatedTimeSeries(vecTime',vecData,vecHouseOn,vecRefT);
				matCond2 = getInterpolatedTimeSeries(vecTime',vecData,vecFaceOn,vecRefT);
				
				%cluster analysis
				vecT = [];
				for intCutOffIdx=1:numel(vecCutOffs)
					dblCutOff = vecCutOffs(intCutOffIdx);
					[dblClustP,sClustPos,sClustNeg] = clustertest(matCond1,matCond2,intReps,[],dblCutOff);
					matClustZ_diff_temp(intElectrode,intCutOffIdx) = -norminv(dblClustP/2);
				end
				
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
				matClustZ_diff = matClustZ_diff_temp;
			else
				matClustZ_diff(:,:,2) = matClustZ_diff_temp;
			end
		end
		matAllCoords = cat(1,matAllCoords,matCoords);
		matAllClustZ_diff = cat(1,matAllClustZ_diff,matClustZ_diff);
		matAllSubject = cat(1,matAllSubject,intSubject*ones(size(matClustZ_diff)));
	end
	% save
	sClustEEG=struct;
	sClustEEG.matAllCoords = matAllCoords;
	sClustEEG.matAllClustZ_diff = matAllClustZ_diff;
	sClustEEG.matAllSubject = matAllSubject;
	sClustEEG.vecCutOffs = vecCutOffs;
	save(fullpath(strDataPath,'ClustEEG'),'sClustEEG');
end

%% plot
%load zeta
sLoad=load(fullpath(strDataPath,'ZetaEEG.mat'),'sZetaEEG');
matAllZetaZ_diff = sLoad.sZetaEEG.matAllZetaZ_diff;
intN = size(matAllZetaZ_diff,1);
matAllZetaP_diff = 2-2*normcdf(matAllZetaZ_diff);
matAllClustZ_diff = sLoad.sZetaEEG.matAllClustZ_diff;
[dblTPR_Z,vecTPR_Z] = binofit( sum(matAllZetaP_diff(:,1)<0.05),intN);
[dblFPR_Z,vecFPR_Z] = binofit( sum(matAllZetaP_diff(:,2)<0.05),intN);
%remove nans
vecClusterAlphas = sClustEEG.vecCutOffs;
intAlphaNum = numel(vecClusterAlphas);
matClustTPRs = nan(3,intAlphaNum);
matClustFPRs = nan(3,intAlphaNum);
matColor = parula(intAlphaNum);
figure;maxfig;
hAx1=subplot(2,3,1);hold on;
hAx4=subplot(2,3,4);hold on;
for intAlphaIdx=1:intAlphaNum
	dblAlpha = vecClusterAlphas(intAlphaIdx);
	matZ = squeeze(sClustEEG.matAllClustZ_diff(:,intAlphaIdx,:));
	matP = 2-2*normcdf(matZ);
	[dblTPR,vecTPRci] = binofit(sum(matP(:,1)<0.05),intN);
	[dblFPR,vecFPRci] = binofit(sum(matP(:,2)<0.05),intN);
	dblSensitivity = dblTPR/dblFPR;
	matClustTPRs(:,intAlphaIdx) = cat(1,dblTPR,vecTPRci');
	matClustFPRs(:,intAlphaIdx) = cat(1,dblFPR,vecFPRci');
	
	%% plot
	% AUCs
	axes(hAx1);
	vecCol=matColor(intAlphaIdx,:);
	
	intCells = size(matP,1);
	vecBothData = cat(1,matP(:,1),matP(:,2));
	vecBothLabels = cat(1,zeros(size(matP(:,1))),ones(size(matP(:,1))));
	vecThresholds = sort(vecBothData);
	vecThresholds(isnan(vecThresholds))=1;
	vecRealP = matP(:,1);
	vecShuffP = matP(:,2);
	%remove if both==1
	indRem = vecRealP==1 & vecShuffP==1;
	vecRealP(indRem)=[];
	vecShuffP(indRem)=[];
	
	vecTP = sum(vecRealP<=vecThresholds',1)/sum(~isnan(vecRealP));
	vecFP = sum(vecShuffP<=vecThresholds',1)/sum(~isnan(vecShuffP));
	dblAlpha = 0.05;
	[dblTPR, dblTPR_se] = binofit(sum(vecRealP<dblAlpha),sum(~isnan(vecRealP)),0.05);
	[dblFPR, dblFPR_se] = binofit(sum(vecShuffP<dblAlpha),sum(~isnan(vecShuffP)),0.05);
	
	plot(vecFP,vecTP,'Color',vecCol);
	[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
	vecAUC(intAlphaIdx) = dblAUC;
	vecAUC_se(intAlphaIdx) = Ase;
	
	%fprs
	axes(hAx4);
	
	vecRandSorted = sort(matP(:,2));
	vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
	%plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
	vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
	vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
	plot(vecAlphas,vecFPR,'Color',vecCol);
	
	
end

%finish
axes(hAx1);
matZ_P = 2-2*normcdf(matAllZetaZ_diff);
matC_P = 2-2*normcdf(matAllClustZ_diff);
%zeta orig
vecBothData = cat(1,matZ_P(:,1),matZ_P(:,2));
vecThresholds = sort(vecBothData);
vecThresholds(isnan(vecThresholds))=1;
vecRealP = matZ_P(:,1);
vecShuffP = matZ_P(:,2);
%remove if both==1
indRem = vecRealP==1 & vecShuffP==1;
vecRealP(indRem)=[];
vecShuffP(indRem)=[];

vecTP = sum(vecRealP<=vecThresholds',1)/sum(~isnan(vecRealP));
vecFP = sum(vecShuffP<=vecThresholds',1)/sum(~isnan(vecShuffP));

plot(vecFP,vecTP,'Color',lines(1));
%clust orig
vecBothData = cat(1,matC_P(:,1),matC_P(:,2));
vecThresholds = sort(vecBothData);
vecThresholds(isnan(vecThresholds))=1;
vecRealP = matC_P(:,1);
vecShuffP = matC_P(:,2);
%remove if both==1
indRem = vecRealP==1 & vecShuffP==1;
vecRealP(indRem)=[];
vecShuffP(indRem)=[];

vecTP = sum(vecRealP<=vecThresholds',1)/sum(~isnan(vecRealP));
vecFP = sum(vecShuffP<=vecThresholds',1)/sum(~isnan(vecShuffP));

plot(vecFP,vecTP,'Color','r');
%finish
cellLegend = cat(2,cellfun(@num2str,vec2cell(vecClusterAlphas),'uniformoutput',false),{'ZETA'},{'Clust'});
xlabel('Fraction of false positives');
ylabel('Fraction of included cells');
legend(cellLegend);

%fprs
axes(hAx4);
%zeta
vecRandSorted = sort(matZ_P(:,2));
vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
plot(vecAlphas,vecFPR,'Color',lines(1));
%clust
vecRandSorted = sort(matC_P(:,2));
vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
plot(vecAlphas,vecFPR,'Color','r');


%finish
xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('Fraction of false positives'));
set(gca,'xscale','log','yscale','log');
dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
plot([dblMinVal 1],[dblMinVal 1],'k--');
hold off;
xlim([1e-2 1]);
ylim([1e-2 1]);


%plot inclusions/fprs
hAx2=subplot(2,3,2);hold on
hAx5=subplot(2,3,5);hold on
[dblTPR_Z,vecTPR_Z] = binofit( sum(matAllZetaP_diff(:,1)<0.05),intN);
[dblFPR_Z,vecFPR_Z] = binofit( sum(matAllZetaP_diff(:,2)<0.05),intN);
axes(hAx2);
plot(vecClusterAlphas,matClustTPRs(1,:),'Color','r');
plot(vecClusterAlphas,matClustTPRs(2,:),'--','Color','r');
plot(vecClusterAlphas,matClustTPRs(3,:),'--','Color','r');
set(hAx2,'xscale','log')
axes(hAx5);
plot(vecClusterAlphas,matClustFPRs(1,:),'Color','r');
plot(vecClusterAlphas,matClustFPRs(2,:),'--','Color','r');
plot(vecClusterAlphas,matClustFPRs(3,:),'--','Color','r');
set(hAx5,'xscale','log')
fixfig([],[],2,16)