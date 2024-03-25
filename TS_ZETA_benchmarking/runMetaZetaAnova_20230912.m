clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

cellUniqueAreas = {...
	'HeteroPoissonPeak',...Area 1
	'TriPhasic',...Area 2
	'QuadriPhasic',...Area 3
	'iidGaussian',...Area 4
	'TriHomoPhasic',...Area 5
	'',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8; DG
	'V1'};%area 9; NM

cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};
cellRepStr = {...
	'RunDriftingGratings','';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
	};

%% prep
cellDatasetNames = {};
matNumCells = [];

intIdxNpx = 0;


intResamps = 500;
boolDirectQuantile = false;
strQ = ['Q' num2str(boolDirectQuantile) ];
strR = ['Resamp' num2str(intResamps)];
for intArea=8%[1:4 8]%8%[1:4 8]%[1:4]%1:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if strcmp(strArea,'Primary visual') %DG
		strStim = 'RunDriftingGratings';
		
	elseif strcmp(strArea,'V1') %NM
		strStim = 'RunNaturalMovie';
	else
		strStim = '';
	end
	matSignifZetaOld = [];
	
	matSignifZetaUniStitch = [];
	matSignifZetaUniNoStitch = [];
	matSignifZetaLinStitch = [];
	matSignifZetaLinNoStitch = [];
	
	matSignifTtest = [];
	matSignifAnova = [];
	
	cellTtestP = [];
	cellAnovaP = [];
	cellNumSpikes = [];
	cellUniStitchP = [];
	cellUniNoStitchP = [];
	cellLinStitchP = [];
	cellLinNoStitchP = [];
	
	cellZetaOldP = [];
	intIdx = 1;
	boolSplit = true;
	for intRandType=1:2
		%set var
		strRand = cellRunRand{intRandType};
		
		%% load data
		strRunType = [strArea strRand strStim];
		if ~strcmp(strArea,'V1')
			%strRunType = ['Anova' strRunType];
			%strRunType = ['Dev' strRunType];
			strRunType = ['JitterSpikes' strRunType];
		end
		sDir=dir([strDataPath 'ZetaData' strRunType 'Resamp' num2str(intResamps) '.mat']);
		intFiles=numel(sDir);
		if intFiles == 0 && boolSplit
			error
		end
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			
			sLoad=load([strDataPath strFile]);
			
			%remove cells with too few spikes
			if isfield(sLoad,'vecNumSpikes')
				indRemCells = sLoad.vecNumSpikes < 3;%median(sLoad.vecNumSpikes);
				cellFields = fieldnames(sLoad);
				for intField=1:numel(cellFields)
					varData = sLoad.(cellFields{intField});
					if numel(varData) ~= numel(indRemCells),continue;end
					sLoad.(cellFields{intField}) = varData(~indRemCells);
				end
				
				cellTtestP{intIdx,intRandType} = sLoad.vecTtestP;
				cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
				
				cellUniStitchP{intIdx,intRandType} = sLoad.vecZetaP_UniStitch;
				cellUniNoStitchP{intIdx,intRandType} = sLoad.vecZetaP_UniNoStitch;
				
				if isfield(sLoad,'vecAnovaP_optimal')
					cellAnovaP{intIdx,intRandType} = sLoad.vecAnovaP_optimal;
				else
					cellAnovaP{intIdx,intRandType} = sLoad.vecAnovaP;
				end
			else
				%only one file
				boolSplit = false;
				indRemCells = any(sLoad.matNumSpikes < 3,2);%median(sLoad.vecNumSpikes);
				
				cellFields = fieldnames(sLoad);
				for intField=1:numel(cellFields)
					varData = sLoad.(cellFields{intField});
					if numel(varData) ~= numel(indRemCells),continue;end
					sLoad.(cellFields{intField}) = varData(~indRemCells,:);
				end
				
				cellTtestP{intIdx,1} = sLoad.matTtestP(:,1);
				cellTtestP{intIdx,2} = sLoad.matTtestP(:,2);
				cellNumSpikes{intIdx,1} = sLoad.matNumSpikes(:,1);
				cellNumSpikes{intIdx,2} = sLoad.matNumSpikes(:,2);
				
				cellUniStitchP{intIdx,1} = sLoad.matZetaP_Stitch(:,1);
				cellUniStitchP{intIdx,2} = sLoad.matZetaP_Stitch(:,2);
				cellUniNoStitchP{intIdx,1} = sLoad.matZetaP_NoStitch(:,1);
				cellUniNoStitchP{intIdx,2} = sLoad.matZetaP_NoStitch(:,2);
				
				cellAnovaP{intIdx,1} = sLoad.matAnovaP_optimal(:,1);
				cellAnovaP{intIdx,2} = sLoad.matAnovaP_optimal(:,2);
				
			end
		end
	end
	
	%flatten
	intUseN = min(sum(cellfun(@numel,cellTtestP),1));
	matMeanZ = -norminv(cell2vec(cellTtestP(:,1))/2);
	matMeanZ = matMeanZ(1:intUseN);
	matMeanZ(:,2) = -norminv(cell2vec(cellTtestP(:,2))/2);
	matMeanP = 2-2*normcdf(matMeanZ);
	
	cellUseZ = cellUniStitchP;%cellUniStitchP cellZetaOldP
	matZetaZ = -norminv(cell2vec(cellUseZ(:,1))/2);
	matZetaZ = matZetaZ(1:intUseN);
	matZetaZ(:,2) = -norminv(cell2vec(cellUseZ(:,2))/2);
	matZetaP = 2-2*normcdf(matZetaZ);
	
	cellUse2Z = cellUniNoStitchP;%cellUniStitchP cellZetaOldP
	matZeta2Z = -norminv(cell2vec(cellUse2Z(:,1))/2);
	matZeta2Z = matZeta2Z(1:intUseN);
	matZeta2Z(:,2) = -norminv(cell2vec(cellUse2Z(:,2))/2);
	matZeta2P = 2-2*normcdf(matZeta2Z);
	
	matAnovaZ = -norminv(cell2vec(cellAnovaP(:,1))/2);
	matAnovaZ = matAnovaZ(1:intUseN);
	matAnovaZ(:,2) = -norminv(cell2vec(cellAnovaP(:,2))/2);
	matAnovaP = 2-2*normcdf(matAnovaZ);
	
	%remove nans
	indRem = any(matMeanZ == 0 | matZetaZ == 0 | matZeta2Z == 0 | matAnovaZ == 0,2);
	matMeanP(indRem,:)=[];
	matZetaP(indRem,:)=[];
	matAnovaP(indRem,:)=[];
	matZeta2P(indRem,:)=[];
	matMeanZ(indRem,:)=[];
	matZetaZ(indRem,:)=[];
	matAnovaZ(indRem,:)=[];
	matZeta2Z(indRem,:)=[];
	
	%plot ROC
	matAUCp = [];
	matAUC_dprime = [];
	if size(cellTtestP,1) >= intIdx && ~isempty(cellTtestP{intIdx,1})
		intIdxNpx = intIdxNpx + 1;
		
		figure
		dblAlpha = 0.05;
		matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
		matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
		matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
		matMeanP(matMeanP(:)==0) = 1e-29;
		matZetaP(matZetaP(:)==0) = 1e-29;
		matAnovaP(matAnovaP(:)==0) = 1e-29;
		h1 =subplot(2,3,1);
		matC = [0.5 0.5 0.5;...
			0 0.8 0;...
			0.8 0 0;...
			0 0 0.8];
		vecColor1 = 1 + (matZetaP(:,1) < dblAlpha & matMeanP(:,1) > dblAlpha) + 2*(matZetaP(:,1) > dblAlpha & matMeanP(:,1) < dblAlpha) + 3*(matZetaP(:,1) < dblAlpha & matMeanP(:,1) < dblAlpha);
		scatter(matMeanZ(:,1),matZetaZ(:,1),100,vecColor1,'.');
		colormap(h1,matC(1:max(vecColor1),:));
		%xlim([0 1]);ylim([0 1]);
		xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
		ylabel('ZETA (\zeta_c)')
		strTit = sprintf('A) Inclusion at %s=%.3f: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,1)<dblAlpha)/numel(matZetaP(:,1)),getGreek('mu'),sum(matMeanP(:,1)<dblAlpha)/numel(matMeanP(:,1)),size(matZetaP,1));
		title(strTit)	%set(gca,'xscale','log','yscale','log');
		
		h2=subplot(2,3,2);
		vecColor2 = 1 + 1*(matZetaP(:,2) > dblAlpha & matMeanP(:,2) < dblAlpha) + 2*(matZetaP(:,2) < dblAlpha & matMeanP(:,2) > dblAlpha) + 3*(matZetaP(:,2) < dblAlpha & matMeanP(:,2) < dblAlpha);
		scatter(matMeanZ(:,2),matZetaZ(:,2),100,vecColor1,'.');
		colormap(h2,matC(1:max(vecColor1),:));
		%xlim([0 1]);ylim([0 1]);
		xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
		ylabel('ZETA (\zeta_c)')
		title(sprintf('B) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,2)<dblAlpha)/numel(matZetaP(:,2)),getGreek('mu'),sum(matMeanP(:,2)<dblAlpha)/numel(matMeanP(:,2))))
		%set(gca,'xscale','log','yscale','log');
		
		intNumN = size(matZetaP,1);
		vecFP_sortedZ = sort(matZetaP(:,2));
		vecFP_sortedA = sort(matAnovaP(:,2));
		dblAlphaAtFpAlphaPercZ = vecFP_sortedZ(round(intNumN*dblAlpha));
		dblAlphaAtFpAlphaPercA = vecFP_sortedA(round(intNumN*dblAlpha));
		dblInclusionZ_at_Alpha = sum(matZetaP(:,1)<dblAlphaAtFpAlphaPercZ)/numel(matZetaP(:,1));
		h4 =subplot(2,3,4);
		matC = [0.5 0.5 0.5;...
			0 0.8 0;...
			0.8 0 0;...
			0 0 0.8];
		vecColor1 = 1 + (matZetaP(:,1) < dblAlpha & matAnovaP(:,1) > dblAlpha) + 2*(matZetaP(:,1) > dblAlpha & matAnovaP(:,1) < dblAlpha) + 3*(matZetaP(:,1) < dblAlpha & matAnovaP(:,1) < dblAlpha);
		scatter(matAnovaZ(:,1),matZetaZ(:,1),100,vecColor1,'.');
		colormap(h4,matC(1:max(vecColor1),:));
		%xlim([0 1]);ylim([0 1]);
		xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
		ylabel('ZETA (\zeta_c)')
		title(sprintf('C) Inclusion at FPR=%.3f: %s=%.3f, %s=%.3f; n=%d',dblAlpha,getGreek('zeta'),dblInclusionZ_at_Alpha,'A',sum(matAnovaP(:,1)<dblAlphaAtFpAlphaPercA)/numel(matAnovaP(:,1)),intNumN))
		%set(gca,'xscale','log','yscale','log');
		
		h5=subplot(2,3,5);
		vecColor2 = 1 + 1*(matZetaP(:,2) > dblAlpha & matAnovaP(:,2) < dblAlpha) + 2*(matZetaP(:,2) < dblAlpha & matAnovaP(:,2) > dblAlpha) + 3*(matZetaP(:,2) < dblAlpha & matAnovaP(:,2) < dblAlpha);
		scatter(matAnovaZ(:,2),matZetaZ(:,2),100,vecColor1,'.');
		colormap(h5,matC(1:max(vecColor1),:));
		%xlim([0 1]);ylim([0 1]);
		xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
		ylabel('ZETA (\zeta_c)')
		title(sprintf('D) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,2)<dblAlpha)/numel(matZetaP(:,2)),'A',sum(matAnovaP(:,2)<dblAlpha)/numel(matAnovaP(:,2))))
		%set(gca,'xscale','log','yscale','log');
		
		%% plot ROC
		cellColor = {lines(1),'r','k','b','m'};
		%vecH(intResampNpx) = subplot(4,3,intResampNpx);
		subplot(2,3,3)
		maxfig;
		hold on;
		
		cellLegend = {};
		hold on;
		for intTest=1:4
			if intTest == 1
				matData = matZetaP;
				cellLegend(end+1) = {'ZETA'};
			elseif intTest == 2
				matData = matAnovaP;
				cellLegend(end+1) = {'ANOVA'};
			elseif intTest == 3
				matData = matMeanP;
				cellLegend(end+1) = {'T-test'};
			elseif intTest == 4
				matData = matZeta2P;
				cellLegend(end+1) = {'ZETA-NS'};
			end
			intCells = size(matData,1);
			vecBothData = cat(1,matData(:,1),matData(:,2));
			vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
			vecThresholds = [0; sort(vecBothData); 1];
			vecRealP = matData(:,1);
			vecShuffP = matData(:,2);
			
			vecTP = sum(vecRealP<=vecThresholds',1)/intCells;
			vecFP = sum(vecShuffP<=vecThresholds',1)/intCells;
			
			plot(vecFP,vecTP,'Color',cellColor{intTest});
			
			[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
			vecAUC(intTest) = dblAUC;
			vecAUC_se(intTest) = Ase;
			
			cellLegend{end} = [cellLegend{end} sprintf(', AUC=%.3f',dblAUC)];
		end
		
		%% run tests on aucs
		AUC_T = vecAUC(3) ;
		AUC_A = vecAUC(2);
		AUC_Z = vecAUC(1);
		AUC_Zns = vecAUC(4);
		Ase_T = vecAUC_se(3) ;
		Ase_A = vecAUC_se(2);
		Ase_Z = vecAUC_se(1);
		Ase_Zns = vecAUC_se(4);
		
		%t vs a
		m0 = AUC_T - AUC_A;
		s0 = (Ase_T + Ase_A)/2;
		zTA = abs(m0/s0);
		AUC_pTA = normcdf(zTA,'upper')*2;
		
		%t v z
		m0 = AUC_T - AUC_Z;
		s0 = (Ase_T + Ase_Z)/2;
		zTZ = abs(m0/s0);
		AUC_pTZ = normcdf(zTZ,'upper')*2;
		
		%a vs z
		m0 = AUC_A - AUC_Z;
		s0 = (Ase_A + Ase_Z)/2;
		zAZ = abs(m0/s0);
		AUC_pAZ = normcdf(zAZ,'upper')*2;
		
		%z vs zna
		m0 = AUC_Zns - AUC_Z;
		s0 = (Ase_Zns + Ase_Z)/2;
		zZZns = abs(m0/s0);
		AUC_pZZns = normcdf(zZZns,'upper')*2;
		
		
		%plot
		hold off;
		xlabel('False positive fraction');
		ylabel('Inclusion fraction');
		title(sprintf('E) ROC; %s %s %s',strQ,strR,strArea));
		legend(cellLegend,'location','best');
		%legend({sprintf('ZETA-test, AUC=%.3f',vecAUC(1)),sprintf('ANOVA, AUC=%.3f',vecAUC(2)),sprintf('t-test, AUC=%.3f',vecAUC(3))},'location','best','interpreter','none');
		
		subplot(2,3,6)
		cellLegend = {};
		hold on;
		for intTest=1:4
			if intTest == 1
				matData = matZetaP;
				cellLegend(end+1) = {'ZETA'};
			elseif intTest == 2
				matData = matAnovaP;
				cellLegend(end+1) = {'ANOVA'};
			elseif intTest == 3
				matData = matMeanP;
				cellLegend(end+1) = {'T-test'};
			elseif intTest == 4
				matData = matZeta2P;
				cellLegend(end+1) = {'ZETA-NS'};
			end
			vecRandSorted = sort(matData(:,2));
			vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
			plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
		end
		xlabel(sprintf('Significance level %s',getGreek('alpha')));
		ylabel(sprintf('P-value cut-off needed for FPR=%s',getGreek('alpha')));
		set(gca,'xscale','log','yscale','log');
		dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
		plot([dblMinVal 1],[dblMinVal 1],'k--');
		cellLegend(end+1) = {'Theoretical norm'};
		hold off;
		legend(cellLegend,'location','best');
		title(sprintf('MW AUC tests; Z-Zns,p=%.1e; Z-T,p=%.1e; Z-A,p=%.1e;',...
			AUC_pZZns,AUC_pTZ,AUC_pAZ));
		fixfig;maxfig;
		
		%% save
		drawnow;
		export_fig(fullpath(strFigPath,['Zeta' strQ strR strArea '.tif']));
		export_fig(fullpath(strFigPath,['Zeta' strQ strR strArea '.pdf']));
	end
end

