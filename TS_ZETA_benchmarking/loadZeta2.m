function [matMeanP,matZetaP,matAnovaP] = loadZeta2(strFileSearch)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%% prep
	if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
		strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	else
		strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	end
	strDataPath = fullfile(strPath,'\Data\');
	
	
	sDir = dir(fullpath(strDataPath,strFileSearch));
	strFile = sDir(1).name;
	sLoad = load(fullpath(sDir(1).folder,strFile));
	
	cellNeuron = sLoad.cellNeuron;
	matAnovaP_b = sLoad.matAnova2_unbalanced;%matAnova2_optimal matAnova2_unbalanced
	matAnovaP = sLoad.matAnova2_optimal;
	matMeanP = sLoad.matTtest2;
	matZetaP = sLoad.matZeta2; %with replacement
	%remove empties
	indRem = any(matZetaP==1,2) | any(matMeanP==1,2) | any(matAnovaP==1,2) ...
		|  any(isnan(matZetaP),2) | any(isnan(matZetaP),2) | any(isnan(matZetaP),2);
	matAnovaP_b(indRem,:) = [];
	matAnovaP(indRem,:) = [];
	matMeanP(indRem,:) = [];
	matZetaP(indRem,:) = [];
	
	%flatten
	matMeanZ = -norminv(matMeanP/2);
	matZetaZ = -norminv(matZetaP/2); %with replacement
	matAnovaZ_b = -norminv(matAnovaP_b/2);
	matAnovaZ = -norminv(matAnovaP/2);
	
	%plot ROC
	dblCapZ = 40;
	dblCapP = normcdf(dblCapZ,'upper')*2;
	
	
	matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ_b(:))));
	matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ_b(:))));
	matAnovaZ_b(isinf(matAnovaZ_b(:))) = max(matAnovaZ_b(~isinf(matAnovaZ_b(:))));
	matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ_b(~isinf(matAnovaZ_b(:))));
	
	matMeanZ(matMeanZ>dblCapZ) = dblCapZ;
	matZetaZ(matZetaZ>dblCapZ) = dblCapZ;
	matAnovaZ_b(matAnovaZ_b>dblCapZ) = dblCapZ;
	matAnovaZ(matAnovaZ>dblCapZ) = dblCapZ;
	
	matMeanP(matMeanP<dblCapP) = dblCapP;
	matZetaP(matZetaP<dblCapP) = dblCapP;
	matAnovaP_b(matAnovaP_b<dblCapP) = dblCapP;
	matAnovaP(matAnovaP<dblCapP) = dblCapP;
end

