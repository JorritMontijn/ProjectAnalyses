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
	matAnovaP = sLoad.matAnova2_optimal;
	matMeanP = sLoad.matTtest2;
	matZetaP = sLoad.matZeta2; %with replacement
	%remove empties
	indRem = any(matZetaP==1,2) | any(matMeanP==1,2) | any(matAnovaP==1,2) ...
		|  any(isnan(matZetaP),2) | any(isnan(matZetaP),2) | any(isnan(matZetaP),2);
	matAnovaP(indRem,:) = [];
	matMeanP(indRem,:) = [];
	matZetaP(indRem,:) = [];
	
	%flatten
	matMeanZ = -norminv(matMeanP/2);
	matZetaZ = -norminv(matZetaP/2); %with replacement
	matAnovaZ = -norminv(matAnovaP/2);
	
	%plot ROC
	dblCapZ = 40;
	dblCapP = normcdf(dblCapZ,'upper')*2;
	
	
	matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
	matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
	matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
	
	matMeanZ(matMeanZ>dblCapZ) = dblCapZ;
	matZetaZ(matZetaZ>dblCapZ) = dblCapZ;
	matAnovaZ(matAnovaZ>dblCapZ) = dblCapZ;
	
	matMeanP(matMeanP<dblCapP) = dblCapP;
	matZetaP(matZetaP<dblCapP) = dblCapP;
	matAnovaP(matAnovaP<dblCapP) = dblCapP;
end

