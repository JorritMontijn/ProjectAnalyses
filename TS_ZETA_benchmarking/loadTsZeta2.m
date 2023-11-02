function [matMeanP,matZetawrP] = loadTsZeta2(strIndicator)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%% prep
	if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
		strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	else
		strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	end
	strDataPath = fullfile(strPath,'\Data\');
	strFigPath = fullfile(strPath,'\Figs\');
	
	sDirAll1=dir([strDataPath 'TsZeta2' strIndicator '*Q0*sesResamp1001.mat']);
	sDirAll2=dir([strDataPath 'TsZeta2' strIndicator '*Q0*ses-RandResamp1001.mat']);
	
	sDirAll = cat(1,sDirAll1,sDirAll2);
	vecRand = contains({sDirAll.name},'Rand');
	sDirReal = sDirAll(~vecRand);
	sDirRand = sDirAll(vecRand);
	cellMeanP = [];
	cellZetaP = [];
	cellAnovaP = [];
	cellZetaP_wr = [];
	
	matZetaP = [];
	matMeanP = [];
	matZetawrP = [];
	
	for intRandType=1:2
		if intRandType == 1
			sDir = sDirReal;
		else
			sDir = sDirRand;
		end
		intFiles=numel(sDir);
		
		
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			sLoad=load([strDataPath strFile]);
			cellMeanP{intRandType}{intFile} = sLoad.vecTtestP;
			cellZetaP{intRandType}{intFile} = sLoad.vecTsZetaP;
			cellAnovaP{intRandType}{intFile} = sLoad.vecAnovaP;
			cellZetaP_wr{intRandType}{intFile} = sLoad.vecTsZetaP_wr;
		end
	end
	
	%% check if recording is above chance
	vecResp = cellfun(@(x) sum(x<0.05),cellMeanP{1});
	vecTot = cellfun(@numel,cellMeanP{1});
	vecRespR = vecResp./vecTot;
	pBino=bonf_holm(myBinomTest(vecResp,vecTot,0.05,'two'));
	indRemRecs = pBino>0.05;
	cellMeanP{1}(indRemRecs) = [];
	cellMeanP{2}(indRemRecs) = [];
	cellZetaP{1}(indRemRecs) = [];
	cellZetaP{2}(indRemRecs) = [];
	cellAnovaP{1}(indRemRecs) = [];
	cellAnovaP{2}(indRemRecs) = [];
	cellZetaP_wr{1}(indRemRecs) = [];
	cellZetaP_wr{2}(indRemRecs) = [];
	
	vecRandMeanP = cell2vec(cellMeanP{2});
	vecRealMeanP = cell2vec(cellMeanP{1});
	vecRandZetaP = cell2vec(cellZetaP{2});
	vecRealZetaP = cell2vec(cellZetaP{1});
	vecRandAnovaP = cell2vec(cellAnovaP{2});
	vecRealAnovaP = cell2vec(cellAnovaP{1});
	vecRandZetawrP = cell2vec(cellZetaP_wr{2});
	vecRealZetawrP = cell2vec(cellZetaP_wr{1});
	
	matMeanP = cat(2,vecRealMeanP,vecRandMeanP)';
	matZetaP = cat(2,vecRealZetaP,vecRandZetaP)';
	matAnovaP = cat(2,vecRealAnovaP,vecRandAnovaP)';
	matZetawrP = cat(2,vecRealZetawrP,vecRandZetawrP)';
	
	%remove nans
	indRem = any(isnan(matZetaP),1) | any(isnan(matMeanP),1);
	matMeanP(:,indRem)=[];
	matZetaP(:,indRem)=[];
	matAnovaP(:,indRem)=[];
	matZetawrP(:,indRem)=[];
end

