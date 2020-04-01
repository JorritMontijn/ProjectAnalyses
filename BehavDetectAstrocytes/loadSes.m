%get block data
if ~exist('cellMultiSes','var')
	if intMouse == 1
		strSes = '20140207';
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intMouse == 2
		strSes = '20140314';
		vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intMouse == 3
		strSes = '20140425';
		vecBlock = [1 1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intMouse == 4
		strSes = '20140507';
		vecBlock = [1 1 1]; %define whether neurons are in same population or not
	elseif intMouse == 5
		strSes = '20140530';
		vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intMouse == 6
		strSes = '20140604';
		vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
	elseif intMouse == 7
		strSes = '20140711';
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intMouse == 8
		strSes = '20140715';
		vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
	end
	if intUseNeuropilSubtraction == 1
		strSes = [strSes 'NPS'];
	elseif intUseNeuropilSubtraction == 2
		strSes = [strSes 'NPSPost'];
	end
	load(['D:\Data\Processed\StimDetectionAgg\dataPreProAggregate' strSes '.mat']);
end

%block data & init
vecBlockTypes = unique(vecBlock);
intNumBlocks = length(vecBlockTypes);
clear sLoad sSesAggregate ses
sParams.strFigDir = ['D:\Data\ResultsAstroAnalysis' filesep strSes filesep];
sParams.boolSavePlots = true;
sParams.boolSaveData = true;
strOldDir = cd(sParams.strFigDir);

%change name for no split
if boolExcludeLocomotor
	strSes = ['NL_' strSes];
end