%clear data
clearvars -except strRunType strRunStim sAggABI;

%set paths
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPathABI = '';
	strDataPathSim = 'Z:\Data\Processed\Simulations\';
	strDataPathSimT0 = 'F:\Data\Processed\PopTimeCoding\';
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPathABI = 'E:\AllenBrainVisualEphys';
	strDataPathSim = 'Z:\Data\Processed\Simulations\';
	strDataPathSimT0 = 'F:\Data\Processed\PopTimeCoding\';
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end
