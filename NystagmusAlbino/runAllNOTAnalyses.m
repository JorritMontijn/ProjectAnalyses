%% retrieve files
cellRem = {[mfilename '.m'],'runHeaderNOT.m','runCalcRFs.m','runCrossCorrelogram.m'};
[strPath,strFile,strExt] = fileparts(mfilename('fullpath'));

sDir=dir(fullpath(strPath,filesep,'*.m'));
cellFiles = {sDir.name};
cellFiles(ismember(cellFiles,cellRem))=[];
cellFailed = {};

%% run all files
for intFile=1:numel(cellFiles)
	drawnow;
	close all;
	clearvars -except intFile cellFiles cellFailed;
	strFile = cellFiles{intFile};
	[strPath,strName,strExt] = fileparts(strFile);
	try
		%run file
		fprintf('Running %s (%d/%d) [%s]\n',strName,intFile,numel(cellFiles),getTime);
		eval(strName);
	catch ME
		dispErr(ME);
		cellFailed(end+1) = {strFile};
	end
end