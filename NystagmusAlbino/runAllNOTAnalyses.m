%% retrieve files
cellFiles = {...
    'runAnalysisGratings.m'             ,...
    'runAnalysisGratingsOriVsDir.m'     ,...
    'runDecodingGratings.m'             ,...
    'runDecodingGratings2.m'            ,...
    'runDecodingGratings4.m'            ,...
    'runSpikeShapeAnalysis.m'           ,...
    'runSpikeShapeAnalysis2.m'};

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
cellFailed'