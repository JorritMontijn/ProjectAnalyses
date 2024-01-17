
%what to load?
strDataPath = 'D:\Data\Processed\Steinmetz\';
sRecs = dir(strDataPath);
cellSubPaths = {sRecs.name};
indIsDir = cell2vec({sRecs.isdir})' & cellfun(@(x) ~strcmp(x,'.') & ~strcmp(x,'..'),cellSubPaths);
cellSubPaths = cellSubPaths(indIsDir);
for intDataset = 1:numel(cellSubPaths)
	
	%load data
	strRec = cellSubPaths{intDataset};
	strSesPath = fullpath(strDataPath,strRec);
	fprintf('\nTransforming rec %d/%d (%s) [%s]\n',intDataset,numel(cellSubPaths),strRec,getTime);
	sAP = loadSteinmetzSession(strSesPath);
	
	% save
	strTargetFile = [strRec '.mat'];
	save(fullpath(strDataPath,strTargetFile),'sAP');
	
	fprintf('Transformed rec %d/%d (%s) [%s]\n',intDataset,numel(cellSubPaths),strRec,getTime);
end