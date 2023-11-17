%%
%https://github.com/meagmohit/EEG-Datasets
%https://purl.stanford.edu/xd109qh3109
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004660

%% 
strDataPath = 'D:\Data\Processed\EEG\fhpred\data';
strCodePath = 'D:\Data\Processed\EEG\fhpred';
sFiles = dir(fullpath(strDataPath,'*.mat'));
cellFiles = {sFiles.name};
cellSubjects = cellfun(@(x) x(1:2),cellFiles,'uniformoutput',false);
strOldPath = cd(strCodePath);
fhpred_master(cellSubjects{1});