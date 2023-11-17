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
subject = cellSubjects{1};
%fhpred_master

%% load erp
load(['data/' subject '_erp_cross_folds'],'*fold*')
intFold = 1;

matTrainFace=f_template_train_fold{intFold}; %train data for face templates
matTrainHouse=h_template_train_fold{intFold}; %train data for house templates
matTrainLabel=train_events_fold{intFold}(:,2); %labels of training data
matTestFace=f_template_test_fold{intFold}(test_events_fold{intFold}(:,1),:); %test data for face templates
matTestHouse=h_template_test_fold{intFold}(test_events_fold{intFold}(:,1),:); %test data for house templates
matTestLabel=test_events_fold{intFold}(:,2); %labels of testing data
