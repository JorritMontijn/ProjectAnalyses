function [vecR2_faces,vecR2_houses,dblR2_cutoff] = getfhpredr2(strSubject,strDataPath,boolRandomize)
	
	% for both discrete and continuous classification
	if ~exist('boolRandomize','var') || isempty(boolRandomize)
		boolRandomize = false;
	end
	
	%subject = 'ca';
	cls='erp';
	sLoad=load(fullpath([strDataPath filesep strSubject],[strSubject '_' cls '_cross_folds']),'*fold*');
	
	clsparms.preselect_r2=.05;
	clsparms.disc_type='linear';
	
	% for continuous classification
	clsparms.minprob=.51; % minimum post probability to be included as event
	clsparms.minTdist=320; %minimum distance in time, units samplesize
	clsparms.sm_pp='y'; % smooth posterior probability prior to finding peaks
	
	cf=1;%:3
	
	traindata_f = sLoad.f_template_train_fold{cf};
	traindata_h = sLoad.h_template_train_fold{cf};
	trainlabels = sLoad.train_events_fold{cf}(:,2);
	testdata_f = sLoad.f_template_test_fold{cf}(sLoad.test_events_fold{cf}(:,1),:);
	testdata_h = sLoad.h_template_test_fold{cf}(sLoad.test_events_fold{cf}(:,1),:);
	testlabels = sLoad.test_events_fold{cf}(:,2);
	
	if boolRandomize
		trainlabels = trainlabels(randperm(numel(trainlabels)));
		testlabels = testlabels(randperm(numel(testlabels)));
	end
	
	%% parameters to use
	preselect_r2=clsparms.preselect_r2; %r2 threshold to keep channels
	disc_type=clsparms.disc_type; %classifier discriminant analysis type
	
	
	%% down-select to discriminable channels only for feature space - decided not to explicitly include face vs house
	for k=1:size(traindata_f,2)
		rf0(k)=rsa(traindata_f(find(trainlabels==2),k),traindata_f(find(trainlabels==0),k)); % 2 is face
		rh0(k)=rsa(traindata_h(find(trainlabels==1),k),traindata_h(find(trainlabels==0),k)); % 1 is house
	end
	f2u=find(abs(rf0)>preselect_r2);
	h2u=find(abs(rh0)>preselect_r2);
	
	%plot
	%figure, plot(rf0,'b*')
	%hold on, plot(rh0,'rs')
	%hold on, plot([0 size(traindata_f,2)],preselect_r2*[1 1],'k-')
	%xlabel('electrode number (concat if both BB and ERP)'), ylabel('signed r2')
	
	%legend('face feat','house feat', 'cutoff','Location','NorthEastOutside')
	%title([subject ' r2 values by electrode, with cutoff, class=' cls ', fold #' num2str(cf)])
	
	vecR2_faces = rf0;
	vecR2_houses = rh0;
	dblR2_cutoff = preselect_r2;