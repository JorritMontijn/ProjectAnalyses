%% koen
% Er zitten 2 folders in (van voor trainen en na trainen) met elk 12 files (6 x Res.mat files en 6 x SPSIG.mat files, 1 van beiden voor iedere muis, dus 6 muizen). Relevante variabelen:
% 
% Res.mat file:
% 
%     Res.CaSigCorrected (frames x trials x ROIs). Als je hier de waarde '1' vanaf trekt heb je in principe dF/F. 
%     Res.CaDeconCorrected (frames x trials x ROIs). Deconvolved data, dus estimated spikes per frame.
%     Res.ax. De as tijden, stimuli waren op t = 0 en duurden 1 seconde.
%     info. Allerlei informatie over de image sessie en ROIs.
%         info.frame. Framestamps van de stimuli (240 trials). (handig voor sigCorrected, zie onder).
%         info.FrameTimes. Timestamps van de frames (handig voor sigCorrected, zie onder).
%         info.Stim.Parameters.Time = stim duratie.
%         info.Stim.Parameters.Wait = ITI (maar is semi-random tussen wait en wait+2). 
%         info.Stim.Parameters.StimSmp = log file. (1,:) = image nr, (2,:) = occluded (1) of niet occluded (0)
% 
% SPSIG.mat file:
% 
%     sigCorrected (frames x ROIs). Als je hier de waarde '1' vanaf trekt heb je in principe dF/F,
%     maar nu dus voor de hele image sessie.  
%     deconCorrected (frames x ROIs). Deconvolved data.
%     freq. Image frequentie.
% 
% 
% Ik denk dat je er zo wel uit moet komen. In principe kun je waarschijnlijk af met
% Res.CaSigCorrected, aangezien daar al de trial-aligned data in staat (tijd as = Res.ax). Ieder
% sessie was 240 trials, 20 herhalingen x 6 images x 2 condities (occl vs nonOccl).  
% 
% Ik heb al een aantal analyses gedaan en zie al een paar hele duidelijke dingen, zoals een
% twee-deling van de populatie (non-occluded responders en occluded responders), een sparsificatie
% van de non-occluded responsen, en nog een aantal interessante dingen (o.a. mbt decoding). Ik heb
% nog mega veel analyses op mijn lijstje staan, maar alle ideeën zijn welkom. Als iets onduidelijk
% is hoor ik het wel.     
% 
%% jorrit
%  Wat me ook nog wel interessant leek is om te kijken hoe de responsen tussen full en occluded
%  verschillen op populatie niveau, en of daar generalizable principles uit te halen zijn. Wat he
%  bijvoorbeeld zou kunnen doen is trainen op full en occluded, en de transformatiematrix tussen die
%  twee uitrekenen. Er bestaat een regression matrix B zo dat   
% 1) W_o = B * W_f, dus:
% 2) B = W_f \ W_o
% Vervolgens kun je dus de occluded responses decoden door de weight matrix B * W_f te gebruiken. De
% crux is nu dat je bij eq 1) ook een reduced rank regression kan doen: 
% 3) W_o* = B * R * R' * W_f
% Hierbij kun je de dimensionality van R (een soort low-dimensional funnel) zelf instellen en voor
% elke rank (i.e., dimensionality) een OLS regression doen. Vervolgens krijg je dus voor elke
% dimensionality een andere W_o*, die je kunt gebruiken om je occluded responsen te decoden. Dan kun
% je dus te weten komen wat de complexiteit van de transformatie tussen W_o en W_f is: als je bij
% een rank van 1 al decoding performance van W_o* krijgt die net zo hoog is als W_o, dan zijn W_o en
% W_f dus eigenlijk identitiek na een lineaire schaling. Aan de andere kant kan W_o* pas net zo goed
% worden als W_o wanneer de dimensionality van W_o hetzelfde is als van min(W_o,W_f). In dat geval
% is er totaal geen correlatie in hoe de occluded en full stimuli worden gecodeerd. Het lijkt me
% vooral interessant als je dit ook voor en na training zou kunnen doen, want ik kan me voorstellen
% dat er na training meer structuur in zit (en dus lagere transformation dimensionality heeft).         

%% plan1:
%do pairwise logistic regression to calculate weight vectors that discriminate a pair of stimuli
%q1: are discriminating dimensions lying in low-dimensional subspace?
%q2: does this subspace grow or shrink after training?
%q3: what is cosine similarity of weight vectors between occluded and full conditions?
%q4: does full/occluded cosine similarity become higher/lower after training? => more/less generalization

%% params
%clear all;
strDataPath = 'D:\Data\Processed\DataKoen';
strDataPathPre = fullpath(strDataPath,'PreTraining');
strDataPathPost= fullpath(strDataPath,'PostTraining');

strFileFormat = '*Res.mat';
sFilesPre = dir(fullpath(strDataPathPre,strFileFormat));
boolRunTrained = 1; %trained (1=[1 2 4 5]) or not trained 0=[3 6])

if boolRunTrained
	vecKeepStim = [1 2 4 5];
else
	vecKeepStim = [3 6];
end
% get variables
for intSubject = 1:numel(sFilesPre)
strNamePre = sFilesPre(intSubject).name;
strFolderPre = sFilesPre(intSubject).folder;
%find post
cellSplit = strsplit(strNamePre,'_');
strSearchPostName = strcat(cellSplit{1},strFileFormat);
sFilesPost = dir(fullpath(strDataPathPost,strSearchPostName));

% pre
sLoadPre = load(fullpath(strFolderPre,strNamePre));
vecStimTypePre = sLoadPre.info.Stim.Parameters.StimSmp(1,:);
indKeepTrialsPre = ismember(vecStimTypePre,vecKeepStim);
vecStimTypePre = label2idx(vecStimTypePre(indKeepTrialsPre));
matActPre = sLoadPre.Res.CaSigCorrected(:,indKeepTrialsPre,:)-1;% (frames x trials x ROIs). Als je hier de waarde '1' vanaf trekt heb je in principe dF/F. 
vecAxT = sLoadPre.Res.ax;%. De as tijden, stimuli waren op t = 0 en duurden 1 seconde.
indStimFrames = vecAxT>0 & vecAxT<1;
%get data
matActNeuronByTrialPre = squeeze(mean(matActPre(indStimFrames,:,:),1))';
vecIsOccludedPre = sLoadPre.info.Stim.Parameters.StimSmp(2,indKeepTrialsPre);
intNeuronsPre = size(matActNeuronByTrialPre,1);

% post
sLoadPost = load(fullpath(sFilesPost(1).folder,sFilesPost(1).name));
vecStimTypePost = sLoadPost.info.Stim.Parameters.StimSmp(1,:);
indKeepTrialsPost = ismember(vecStimTypePost,vecKeepStim);
vecStimTypePost = label2idx(vecStimTypePost(indKeepTrialsPost));
matActPost = sLoadPost.Res.CaSigCorrected(:,indKeepTrialsPost,:)-1;% (frames x trials x ROIs). Als je hier de waarde '1' vanaf trekt heb je in principe dF/F. 
vecAxT = sLoadPost.Res.ax;%. De as tijden, stimuli waren op t = 0 en duurden 1 seconde.
indStimFrames = vecAxT>0 & vecAxT<1;
%get data
matActNeuronByTrialPost = squeeze(mean(matActPost(indStimFrames,:,:),1))';

vecIsOccludedPost = sLoadPost.info.Stim.Parameters.StimSmp(2,indKeepTrialsPost);
intNeuronsPost = size(matActNeuronByTrialPost,1);

intUseNeurons = min([intNeuronsPre intNeuronsPost]) - 10;
vecUseNeuronsPre = randperm(intNeuronsPre,intUseNeurons);
vecUseNeuronsPost = randperm(intNeuronsPost,intUseNeurons);

%% calculate all pairwise weight vectors
intStimNum = numel(unique(vecStimTypePre));
[varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(vecStimTypePre + vecIsOccludedPre*intStimNum); %1-6 is full, 7-12 is occ
intCombNum = (intStimNum^2-intStimNum)/2;
intTypeCV = 2;
dblLambda = 1/10;

matDataFullPre = matActNeuronByTrialPre(vecUseNeuronsPre,vecIsOccludedPre==0)';
vecStimTypeFullPre = vecStimTypePre(vecIsOccludedPre==0);
[dblPerformanceCV_FullPre,vecDecodedIndexCV_All,matPosteriorProbability_All,dblMeanErrorDegs_All,matConfusion_All,matWeights_All] = ...
	doCrossValidatedDecodingLR(matDataFullPre,vecStimTypeFullPre,intTypeCV,[],dblLambda);

matDataFullPost = matActNeuronByTrialPost(vecUseNeuronsPost,vecIsOccludedPost==0)';
vecStimTypeFullPost = vecStimTypePost(vecIsOccludedPost==0);
[dblPerformanceCV_FullPost,vecDecodedIndexCV_All,matPosteriorProbability_All,dblMeanErrorDegs_All,matConfusion_All,matWeights_All] = ...
	doCrossValidatedDecodingLR(matDataFullPost,vecStimTypeFullPost,intTypeCV,[],dblLambda);

matDataOccPre = matActNeuronByTrialPre(vecUseNeuronsPre,vecIsOccludedPre==1)';
vecStimTypeOccPre = vecStimTypePre(vecIsOccludedPre==1);
[dblPerformanceCV_OccPre,vecDecodedIndexCV_All,matPosteriorProbability_All,dblMeanErrorDegs_All,matConfusion_All,matWeights_All] = ...
	doCrossValidatedDecodingLR(matDataOccPre,vecStimTypeOccPre,intTypeCV,[],dblLambda);

matDataOccPost = matActNeuronByTrialPost(vecUseNeuronsPost,vecIsOccludedPost==1)';
vecStimTypeOccPost = vecStimTypePost(vecIsOccludedPost==1);
[dblPerformanceCV_OccPost,vecDecodedIndexCV_All,matPosteriorProbability_All,dblMeanErrorDegs_All,matConfusion_All,matWeights_All] = ...
	doCrossValidatedDecodingLR(matDataOccPost,vecStimTypeOccPost,intTypeCV,[],dblLambda);

%% pre-allocate
matPairPerf_FullPre = nan(intStimNum,intStimNum);
matW_FullPre = nan(intCombNum,intUseNeurons);
vecBias_FullPre = nan(intCombNum,1);

matPairPerf_FullPost = nan(intStimNum,intStimNum);
matW_FullPost = nan(intCombNum,intUseNeurons);
vecBias_FullPost = nan(intCombNum,1);

matPairPerf_OccPre = nan(intStimNum,intStimNum);
matW_OccPre = nan(intCombNum,intUseNeurons);
vecBias_OccPre = nan(intCombNum,1);

matPairPerf_OccPost = nan(intStimNum,intStimNum);
matW_OccPost = nan(intCombNum,intUseNeurons);
vecBias_OccPost = nan(intCombNum,1);

%% run
intComb = 0;
for intS1=1:intStimNum
	for intS2=(intS1+1):intStimNum
		intComb = intComb + 1;
		%% pairwise weight vectors
		%full-pre
		indTrialsFullPre = vecStimTypeFullPre==intS1|vecStimTypeFullPre==intS2;
		[dblPerformanceFullPre,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matDataFullPre(indTrialsFullPre,:),vecStimTypeFullPre(indTrialsFullPre),intTypeCV,[],1);
		matW_FullPre(intComb,:) = matWeights(1:(end-1),1); %remove bias term
		vecBias_FullPre(intComb) = matWeights(end,1);
		matPairPerf_FullPre(intS1,intS2) = dblPerformanceFullPre;
		
		%occ-pre
		indTrialsOccPre = vecStimTypeOccPre==intS1|vecStimTypeOccPre==intS2;
		[dblPerformanceOccPre,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matDataOccPre(indTrialsOccPre,:),vecStimTypeOccPre(indTrialsOccPre),intTypeCV,[],1);
		matW_OccPre(intComb,:) = matWeights(1:(end-1),1); %remove bias term
		vecBias_OccPre(intComb) = matWeights(end,1);
		matPairPerf_OccPre(intS1,intS2) = dblPerformanceOccPre;
		
		%full-pre
		indTrialsFullPost = vecStimTypeFullPost==intS1|vecStimTypeFullPost==intS2;
		[dblPerformanceFullPost,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matDataFullPost(indTrialsFullPost,:),vecStimTypeFullPost(indTrialsFullPost),intTypeCV,[],1);
		matW_FullPost(intComb,:) = matWeights(1:(end-1),1); %remove bias term
		vecBias_FullPost(intComb) = matWeights(end,1);
		matPairPerf_FullPost(intS1,intS2) = dblPerformanceFullPost;
		
		%full-pre
		indTrialsOccPost = vecStimTypeOccPost==intS1|vecStimTypeOccPost==intS2;
		[dblPerformancOccPost,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matDataOccPost(indTrialsOccPost,:),vecStimTypeOccPost(indTrialsOccPost),intTypeCV,[],1);
		matW_OccPost(intComb,:) = matWeights(1:(end-1),1); %remove bias term
		vecBias_OccPost(intComb) = matWeights(end,1);
		matPairPerf_OccPost(intS1,intS2) = dblPerformancOccPost;
		
	end
end

%% q1: are discriminating dimensions lying in low-dimensional subspace?
%max dim here is probably n-1 = 5 with 6 stimuli

%matW_FullPre matW_OccPre matW_FullPost matW_OccPost
matX = matW_FullPre;
matY = matW_OccPre;
vecR2A = nan(1,intCombNum);
dblLambda  =1/1000;
for intRank = 1:intCombNum
[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR(matX,matY,intRank,dblLambda);
vecR2C(intRank) = sSuppOut.dblR2;

[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR2(matX,matY,intRank,dblLambda);
vecR2D(intRank) = sSuppOut.dblR2;

end

%shuffle
matX_shuff = matX;
matY_shuff = matY;
for intP=1:size(matX,2)
	matX_shuff(:,intP) = circshift(matX(:,intP),[intP-1 0]);
	matY_shuff(:,intP) = circshift(matY(:,intP),[intP-1 0]);
end
for intRank = 1:intCombNum
[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR(matX_shuff,matY_shuff,intRank,dblLambda);
vecR2CS(intRank) = sSuppOut.dblR2;

[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR2(matX_shuff,matY_shuff,intRank,dblLambda);
vecR2DS(intRank) = sSuppOut.dblR2;

end

figure;
hold on
plot(vecR2C,'b')
plot(vecR2CS,'b--')
plot(vecR2D,'r')
plot(vecR2DS,'r--')
hold off
title(sprintf('Animal %d, Pre',intSubject));

%% q2: does this subspace grow or shrink after training?

%matW_FullPre matW_OccPre matW_FullPost matW_OccPost
matX = matW_FullPost;
matY = matW_OccPost;
vecR2A = nan(1,intCombNum);
dblLambda  =1/1000;
for intRank = 1:intCombNum
[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR(matX,matY,intRank,dblLambda);
vecR2C(intRank) = sSuppOut.dblR2;

[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR2(matX,matY,intRank,dblLambda);
vecR2D(intRank) = sSuppOut.dblR2;

end

%shuffle
matX_shuff = matX;
matY_shuff = matY;
for intP=1:size(matX,2)
	matX_shuff(:,intP) = circshift(matX(:,intP),[intP-1 0]);
	matY_shuff(:,intP) = circshift(matY(:,intP),[intP-1 0]);
end
for intRank = 1:intCombNum
[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR(matX_shuff,matY_shuff,intRank,dblLambda);
vecR2CS(intRank) = sSuppOut.dblR2;

[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR2(matX_shuff,matY_shuff,intRank,dblLambda);
vecR2DS(intRank) = sSuppOut.dblR2;

end

figure;
hold on
plot(vecR2C,'b')
plot(vecR2CS,'b--')
plot(vecR2D,'r')
plot(vecR2DS,'r--')
hold off
title(sprintf('Animal %d, Post',intSubject));


%% q3: what is cosine similarity of weight vectors between occluded and full conditions?
matX = matW_FullPost;
matY = matW_OccPost;
vecAngSim_Pre = nan(1,intCombNum);
vecAngSim_Post = nan(1,intCombNum);
for intComb=1:intCombNum
	[dblCS,dblAS]=cossim(matW_FullPre(intComb,:),matW_OccPre(intComb,:));
	vecAngSim_Pre(intComb) = dblAS;
	
	[dblCS,dblAS]=cossim(matW_FullPost(intComb,:),matW_OccPost(intComb,:));
	vecAngSim_Post(intComb) = dblAS;
end
[h,pSimPre]=ttest(vecAngSim_Pre);
[h,pSimPost]=ttest(vecAngSim_Post);

%% q4: does full/occluded cosine similarity become higher/lower after training? => more/less generalization
vecChangeSim = vecAngSim_Post - vecAngSim_Pre;
dblMeanChange = mean(vecChangeSim);
[h,p]=ttest(vecChangeSim);

vecAggMeanSimPre(intSubject) = mean(vecAngSim_Pre)
vecAggMeanSimPost(intSubject) = mean(vecAngSim_Post)
vecAggChangeSim(intSubject) = dblMeanChange
end

[h,pPre]=ttest(vecAggMeanSimPre)
[h,pPost]=ttest(vecAggMeanSimPost)
[h,pDiff]=ttest2(vecAggMeanSimPost,vecAggMeanSimPre)
%% plot
figure
scatter(vecAggMeanSimPre,vecAggMeanSimPost);
hold on
plot([-0.1 0.1],[-0.1 0.1],'--','color',[0.7 0.7 0.7])
plot([0 0],[-0.1 0.1],'color',[0.7 0.7 0.7])
plot([-0.1 0.1],[0 0],'color',[0.7 0.7 0.7])
hold off
xlabel('Pre-training pop code sim full/occluded ')
ylabel('Post-training pop code sim  full/occluded ')
title(sprintf('Pre-training vs 0, p=%.3f; post, p=%.3f; diff,p=%.3f',pPre,pPost,pDiff))
ylim([-0.1 0.1]);
xlim([-0.1 0.1]);
fixfig;grid off

export_fig('MultiDimOcclusionAnalysis.jpg')