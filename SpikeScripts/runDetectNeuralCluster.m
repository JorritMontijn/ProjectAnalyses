%% runDetectNeuralCluster
%{
search for frames where multiple neurons are
active; or where other neurons are active in the subsequent frame;
keep adding frames until no neurons are active at a certain point;
then call this an assembly event and save the event with which
neurons are active at which point

input: matSpikeCountsHD [neuron x time-window] spiking matrix
%}
%% generate fake data
clear all;
boolSavePlots = false;
%neural data
intNeurons = 50;
intFrames = 10000;
dblSamplingFreq = 25;
dblTotLengthSecs = intFrames/dblSamplingFreq;
matSpikeCountsHD=(rand(intNeurons,intFrames) - 0.99) > 0;
%stimuli

%% run analysis
fprintf('Retrieving population spiking data for assembly formation analysis...   [%s]\n',getTime);

%get number of active neurons per frame
vecNrActive=sum(matSpikeCountsHD>0,1);

%select which frames are potentially interesting
vecPotentialAssemblyFrames = vecNrActive > 0;

%remove frames where single neurons fire single spikes in temporal isolation
vecSingleSpikes = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))<2 & vecNrActive(2:(end-1))>0 & vecNrActive(3:end)==0;
vecPotentialAssemblyFrames([false vecSingleSpikes false]) = 0;

%get frames where multiple neurons are firing in temporal isolation
vecTransientAssemblies = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))>2 & vecNrActive(3:end)==0;
vecPotentialAssemblyFrames([false vecTransientAssemblies false]) = 1;

%define start frames of assemblies
vecAssemblyStarts = [false vecPotentialAssemblyFrames(1:(end-1))==0 & vecPotentialAssemblyFrames(2:(end))>0];
vecAssemblyStops = [false vecPotentialAssemblyFrames(1:(end-1))>0 & vecPotentialAssemblyFrames(2:end)==0];

%combine
vecAssemblyStarts = find(vecAssemblyStarts);
vecAssemblyStops = find(vecAssemblyStops);
if vecAssemblyStarts(1) > vecAssemblyStops(1)
	vecAssemblyStarts = [1 vecAssemblyStarts];
end
[dummy,vecReorder] = sort(vecAssemblyStarts,'ascend');
if length(vecAssemblyStops) < length(vecAssemblyStarts)
	vecAssemblyStops(end+1) = length(vecPotentialAssemblyFrames);
end
vecAssemblyStarts = vecAssemblyStarts(vecReorder);
vecAssemblyStops = vecAssemblyStops(vecReorder);

%get assemblies from spike matrix + assembly stop/starts
matAssemblies = getAssemblies(matSpikeCountsHD,vecAssemblyStarts,vecAssemblyStops)>0;

%remove assemblies with 2 or fewer neurons
indRealAssemblies = sum(matAssemblies,1)>=5 & sum(matAssemblies,1) < intNeurons;
matAssemblies = matAssemblies(:,indRealAssemblies);
vecAssemblyStarts = vecAssemblyStarts(indRealAssemblies);
vecAssemblyStops = vecAssemblyStops(indRealAssemblies);

%get consistencies
matNeuronCorrelations = corr(matAssemblies');
matAssemblyCorrelations = corr(matAssemblies);

%cluster assemblies
matSelect = tril(true(size(matAssemblyCorrelations)),0) & triu(true(size(matAssemblyCorrelations)),0);
matDistAN = 1-matAssemblyCorrelations;
matDistAN(matSelect) = 0;
[intNumberOfAssemblies,vecSilhouetteDistances] = doFastClustering(matDistAN,min([100 round(sqrt(length(matDistAN))*2)]));
fprintf('Optimal number of assemblies: %d. Proceeding with assembly assignment... [%s]\n',intNumberOfAssemblies,getTime);

%re-cluster with detected number of assemblies
matLinkageA = linkage(matDistAN,'ward');
vecA = cluster(matLinkageA,'maxclust',intNumberOfAssemblies);
[vecAssemblyIdentity,vecReorderA] = sort(vecA,'ascend');

%return
intT = length(0:(1/dblSamplingFreq):(dblTotLengthSecs+5));
matAssemblyActivity = zeros(intNumberOfAssemblies,intT);
matAssemblyActivity = putAssemblies(matAssemblyActivity,vecA,vecAssemblyStarts,vecAssemblyStops);

%% plot
fprintf('Assemblies assigned. Creating figure and dendrogram... [%s]\n',getTime);

figure
subplot(2,2,1)
matLinkageAN = linkage(matDistAN,'ward');

[H,T,vecClusterReorderAN] = dendrogram(matLinkageAN,0);
drawnow;

subplot(2,2,2)
vecFiltSil = conv(vecSilhouetteDistances,normpdf(-2:2,0,0.8),'same');
plot(vecSilhouetteDistances);
xlabel('Number of clusters');
ylabel('Silhouette distance');
title(sprintf('Number of clusters that is optimal; %d [val=%.3f]',intNumberOfAssemblies,vecSilhouetteDistances(intNumberOfAssemblies)));

subplot(2,2,3)
imagesc(matAssemblyCorrelations(vecReorderA,vecReorderA),[-1 1]);colormap('redblue');freezeColors;
title('Assembly correlations, ordered by cluster')

subplot(2,4,7)
imagesc(vecAssemblyIdentity);colormap('jet');freezeColors;
title('Cluster membership')

subplot(2,4,8)
matM = getBlockMeans(matAssemblyCorrelations(vecReorderA,vecReorderA),vecAssemblyIdentity);
vecAssemblyAutoCorr = diag(matM);
scatter(1:length(vecAssemblyAutoCorr),vecAssemblyAutoCorr');
ylim([0 1])
title('Assembly consistency')

drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
if boolSavePlots
	strFig = sprintf('%sagg_assembly_correlations_pop%d_raw',strSes,intPopulation);
	export_fig([sParams.strFigDir strFig '.tif']);
	export_fig([sParams.strFigDir strFig '.pdf']);
end
