function [intNumberOfClustersThatIsOptimal,vecSilhouetteDistances] = doClustering(matDist,intN)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	if nargin < 2 || intN < 2 || intN > (size(matDist,1)-1)
		intN = 10;
	end
	
	matLinkage = linkage(matDist,'ward');
	vecSilhouetteDistances = zeros(1,intN);
	
	for intNumberOfClusters=2:intN
		fprintf('Cluster size %d\n',intNumberOfClusters);
		tic
		vecT = cluster(matLinkage,'maxclust',intNumberOfClusters);
		toc
		tic
		vecSilhouetteDistances(intNumberOfClusters) = mean(silhouette(matDist,vecT));
		toc
	end
	
	[dblSilhouetteValueOfTheNumberOfClustersThatIsOptimal,intNumberOfClustersThatIsOptimal] = max(vecSilhouetteDistances);
end

