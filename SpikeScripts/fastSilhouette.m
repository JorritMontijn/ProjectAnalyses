function [vecDistClusterSilhouettes,matDistC] = fastSilhouette(matDist,vecC)
	%fastSilhouette Calculates fast silhouette distance
	%   [vecDistClusterSilhouettes,matDistC] = fastSilhouette(matDist,vecC)
	
	%get cluster-based distance; [intC x intC] matrix, for all within and across
	if size(vecC,2) == 1,vecC = vecC';end
	intMaxC = length(unique(vecC));
	matDistC = getBlockMeans(matDist,vecC);
	
	%cluster-based distance; [intC x intC] matrix, for all within and across
	%cluster distances
	vecWithinDists = diag(matDistC);
	vecDistClusterSilhouettes = nan(1,intMaxC);
	vecAllC = true(1,intMaxC);
	for intC=1:intMaxC
		vecSelect = vecAllC;
		vecSelect(intC) = false;
		vecDistClusterSilhouettes(intC) = (min(matDistC(intC,vecSelect)) - vecWithinDists(intC)) ./ max(vecWithinDists(intC),min(matDistC(intC,vecSelect)));
	end
	%where a(i) is the average distance within the ith cluster, and
	%matDistC(i,k) is  the average distance from the ith cluster to another
	%cluster k.
end

