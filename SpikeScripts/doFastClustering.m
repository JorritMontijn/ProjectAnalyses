function [intNumberOfClustersThatIsOptimal,vecSilhouetteDistances] = doFastClustering(matDist,intN)
	%doFastClustering Performs fast clustering
	%   [intNumberOfClustersThatIsOptimal,vecSilhouetteDistances] = doFastClustering(matDist,intN)
	%Inputs:
	% - matDist: matrix of distances
	% - intN: maximum number of clusters to calculate
	
	if nargin < 2 || intN < 2 || intN > (size(matDist,1)-1)
		intN = 15;
	end
	
	%get silhouette distances
	matLinkage = linkage(matDist,'ward');
	vecSilhouetteDistances = zeros(1,intN);
	
	for intNumberOfClusters=2:intN
		%fprintf('Cluster size %d\n',intNumberOfClusters);
		vecT = cluster(matLinkage,'maxclust',intNumberOfClusters);
		vecSilhouetteDistances(intNumberOfClusters) = nanmean(fastSilhouette(matDist,vecT));
	end
	
	%smooth and get largest relative peak
	vecFiltSil = conv(vecSilhouetteDistances,normpdf(-2:2,0,1),'same');
	vecDiff = diff(vecFiltSil);
	
	%get troughs/peaks
	[vecTroughs,vecTroughIdx] = findpeaks(-vecDiff);
	[vecPeaks,vecPeakIdx] = findpeaks(vecDiff);
	vecLocalBasins = nan(size(vecPeaks));
	for intTrough = 1:length(vecTroughs)
		intTroughIdx = vecTroughIdx(intTrough);
		vecLeftPeaks = vecPeakIdx(vecPeakIdx<intTroughIdx);
		if isempty(vecLeftPeaks)
			dblLeftPeak = vecDiff(1);
		else
			dblLeftPeak = vecDiff(vecLeftPeaks(end));
		end
		vecRightPeaks = vecPeakIdx(vecPeakIdx>intTroughIdx);
		if isempty(vecRightPeaks)
			dblRightPeak = vecDiff(end);
		else
			dblRightPeak = vecDiff(vecRightPeaks(1));
		end
		vecLocalBasins(intTrough) = min([dblLeftPeak dblRightPeak]-vecDiff(intTroughIdx));
	end
	[dblVal,intIdx] = max(vecLocalBasins);
	intNumberOfClustersThatIsOptimal = vecTroughIdx(intIdx)-1;
end

