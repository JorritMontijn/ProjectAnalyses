function [vecPopSync,vecTimestamps,intJitterBins] = getPopulationSynchronization(cellSpikeTimes,dblJitter_ms,dblBinSize_ms_OR_vecTimestamps)
	%getPopulationSynchronization Detects population-wide synchronization events
	%   Syntax: matSynchrony = calcSynchrony(cellSpikeTimes)
	%Input: 
	%	cellSpikeTimes: 1-dimensional cell array with N elements, where 
	%		one element is a vector of point event times
	%	[dblTau_ms]: temporal integration kernel (default: 5 ms)
	%	[dblBinSize_ms_OR_vecTimestamps]: temporal bin size (default: 1 ms)	/ timestamp vector
	%
	%Output: 
	%	vecPopSync: vector with number of active neurons per time bin
	%	vecTimestamps: vector with time stamps per bin (centers)
	%	intJitterBins: number of bins within which an event gets merged
	%
	%Version History:
	%2019-07-10 Created getPopulationSynchronization function [by Jorrit Montijn]
	%			Inspired by Shahidi, Andrei, Hu & Dragoi (2019), Nat Neurosci
	%			Inspired by Pipa, Wheeler, Singer & Nikolic (2008), J Comp Neurosci
	
	%% set parameters
	if ~exist('dblBinSize_ms_OR_vecTimestamps','var') || isempty(dblBinSize_ms_OR_vecTimestamps),dblBinSize_ms_OR_vecTimestamps = 1;end
	if ~exist('dblJitter_ms','var') || isempty(dblJitter_ms),dblJitter_ms = 5;end
	
	%% transform parameters
	ptrTic = tic;
	dblJitterSecs = dblJitter_ms/1000;
	if numel(dblBinSize_ms_OR_vecTimestamps) == 1
		%bin size
		dblBinSizeSecs = dblBinSize_ms_OR_vecTimestamps/1000;
		
		%build timestamps
		dblFirstSpike = min(cellfun(@min,cellSpikeTimes));
		dblLastSpike = max(cellfun(@max,cellSpikeTimes));
		dblStart = dblFirstSpike - dblBinSizeSecs - dblJitterSecs;
		dblStop = dblLastSpike + dblBinSizeSecs + dblJitterSecs*2;
		vecBinEdges = dblStart:dblBinSizeSecs:dblStop;
		vecTimestamps = vecBinEdges(2:end) - dblBinSizeSecs/2;
	else
		%timestamps
		vecTimestamps = dblBinSize_ms_OR_vecTimestamps(:)';
		dblBinSizeSecs = median(diff(vecTimestamps));
		
		%build edges
		vecBinEdges = [vecTimestamps(1)-(dblBinSizeSecs/2) (vecTimestamps + (dblBinSizeSecs/2))];
	end
	
	%% build spiking index matrix
	intMaxN = numel(cellSpikeTimes);
	intJitterBins = round(dblJitterSecs/dblBinSizeSecs)*2 + 1;
	vecKernel = true(1,intJitterBins);
	matSpikes = false(intMaxN,numel(vecTimestamps));
	for intN=1:intMaxN
		if toc(ptrTic) > 10
			fprintf('Now at neuron %d/%d [%s]\n',intN,intMaxN,getTime);
			ptrTic = tic;
		end
		matSpikes(intN,:) = conv(logical(histcounts(cellSpikeTimes{intN},vecBinEdges)),vecKernel,'same');
	end
	vecPopSync = sum(matSpikes,1);
end

