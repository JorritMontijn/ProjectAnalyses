function matSynchrony = calcSynchrony(cellSpikeTimes,dblTau)
	%calcSynchrony Calculates pairwise synchrony matrix from point events
	%   Syntax: matSynchrony = calcSynchrony(cellSpikeTimes)
	%Input: 
	%	cellSpikeTimes: 1-dimensional cell array with N elements, where 
	%		one element is a vector of point event times
	%	dblTau:	temporal kernel width for calculating the synchrony metric;
	%		Default: dblTau = 0.01
	%Output: 
	%	matSynchrony: [N x N] synchrony matrix for pairwise temporal
	%		synchronization values. The diagonal is proportional to the
	%		number of events (e.g., spike count) by a constant that depends
	%		on the kernel width dblTau
	%
	%	By Jorrit Montijn, 06-02-18 (dd-mm-yy; Universite de Geneve)
	
	%% define defaults
	if ~exist('dblTau','var') || isempty(dblTau)
		dblTau=0.01;
	end
	intNumN = numel(cellSpikeTimes);
	
	%% pre-allocate
	matSynchrony = nan(intNumN,intNumN);
	dblDivFac = normpdf(0,0,dblTau);
	hTic = tic;
	%dblLastMsg=-5;
	for intN1=1:1%intNumN
		vecTimes1 = cellSpikeTimes{intN1}';
		for intN2=intN1:1%intNumN
			vecTimes2 = cellSpikeTimes{intN2};
			intNum2 = numel(vecTimes2);
			dblSync = 0;
			parfor intSpike2=1:2000%intNum2
				
				dblTime2 = vecTimes2(intSpike2);
				vecDiffT = vecTimes1 - dblTime2;
				vecDiffT = vecDiffT(abs(vecDiffT)<(dblTau*3));
				vecDens = normpdf(vecDiffT,0,dblTau);
				dblSync = dblSync + sum(vecDens);
				
				
			end
			%if (toc(hTic) - dblLastMsg) > 5
				dblLastMsg = toc(hTic);
				fprintf('N1=%d/%d;N2=%d/%d [%s]\n',intN1,intNumN,intN2-intN1+1,intNumN-intN1+1,getTime);
			%end
			matSynchrony(intN1,intN2) = dblSync;
			matSynchrony(intN2,intN1) = dblSync;
				
		end
	end
	matSynchrony = matSynchrony./dblDivFac;
end
