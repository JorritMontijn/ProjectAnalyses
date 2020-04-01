function matCoincidence = calcCoincidence(cellSpikeTimes,vecWindow,numTau,dblStepSize)
	%calcCoincidence Calculates pairwise coincidence matrix from point events
	%   Syntax: matCoincidence = calcCoincidence(cellSpikeTimes)
	%Input: 
	%	cellSpikeTimes: 1-dimensional cell array with N elements, where 
	%		one element is a vector of point event times
	%	vecWindow: 2-element vector specifying start and stop time of
	%		calculation window
	%	[numTau]: temporal kernel width for calculating the coincidence metric;
	%		Default: dblTau = 0.01
	%	[dblStepSize]: time step size; only relevant in case of integer input;
	%		Default: dblStepSize = 0.0005
	%Output: 
	%	matCoincidence: [N x N] coincidence matrix for pairwise temporal
	%		coincidence values. The diagonal is normalized to be exactly
	%		the number of events (e.g., spike count) iff all events are 
	%		regularly interspersed in time such that their temporal
	%		summation becomes negligeable in the case of kernel width numTau
	%
	%Note: this function can use two different implementations, depending
	%on whether you supplied integer-type input or floating point input.
	%The integer-based version is somewhat faster than the float-based
	%version (around 10-50% depending on the number of trials and neurons). 
	%
	%	By Jorrit Montijn (Alex Pouget lab), 06-02-18 (dd-mm-yy; Universite de Geneve)
	
	%% define defaults
	boolInt = all(isint(cellSpikeTimes{1}));
	if ~exist('dblStepSize','var')
		dblStepSize = 0.5e-3;
	end
	if ~exist('numTau','var') || isempty(numTau)
		numTau=0.01;
		if boolInt
			numTau = cast(numTau/dblStepSize,class(cellSpikeTimes{1})); %#ok<ZEROLIKE>
		end
	end
	if boolInt && (~isint(numTau) || ~all(isint(vecWindow)))
		fprintf('Spike times are integer, but parameters (Tau and/or vecWindow) are not\n');
		error([mfilename ':IncompatibleTypes'],['Tau and spike times are inconsistent']);
	end
	intNumN = numel(cellSpikeTimes);
	
	%% pre-allocate
	matCoincidence = nan(intNumN,intNumN);
	dblDivFac = normpdf(0,0,double(numTau));
	
	%% pre-process event times; remove everything outside window
	cellSpikeTimesSmall = cell(size(cellSpikeTimes));
	for intN1=1:intNumN
		vecTimes1 = cellSpikeTimes{intN1};
		cellSpikeTimesSmall{intN1} = vecTimes1(vecTimes1 > vecWindow(1) & vecTimes1 < vecWindow(end));
	end
	
	%% run loop
	%hTic = tic;
	%dblLastMsg=-5;
	if boolInt
		%% run integer loop
		vecLookupTable = normpdf(double(0:(range(vecWindow)-1)),0,double(numTau));
		for intN1=1:intNumN
			vecTimes1 = cellSpikeTimesSmall{intN1}';
			vecN2s = intN1:intNumN;
			intNumN2 = numel(vecN2s);
			vecCoincidence = nan(1,intNumN2);
			for intN2=1:intNumN2
				vecTimes2 = cellSpikeTimesSmall{intN2+intN1-1};
				matDiffs = bsxfun(@minus,vecTimes1,vecTimes2);
				dblCoincidence = sum(sum(vecLookupTable(matDiffs(matDiffs>-1)+1)));
				vecCoincidence(intN2) = dblCoincidence;
			end
			matCoincidence(intN1,vecN2s) = vecCoincidence;
			matCoincidence(vecN2s,intN1) = vecCoincidence;
			
			%if (toc(hTic) - dblLastMsg) > 5
			%	dblLastMsg = toc(hTic);
			%	fprintf('N1=%d/%d [%s]\n',intN1,intNumN,getTime);
			%end
			
		end
	else
		%% run floating point loop
		for intN1=1:intNumN
			vecTimes1 = cellSpikeTimesSmall{intN1}';
			vecN2s = intN1:intNumN;
			intNumN2 = numel(vecN2s);
			vecCoincidence = nan(1,intNumN2);
			for intN2=1:intNumN2
				vecTimes2 = cellSpikeTimesSmall{intN2+intN1-1};
				dblCoincidence = sum(sum(normpdf(bsxfun(@minus,vecTimes1,vecTimes2),0,numTau)));
				vecCoincidence(intN2) = dblCoincidence;
			end
			matCoincidence(intN1,vecN2s) = vecCoincidence;
			matCoincidence(vecN2s,intN1) = vecCoincidence;
			
			%if (toc(hTic) - dblLastMsg) > 5
			%	dblLastMsg = toc(hTic);
			%	fprintf('N1=%d/%d [%s]\n',intN1,intNumN,getTime);
			%end
			
		end
	end
	%toc(hTic)
	matCoincidence = matCoincidence./dblDivFac;
end
%{
	%% run loop v0
	hTic = tic;
	dblLastMsg=-5;
	for intN1=1:intNumN
		vecTimes1 = cellSpikeTimes{intN1}';
		vecTimes1 = vecTimes1(vecTimes1 > dblStart & vecTimes1 < dblStop);
		for intN2=intN1:intNumN
			vecTimes2 = cellSpikeTimes{intN2};
			vecTimes2 = vecTimes2(vecTimes2 > dblStart & vecTimes2 < dblStop);
			intNum2 = numel(vecTimes2);
			dblSync = 0;
			for intSpike2=1:intNum2
				
				dblTime2 = vecTimes2(intSpike2);
				vecDiffT = vecTimes1 - dblTime2;
				vecDiffT = vecDiffT(abs(vecDiffT)<(dblTau*4));
				vecDens = normpdf(vecDiffT,0,dblTau);
				dblSync = dblSync + sum(vecDens);
				
				
			end
			if (toc(hTic) - dblLastMsg) > 5
				dblLastMsg = toc(hTic);
				fprintf('N1=%d/%d;N2=%d/%d [%s]\n',intN1,intNumN,intN2-intN1+1,intNumN-intN1+1,getTime);
			end
			matSynchrony(intN1,intN2) = dblSync;
			matSynchrony(intN2,intN1) = dblSync;
				
		end
	end
	toc(hTic)
	matSynchrony = matSynchrony./dblDivFac;

%}
