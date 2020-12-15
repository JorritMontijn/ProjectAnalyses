function [vecNoiseParallel,vecNoiseOrthogonal] = getNoiseParaOrtho(matSpikeCounts,vecOrientation)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
		
	[vecUniqueOris,dummy,vecOriIdx] = unique(vecOrientation);
	intTrials = numel(vecOrientation);
	intOriNum = numel(vecUniqueOris);

	vecNoiseParallel = nan(1,intTrials);
	vecNoiseOrthogonal = nan(1,intTrials);
	for intOriIdx=1:intOriNum
		%stim1
		vecTrialsOri1 = find(vecOriIdx==intOriIdx);
		intRepNum = numel(vecTrialsOri1);
		matAct1 = matSpikeCounts(:,vecTrialsOri1);
		vecMu1 = mean(matAct1,2);
		
		%stim2
		intOriIdx2 = modx(intOriIdx+1,intOriNum);
		matAct2 = matSpikeCounts(:,vecOriIdx==intOriIdx2);
		vecMu2 = mean(matAct2,2);
		
		%get ref & project
		vecRef = vecMu2 - vecMu1;
		vecRef = vecRef / norm(vecRef);
		matRef = matSpikeCounts  - vecMu1;
		
		%pre-alloc
		vecTempNoiseParallel = nan(1,intRepNum);
		vecTempNoiseOrthogonal = nan(1,intRepNum);
		parfor intRep = 1:intRepNum
			intTrial = vecTrialsOri1(intRep);
			vecPoint = matRef(:,intTrial);
			
			%project vecPoint onto vecRef and decompose in parallel & orthogonal
			[vecParallel,vecOrtho] = getProjectPointOnVector(vecRef,vecPoint);
			vecTempNoiseParallel(intRep) = norm(vecParallel);
			vecTempNoiseOrthogonal(intRep) = norm(vecOrtho);
		end
		
		%assign
		vecNoiseParallel(vecTrialsOri1) = vecTempNoiseParallel;
		vecNoiseOrthogonal(vecTrialsOri1) = vecTempNoiseOrthogonal;
	end
end

