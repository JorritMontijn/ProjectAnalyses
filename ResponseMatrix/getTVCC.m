function matTVCOV = getTVCC(matData)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	[intTimePoints,intNeurons] = size(matData);

	%initial parameters
	vecKernel = ones(1,intTimePoints)/intTimePoints;
	dblMu2 = sum(vecKernel .* (1:intTimePoints).^2);
	dbPsi = sum(vecKernel.^2);
	vecWeightFunction = ones(1,intTimePoints)/intTimePoints;
	
	%get cov matrix at all points t
	matCov = nan(intNeurons,intNeurons,intTimePoints);
	for intT=1:intTimePoints
		matCov(:,:,intT) = matData(intT,:)' * matData(intT,:);
	end
	matTVCOV = nan(intNeurons,intNeurons,intTimePoints);
	
	%run adaptive fitting
	for intT=1:intTimePoints
		for intNeuron1=1:(intNeurons-1)
			for intNeuron2=(intNeuron1+1):intNeurons
				%% get optimal bandwidth
				%initial parameters
				dblH = 1/intTimePoints;
				
				
				for i=1:10
					matTVCOV(intNeuron1,intNeuron2,intT) = (1/dblH) * sum(( vecKernel(round(((1:intTimePoints)-intT)/dblH) ) ) * (matCov(intNeuron1,intNeuron2,intT)));
					
					dblVar = ( 1 / (intTimePoints - 1) ) * sum( (matCov(intNeuron1,intNeuron2,2:end) -  matCov(intNeuron1,intNeuron2,1:(end-1)))^.2 );

					

					vecRii = (1/(intTimePoints * (dblH* intTimePoints^(1/10)).^3 )) * sum(vecKernel(((1:intTimePoints)-intT)/(dblH * intTimePoints.^(1/10))) * matCovStat(intNeuron1,intNeuron2));

					

					dblH =( ( dblVar * dbPsi * sum(vecWeightFunction) ) / ...
							( intTimePoints * dblMu2 .^2 * sum( vecWeightFunction .* vecRii.^2) )...
								).^(1/5);
				end
				
				
				
			end
		end
	end
	
	
	
end

