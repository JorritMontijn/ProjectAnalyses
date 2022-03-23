function [vecProjectedLocation,matProjectedPoints] = getProjOnLine(matPoints,vecRef,intIters)
	%getProjOnLine Summary of this function goes here
	%   [vecProjectedNorm,matProjectedPoints] = getProjOnLine(matPoints,vecRef,intIters)
	
	%get data
	intD=size(matPoints,1);
	intPoints = size(matPoints,2);
	if intPoints < intD
%		error([mfilename ':WrongDims'],'Number of dimensions is larger than number of points; please make sure matrix is in form [Trials x Neurons]');
	end
	
	%recenter
	matProj = ((vecRef*vecRef')/(vecRef'*vecRef));
	vecNormRef = vecRef/norm(vecRef);
	
	%calculate projected points
	matProjectedPoints = nan(size(matPoints));
	vecProjectedLocation = nan(size(matPoints,2),1);
	for intTrial=1:size(matPoints,2)
		vecPoint = matPoints(:,intTrial);
		vecOrth = matProj*vecPoint;
		matProjectedPoints(:,intTrial) = vecOrth;
		vecProjectedLocation(intTrial) = vecOrth'/vecNormRef';
	end
	
end
