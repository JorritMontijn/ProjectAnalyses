function [dblCorr,dblCorr1,vecProjectedLocation,matProjectedPoints,dblCorr2,vecProjectedLocation2,matProjectedPoints2] = getMeanRFcorrWithAngleConstraint(dblAngle,vecX,vecY,vecV1,vecV2)
	%getMeanRFcorrWithAngleConstraint Correlation of vecV with location of X/Y projected unto a line with angle dblAngle
	%   dblCorr = getMeanRFcorrWithAngleConstraint(dblAngle,vecX,vecY,vecV1,vecV2)
	
	%center x/y
	%vecX=zscore(vecX(:));
	%vecY=zscore(vecY(:));
	matXY=[vecX vecY]';
	%rotate reference vector
	matRot = [cos(dblAngle) sin(dblAngle);...
		-sin(dblAngle) cos(dblAngle)];
	vecRefVector=[1;0];
	vecRotRef = matRot * vecRefVector;
	
	%calc corr
	[vecProjectedLocation,matProjectedPoints] = getProjOnLine(matXY,vecRotRef);
	dblCorr1 = corr(vecProjectedLocation,vecV1);
	
	%rotate 2nd reference vector
	dblAngle2 = dblAngle + pi/2;
	matRot2 = [cos(dblAngle2) sin(dblAngle2);...
		-sin(dblAngle2) cos(dblAngle2)];
	vecRefVector=[1;0];
	vecRotRef2 = matRot2 * vecRefVector;
	
	%calc 2nd corr
	[vecProjectedLocation2,matProjectedPoints2] = getProjOnLine(matXY,vecRotRef2);
	
	%calc correlation
	dblCorr2 = corr(vecProjectedLocation2,vecV2);
	
	%tot corr
	dblCorr = (dblCorr1 + sign(dblCorr1)*abs(dblCorr2))/2;
end
