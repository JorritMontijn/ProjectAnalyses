function matEllipse = ellipseArea(matX,vecParams)
	%ellipseArea Creates logical matrix where all points within ellipse are
	%true and all others false
	%   Syntax: matY = ellipseArea(matX,vecParams)
	%	Inputs: 
	%	- matX: matrix of specified area size, values are discarded
	%	- vecParams: list of parameters;
	%		dblCenterX = vecParams(1);
	%		dblCenterY = vecParams(2);
	%		dblMajorAxisAngle = vecParams(3); (in radians)
	%		dblMajorRadius = vecParams(4);
	%		dblMinorRadius = vecParams(5);
	%
	%	Output:
	%	- matEllipse: logical matrix same size as matX, true if inside ellipse

	%get size of grid
	vecSize=size(matX)*2;
	
	%get parameters
	intCenterX = floor(max(min(vecParams(1),vecSize(2)),1));
	intCenterY = floor(max(min(vecParams(2),vecSize(1)),1));
	dblMajorAxisAngle = vecParams(3);
	dblMajorRadius = max(vecParams(4),1);
	dblMinorRadius = max(vecParams(5),1);
	
	%get size of grid & create grid
	[matGridX,matGridY] = meshgrid(1:vecSize(2),1:vecSize(1));
	matGridX = matGridX - vecSize(2)/2;
	matGridY = matGridY - vecSize(1)/2;
	
	matGridX = matGridX * (dblMinorRadius/dblMajorRadius);
	matEllipse = sqrt(matGridY.^2+matGridX.^2);
	
	
	
	%rotate & crop
	matEllipse = imrotate(matEllipse,rad2ang(dblMajorAxisAngle),'bilinear','crop');
	matEllipse = circshift(matEllipse,[intCenterY intCenterX]);
	matEllipse = matEllipse((end-vecSize(1)/2+1):end,(end-vecSize(2)/2+1):end);

	%binarize
	matEllipse = matEllipse <= dblMinorRadius;
	
	%subplot(2,1,2),imagesc(matEllipse),drawnow
end

