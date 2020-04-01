function dblPupilLuminance = getPupilColor(intFrame,intPupFrameRange,dblPosX,dblPosY,matMovie)
	
	%get frame range
	intStartFrame = round(intFrame-(0.5*intPupFrameRange));
	intStopFrame = round(intFrame+(0.5*intPupFrameRange)-1);
	intMaxFrame = size(matMovie,3);
	if intStartFrame<1
		intStartFrame=1;
		intStopFrame=intFrame+intPupFrameRange-1;
	elseif intStopFrame>intMaxFrame
		intStartFrame=intMaxFrame-intPupFrameRange+1;
		intStopFrame=intMaxFrame;
	end
	
	%get pixel locations surrounding middle of pupil
	matSelect = false(size(matMovie,1),size(matMovie,2));
	matSelect(round(dblPosY-2):round(dblPosY+2),round(dblPosX-2):round(dblPosX+2)) = true;
	matVals = nan(25,length(intStartFrame:intStopFrame));
	intCounter = 0;
	for intFrame=intStartFrame:intStopFrame
		intCounter = intCounter + 1;
		matFrame = matMovie(:,:,intFrame);
		matVals(:,intCounter) = matFrame(matSelect);
	end
	dblPupilLuminance = abs(mean(matVals(:)));
end