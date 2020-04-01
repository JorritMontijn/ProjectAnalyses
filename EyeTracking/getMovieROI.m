function vecRect = getMovieROI(strFile,vecPrevRect)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%get video data
	sVideo1 = mmread(strFile,1000);
	sVideo2 = mmread(strFile,15000);

	%display figure and ask for rectangle
	h=figure;
	matRGB = imadjust(sVideo1.frames.cdata(:,:,1));
	matRGB(:,:,2) = imadjust(sVideo2.frames.cdata(:,:,1));
	matRGB(:,:,3) = zeros(size(matRGB(:,:,1)));
	imshow(matRGB);
	
	if isempty(vecPrevRect)
		h2 = imrect;
	else
		h2 = imrect(gca,vecPrevRect);
	end
	
	vecRect = wait(h2);
	
	close(h);
	drawnow;
end

