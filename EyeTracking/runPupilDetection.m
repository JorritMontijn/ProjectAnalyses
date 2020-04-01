%% PUPIL DETECTION LOG
%{
offsets present:
20140207xyt01,xyt03
20140530xyt04
20140604xyt01-05
20140711xyt02-08
20140715xyt01-04

offsets missing:
20140207xyt02,xyt04-08
20140314xyt02-07
20140425xyt01-08
20140507xyt01-03
20140530xyt01-03,xyt05-07

pupil detection
20140207
xyt01-3: good
xyt04-8: bad

20140314
xyt02-07: good

20140425
xyt01-08: good

20140507
xyt01-03: good

20140530
xyt01-07: good

20140604
xyt01-03: good
xyt04-05: okay

20140711
xyt02-08: good

20140715
xyt01-04: good
%}

%% time selection
%define data
strSourceDir = 'D:\Data\Processed\imagingvideo\';
%sPupil = [];
for intMouse = []
	if intMouse == 1
		strSes = '20140207';
		vecRecordings = 1:8;%'*xyt01.avi'
	elseif intMouse == 2
		strSes = '20140314';
		vecRecordings = 1:7;
	elseif intMouse == 3
		strSes = '20140425';
		vecRecordings = 1:8;
	elseif intMouse == 4
		strSes = '20140507';
		vecRecordings = 1:3;
	elseif intMouse == 5
		strSes = '20140530';
		vecRecordings = 1:7;
	elseif intMouse == 6
		strSes = '20140604';
		vecRecordings = 1:5;
	elseif intMouse == 7
		strSes = '20140711';
		vecRecordings = 2:8;
	elseif intMouse == 8
		strSes = '20140715';
		vecRecordings = 1:4;
	end
	
	for intRecording=vecRecordings
		%% prep recording
		%load movie data
		strSourceFile = sprintf('croppedvideo_%sxyt%02d.mat',strSes,intRecording);
		sLoad = load([strSourceDir strSes filesep strSourceFile]);
		sVideoAll = sLoad.sVideoAll;
		clear sLoad;
		
		%detect start + end; define step size
		%ask which cut-off to use
		
		%display figure and ask for rectangle
		vecMean = squeeze(mean(mean(sVideoAll.matMovie,1),2));
		h=figure;
		plot(vecMean);
		xlim([0 length(vecMean)+1]);
		title(sprintf('Time-point selection for Recording %sxyt%02d',strSes,intRecording));
		
		%pre-alloc
		dblCutOff = -1;
		hX1 = [];
		hX2 = [];
		hX3 = [];
		boolAccept = false;
		while ~boolAccept
			
			%ask for cut off
			[dummy,dblCutOff,intButton] = ginput(1);
			
			%check if done
			if intButton ~= 1 && dblCutOff ~= -1
				boolAccept = true;
				break
			end
			
			%take point in-between two means as cut-off point, tranform to logical and
			%get time points of largest bright epoch
			indBright = vecMean > dblCutOff;
			sCC = bwconncomp(indBright);
			[intSize,intIndex] = max(cellfun(@numel,sCC.PixelIdxList));
			intStart = sCC.PixelIdxList{intIndex}(1);
			intStop = sCC.PixelIdxList{intIndex}(end);
			
			%remove old lines
			if ~isempty(hX1)
				delete(hX1);
				delete(hX2);
				delete(hX3);
				hX1 = [];
				hX2 = [];
				hX3 = [];
			end
			
			%plot selection
			hX1 = line(get(gca,'XLim'),[dblCutOff dblCutOff],'Color','k','LineWidth',1,'LineStyle','--');
			hX2 = line([intStart intStart]-0.5,get(gca,'YLim'),'Color','g','LineWidth',1,'LineStyle','--');
			hX3 = line([intStop intStop]+0.5,get(gca,'YLim'),'Color','r','LineWidth',1,'LineStyle','--');
		end
		
		%crop movie to selected time points
		matMovie = uint8(sVideoAll.matMovie(:,:,intStart:intStop));
		clear sVideoAll;
		
		%export figure
		strFig = [strSourceDir strSes filesep sprintf('autoTimeSelection%sxyt%02d',strSes,intRecording)];
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
		
		
		
		%% pupil center input
		%request pupil center for initial pupil color definition
		intPos = numel(sPupil);
		h2 = figure;
		image(matMovie(:,:,100));
		colormap(bone(255));
		vecBaseLocation = round(ginput(1));
		intPupilLuminance = matMovie(vecBaseLocation(2),vecBaseLocation(1),100);
		
		%save ROI data
		strMovieData = [strSourceDir strSes filesep sprintf('moviedata_%sxyt%02d',strSes,intRecording)];
		sPupil(intPos+1).strSes = strSes;
		sPupil(intPos+1).strRec = sprintf('xyt%02d',intRecording);
		sPupil(intPos+1).intRec = intRecording;
		sPupil(intPos+1).strMovieData = [strMovieData '.mat'];
		sPupil(intPos+1).vecBaseLocation = vecBaseLocation;
		sPupil(intPos+1).intPupilLuminance = intPupilLuminance;
		
		%close figures
		close(h);
		close(h2);
		drawnow;
		
		%save cut movie % meta data
		save(strMovieData,'matMovie','-v7.3');
		strFilePupil = ['MoviePupilData_' date];
		save(sprintf('%s%s',strSourceDir,strFilePupil),'sPupil');
		fprintf('\nSaved movie & meta data of %sxyt%02d to %s\n',strSes,intRecording,strSourceDir);
		
	end
end

%% pupil detection
clearvars -except sPupil strSourceDir
for intMovie=48%:numel(sPupil)
	%% run per recording
	%get data
	strSes = sPupil(intMovie).strSes;
	strRec = sPupil(intMovie).strRec;
	intRec	= sPupil(intMovie).intRec;
	strMovieData = sPupil(intMovie).strMovieData;
	vecBaseLocation = sPupil(intMovie).vecBaseLocation;
	dblPupilLuminance = im2double(sPupil(intMovie).intPupilLuminance);
	sLoad = load(strMovieData);
	matMovie = im2double(sLoad.matMovie);
	clear sLoad;
	
	%pre-allocate output
	intFrameMax = size(matMovie,3);
	vecPupilLuminance = nan(1,intFrameMax);
	vecWeight = nan(1,intFrameMax);
	vecRoundness = nan(1,intFrameMax);
	vecPosX = nan(1,intFrameMax);
	vecPosY = nan(1,intFrameMax);
	vecArea = zeros(1,intFrameMax);
	matMask = false(size(matMovie));
	
	%define loop variables
	dblMaxDist = sqrt(size(matMovie,1)^2+size(matMovie,2)^2);
	dblPrevPupLum = dblPupilLuminance;
	dblThisPupLum = dblPupilLuminance;
	dblNextPupLum = dblPupilLuminance;
	dblPosX = vecBaseLocation(1);
	dblPosY = vecBaseLocation(2);
	vecPosX(1) = vecBaseLocation(1);
	vecPosY(2) = vecBaseLocation(2);
	intPupFrameRange = 15;
	
	%blur all frames
	fprintf('Blurring movie from %sxyt%02d [%s]\n',strSes,intRec,getTime);
	matMovieBlurred = nan(size(matMovie));
	matFilter = ones(3,3)/9;
	for intFrame=1:intFrameMax
		matMovieBlurred(:,:,intFrame) = conv2(conv2(conv2(matMovie(:,:,intFrame),matFilter,'same'),matFilter,'same'),matFilter,'same');
	end
	fprintf('Starting processing of %sxyt%02d [%s]\n',strSes,intRec,getTime);
	
	%loop through frames for detection
	vecPupilLuminance(1) = getPupilColor(1,intPupFrameRange,dblPosX,dblPosY,matMovieBlurred);
	for intFrame=1:intFrameMax
		%msg
		if mod(intFrame,100) == 0
			fprintf('Now at frame %d/%d for %sxyt%02d [%s]\n',intFrame,intFrameMax,strSes,intRec,getTime);
		end
		
		%update pupil color over sliding window
		if isnan(dblPosX)
			dblPosX = nanmean(vecPosX);
			dblPosY = nanmean(vecPosY);
		end
		
		%get pupil color
		dblNextPupLum = getPupilColor(intFrame,intPupFrameRange,dblPosX,dblPosY,matMovieBlurred);
		
		if intFrame == 1
			dblThisPupLum = vecPupilLuminance(1);
		else
			dblThisPupLum = mean(vecPupilLuminance(1:(intFrame-1)));
		end
		
		
		%blur image thrice
		matThis = mean(cat(4,matMovieBlurred(:,:,max([1 intFrame-1])),matMovieBlurred(:,:,intFrame),matMovieBlurred(:,:,min([intFrameMax intFrame+1]))),4);
		
		%select potential eye pixels
		matEye = double(matThis>= 0.65*dblThisPupLum & matThis<= 1.6*dblThisPupLum);
		
		%remove small specks
		objSE = strel('ball',20,20);
		matEye = mat2gray(imerode(imerode(matEye,objSE),objSE));
		matEyeBW = bwareaopen(matEye>0,round(numel(matEye)*0.033));
		
		matEyeBW = bwmorph(matEyeBW,'skel',3);
		matEyeBW = bwmorph(matEyeBW,'majority',15);
		matEyeBW = bwareaopen(matEyeBW,round(numel(matEye)*0.01));
		
		%get locations of white pixels, remove all that are >3sds away from
		%mean, , then get
		%mean x+y and total number of pixels
		[vecLocY,vecLocX]=find(matEyeBW);
		dblMeanY = mean(vecLocY);
		dblMeanX = mean(vecLocX);
		dblSDY = std(vecLocY);
		dblSDX = std(vecLocX);
		[matX,matY] = meshgrid(1:size(matEyeBW,2),1:size(matEyeBW,1));
		matSDX = (matX-dblMeanX)/dblSDX;
		matSDY = (matY-dblMeanY)/dblSDY;
		matTotSD = sqrt(matSDX.*matSDX + matSDY.*matSDY);
		matEyeBW = matEyeBW & matTotSD<3;
		
		%select area that is large and close to previous location
		sCC = bwconncomp(matEyeBW);
		%get sizes
		vecSizes = cellfun(@numel,sCC.PixelIdxList);
		if numel(vecSizes) > 0
			vecBlobX = zeros(size(vecSizes));
			vecBlobY = zeros(size(vecSizes));
			
			%get roundness
			sP = regionprops(sCC, 'Perimeter');
			vecPerimeter = cell2mat(struct2cell(sP));
			vecBlobRoundness = 4*pi*(vecSizes./(vecPerimeter.^2));
			
			%get locations
			for intBlob=1:length(vecSizes)
				%select area
				matEyeBW = false(size(matEyeBW));
				matEyeBW(sCC.PixelIdxList{intBlob}) = true;
				
				%find pupil location based on means of largest contiguous area
				[vecLocY,vecLocX]=find(matEyeBW);
				vecBlobX(intBlob) = mean(vecLocX);
				vecBlobY(intBlob) = mean(vecLocY);
			end
			
			%get last location
			intLastDetect = find(~isnan(vecPosX),1,'last');
			dblLastPosX = vecPosX(intLastDetect);
			dblLastPosY = vecPosY(intLastDetect);
			vecBlobDist = sqrt((vecBlobX-dblLastPosX).^2+(vecBlobY-dblLastPosY).^2);
			vecPosW = 1-(vecBlobDist/dblMaxDist);
			
			%get size weight
			dblMu = mean(vecArea(vecArea>0));
			dblSigma = std(vecArea(vecArea>0));
			if isempty(dblSigma) || isnan(dblSigma) || dblSigma == 0
				vecSizeW = vecSizes/max(vecSizes);
			else
				dblMu = mean(vecArea(vecArea>0));
				dblSigma = std(vecArea(vecArea>0));
				vecP = normcdf(vecSizes,dblMu,dblSigma*5);
				vecSizeW = abs(abs(0.5-vecP)-0.5)*2;
			end
			
			%select blob
			[dblWeight,intIndex] = max(((vecSizeW/2+0.5).*(vecPosW/2+0.5).*vecBlobRoundness));
			dblRoundness = vecBlobRoundness(intIndex);
			
			%select area
			matEyeBW = false(size(matEyeBW));
			matEyeBW(sCC.PixelIdxList{intIndex}) = true;
			intSize = sum(matEyeBW(:));
			
			%find pupil location based on means
			[vecLocY,vecLocX]=find(matEyeBW);
			dblPosX = mean(vecLocX);
			dblPosY = mean(vecLocY);
			
			%update pupil luminance values
			dblPrevPupLum = dblThisPupLum;
		else
			%none found, assign dummy values
			dblWeight = 0;
			dblRoundness = 0;
			dblPosX = nan;
			dblPosY = nan;
			intSize = 0;
			dblThisPupLum = dblPrevPupLum;
		end
		
		%assign output
		vecWeight(intFrame) = dblWeight;
		vecRoundness(intFrame) = dblRoundness;
		vecPupilLuminance(intFrame) = dblNextPupLum;
		vecPosX(intFrame) = dblPosX;
		vecPosY(intFrame) = dblPosY;
		vecArea(intFrame) = intSize;
		matMask(:,:,intFrame) = matEyeBW;
		
		%save example output
		if mod(intFrame,1000) == 0
			%define location
			strFigName = sprintf('%sxyt%02d_F%05d.tif',strSes,intRec,intFrame);
			strTotName = sprintf('%s%s%s%s%s%s',strSourceDir,strSes,filesep,'exampleFits',filesep,strFigName);
			
			%make fig
			matRGB = imadjust(matMovie(:,:,intFrame));
			matRGB(:,:,2) = imadjust(matEye);
			matRGB(:,:,3) = double(matEyeBW);
			imshow(matRGB);
			title(strFigName);
			drawnow;
			
			%save
			imwrite(matRGB,strTotName,'tiff','Compression','lzw');
			
			%msg
			fprintf(' > Saved example fit to %s [%s]\n',strFigName,strTotName);
		end
		
		%check pupil position
		if intFrame == 1 && isnan(dblPosX)
			vecPosX(intFrame) = vecBaseLocation(1);
			vecPosY(intFrame) = vecBaseLocation(2);
		end
	end
	
	%detect blinks, saccades, etc
	vecEvents = abs(zscore(vecPupilLuminance).*zscore(vecArea));
	indEvents = vecEvents > 4;
	
	%put in output structure
	sEyeTracking = struct;
	sEyeTracking.strSes = strSes;
	sEyeTracking.strRec = strRec;
	sEyeTracking.intRec	= intRec;
	sEyeTracking.vecWeight = vecWeight;
	sEyeTracking.vecRoundness = vecRoundness;
	sEyeTracking.vecPupilLuminance = vecPupilLuminance;
	sEyeTracking.vecPosX = vecPosX;
	sEyeTracking.vecPosY = vecPosY;
	sEyeTracking.vecArea = vecArea;
	sEyeTracking.matMask = matMask;
	sEyeTracking.vecEvents = vecEvents;
	sEyeTracking.indEvents = indEvents;
	
	%save
	strFileEyeTracking = sprintf('EyeTrackData_%sxyt%02d',strSes,intRec);
	save(sprintf('%s%s%s%s',strSourceDir,strSes,filesep,strFileEyeTracking),'sEyeTracking','-v7.3');
	fprintf('\nSaved eye-tracking data of %sxyt%02d to %s [%s]\n',strSes,intRec,strSourceDir,getTime);
end