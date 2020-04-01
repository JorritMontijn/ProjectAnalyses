%step 1; convert movies
%requires user input at beginning
%runEyeDetectPrePro crops and reformats movies

%step 2; detect pupil location
%requires user input at beginning
%runPupilDetection requests onset/offset boundaries and then detects the pupil

%step 3; resampling
%getResampledEyetrackingData resamples all data to calcium imaging frame

%source/target data
%clear all;
strTempPath = 'D:\Data\Processed\imagingvideo\temp';
strTargetPath = 'D:\Data\Processed\imagingvideo';
strSourcePath = 'G:\Data\Raw\imagingvideo';


%pre-allocate ROI selection structure
sROI = [];
for intMouse=4
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
		vecRecordings = 3;
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
	
	%transform avi to mat
	
	%make target dir
	if ~isdir([strTargetPath filesep strSes])
		mkdir(strTargetPath,strSes);
	end
	
	%loop through recordings
	for intRecording=vecRecordings
		%get source file
		sFile = dir([strSourcePath filesep strSes filesep sprintf('*xyt%02d.avi',intRecording)]);
		strSourceFile = [strSourcePath filesep strSes filesep sFile(1).name];
		
		%get ROI
		intPos = numel(sROI);
		if intPos > 0
			vecPrevRect = sROI(intPos).vecRect;
		else
			vecPrevRect = [];
		end
		vecRect = getMovieROI(strSourceFile,vecPrevRect);
		
		%save ROI data
		sROI(intPos+1).strSes = strSes;
		sROI(intPos+1).strRec = sprintf('xyt%02d',intRecording);
		sROI(intPos+1).intRec = intRecording;
		sROI(intPos+1).strSourceFile = strSourceFile;
		sROI(intPos+1).vecRect = vecRect;
	end
	
	%save ROI
	strFileROI = ['MovieROIs_' date];
	save(sprintf('%s%s%s',strTargetPath,filesep,strFileROI),'sROI');
	fprintf('\nSaved ROIs to %s%s%s.mat\n',strTargetPath,filesep,strFileROI);
end


%loop through recordings and process all movies
intMaxFrameLoad = 1000;
for intFile=1:numel(sROI)
	strMovieFile = sROI(intFile).strSourceFile;
	strSes = sROI(intFile).strSes;
	strRec = sROI(intFile).strRec;
	vecRect = round(sROI(intFile).vecRect);
	
	%get initial video data
	sVideoInit = mmread(strMovieFile,1:2);
	intFramesTot = abs(sVideoInit.nrFramesTotal);
	
	%calc how many parts to split into, given certain max frame load
	intTotSteps = ceil((abs(sVideoInit.nrFramesTotal))/intMaxFrameLoad);
	vecMaxFrame = linspace(1,intFramesTot,intTotSteps);
	vecMaxFrame = round(vecMaxFrame);
	
	%loop through video
	intTotParts = (length(vecMaxFrame)-1);
	for intPart=1:intTotParts
		%get video
		sVideo = mmread(strMovieFile,vecMaxFrame(intPart):vecMaxFrame(intPart+1));
		sVideo.matMovie = zeros(vecRect(4)+1,vecRect(3)+1,numel(sVideoInit.frames),'uint8');
		
		%crop video
		for intFrame = 1:length(sVideo.frames)
			sVideo.matMovie(:,:,intFrame) = mean(imcrop(sVideo.frames(intFrame).cdata,vecRect),3);
			sVideo.frames(intFrame).cdata = [];
		end
		
		%save data to temp file
		strSaveTempFile = [strTempPath filesep 'TMP_videopart_' sprintf('%03d',intPart) '.mat'];
		save(strSaveTempFile, 'sVideo');
		fprintf('Cropped part %d of %d [%s] for file %d/%d (%s%s%s)\n',intPart,intTotParts,getTime,intFile,numel(sROI),strSes,'\',strRec);
	end
	
	%concatenate everything
	%get files
	sDir = dir(sprintf('%s%sTMP_videopart_*.mat',strTempPath,filesep));
	sVideoAll = rmfield(sVideo,'frames');
	clear sVideo;
	sVideoAll.matMovie = [];
	
	%loop through files and put results in sVideoAll
	for intPart = 1:intTotParts
		%get file, load it & append to end of aggregate
		strFile = [strTempPath filesep sDir(intPart).name];
		sTemp = load(strFile);
		intTempFrames = size(sTemp.sVideo.matMovie,3);
		sVideoAll.matMovie(:,:,(end+1):(end+intTempFrames)) = sTemp.sVideo.matMovie;
		
		%delete temp file
		delete(strFile);
		
		%msg
		fprintf('Concatenated part %d of %d [%s] for file %d/%d (%s%s%s)\n',intPart,intTotParts,getTime,intFile,numel(sROI),strSes,'\',strRec);
	end
	
	%remove empty frames
	vecM=squeeze(mean(mean(sVideoAll.matMovie,1),2));
	sVideoAll.matMovie = sVideoAll.matMovie(:,:,vecM~=0);
	
	%save aggregate cropped video data
	strOutputFile = sprintf('%s_%s%s','croppedvideo',strSes,strRec);
	save([strTargetPath filesep strSes filesep strOutputFile], 'sVideoAll','-v7.3');
	
	%msg
	fprintf('\nFinished processing %s%s%s! [%s]\n\n',strSes,'\',strRec,getTime)
end