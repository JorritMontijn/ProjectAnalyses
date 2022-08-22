%step 1; convert movies
%framecapvec transforms avi to .mat [heel_video]

%step 2; detect pupil location
%requires user input at beginning
%plaatesmakenmetellipse gives as output x-y pixel locations per frame [_eyetracking.mat]

%step 3; blink detector
%blinkdet outputs frame_blinks.mat

%step 4; resampling
%getEyetrackingData resamples all data to calcium imaging frame
%acquisitions (dFoF time points)

%source/target data
strTargetPath = 'D:\Data\Processed\imagingvideo';
strSourcePath = 'G:\Data\Raw\imagingvideo';
intMouse = 1;
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
	vecRecordings = 1:8;
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
	
	
end