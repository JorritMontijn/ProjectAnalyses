%Directory of map session file
strDir = 'D:\Data\Processed\imagingdata\20140430\';
%strDir = 'D:\Data\Processed\imagingdata\20140314\';
%strDir = 'C:\Users\willem\Desktop\laatste opdracht\';
%Load session
if ~exist('ses','var')
	%ses_strFile = '20140314xyt02_ses.mat';
	ses_strFile = '20140430xyt01_ses.mat';
	load([strDir ses_strFile]);
end

strDir = 'D:\Data\Processed\video\';
%strDir = 'C:\Users\willem\Desktop\laatste opdracht\';

%Load eyetrack file
if ~exist('eyetracking','var')
	%eye_strFile = '20140314xyt02_eyetracking.mat';
	eye_strFile = '20140430xyt01_vid_eyetracking.mat';
	load([strDir eye_strFile]);
end

%Load blinkdata
if ~exist('blinkdata','var')
	%video_strFile = '20140314xyt02_blinkdata.mat';
	video_strFile = '20140430xyt01_vid_eye_blinkdata.mat';
	load([strDir video_strFile]);
end

%Load videoheel
if ~exist('video_heel','var')
	%video_strFile = '20140314xyt02_video_heel.mat';
	video_strFile = '20140430xyt01_video_heel.mat';
	load([strDir video_strFile]);
end

%Make eyetracking vector composed out of zeros and ones for blink epoch.
blinkvectorRAW = zeros(1,length(eyetracking.ecclipse_area));
blinkvector = zeros(1,length(eyetracking.ecclipse_area));
index = find(blinkdata{:,4} == 1);
rawindex = [blinkdata{:,5}{:,:}];

for i=1:length(index)
	
	blinkindextemp{i} = blinkdata{:,5}{index(i),:};
	
end

blinkindex = cell2mat(blinkindextemp);
if sum(blinkindex)>0
	blinkvector(blinkindex) = 1;
else
end

if sum(rawindex)>0
	blinkvectorRAW(rawindex) = 1;
else
end

%calculate infrared total video time.
videotime = video_heel.times(end)-video_heel.times(1);

%calculate frame on and offset of calcium imaging in the video recordings
run eyetrackdetectstart
%F1 = eyetracking.frame(1);
MicTimeOn = video_heel.times(f1);
%and offset
Fend = camerastop(video_heel);
MicTimeOff = video_heel.times(Fend);

if Fend > length(eyetracking.ecclipse_area)
	
	eyetracking.ecclipse_area(Fend) = 0;
	eyetracking.dblCenterX(Fend) = 0;
	eyetracking.dblCenterY(Fend) = 0;
	blinkvector(Fend) = 0;
	blinkvectorRAW(Fend) = 0;
	
	%align eyetracking data to microscope recordings vector
	blinkvector=blinkvector(f1:Fend);
	blinkvectorRAW=blinkvectorRAW(f1:Fend);
	ecclipse_area = eyetracking.ecclipse_area(f1:Fend);
	dblCenterX = eyetracking.dblCenterX(f1:Fend) ;
	dblCenterY = eyetracking.dblCenterY(f1:Fend) ;
	
else
	
	blinkvector=blinkvector(f1:Fend);
	blinkvectorRAW=blinkvectorRAW(f1:Fend);
	ecclipse_area = eyetracking.ecclipse_area(f1:Fend);
	dblCenterX = eyetracking.dblCenterX(f1:Fend) ;
	dblCenterY = eyetracking.dblCenterY(f1:Fend) ;
	
end

%find starting value in
startVal = find(dblCenterX,1,'first');

%fill vectors with repeats
for fR = startVal:3:length(blinkvector)
	
	if dblCenterX(fR)> 0
		dblCenterX(fR+1)=dblCenterX(fR);
		dblCenterX(fR+2)=dblCenterX(fR);
	else
	end
	
	if dblCenterY(fR)> 0
		dblCenterY(fR+1)=dblCenterY(fR);
		dblCenterY(fR+2)=dblCenterY(fR);
	else
	end
	
	if ecclipse_area(fR)> 0
		ecclipse_area(fR+1)=ecclipse_area(fR);
		ecclipse_area(fR+2)=ecclipse_area(fR);
	else
	end
	
end

%frames to realtime seconds
rt = f1/29.9699970030003;
%sampling rate microscope
micHz = ses.samplingFreq;
%microscope start frame
micsF = rt*ses.samplingFreq;
%Total duration
micDur = ses.time.dur;
micFrames = micHz*micDur;

%resample data
PQ = length(ses.neuron(1, 1).vecSpikes)/length(blinkvector);
[P Q] = rat(PQ);
resampled = resample(blinkvector,P,Q);
resampledRAW = resample(blinkvectorRAW,P,Q);

RS = length(ses.neuron(1, 1).vecSpikes)/length(dblCenterX);
[R S] = rat(RS);
Rs_dblCenterX=resample(dblCenterX,R,S);
Rs_dblCenterY=resample(dblCenterY,R,S);
Rs_ecclipse_area=resample(ecclipse_area,R,S);
% [P Q] = rat(length(ses.neuron(1, 1).vecSpikes)/length(resampled));
% resampled2 = resample(resampled,P,Q);
% [P Q] = rat(length(ses.neuron(1, 1).vecSpikes)/length(resampledRAW));
% resampled2RAW = resample(resampledRAW,P,Q);

BlinksResampled=resampled;
BlinksResampledRaw=resampledRAW;
BlinksResampled(resampled>=0.5)=1;
BlinksResampled(resampled<0.5)=0;
BlinksResampledRaw(resampledRAW>=0.5)=1;
BlinksResampledRaw(resampledRAW<0.5)=0;

%Put everything in Output ses
ses.eyetracking.EyeCenter_X_Coordinate = Rs_dblCenterX;
ses.eyetracking.EyeCenter_Y_Coordinate = Rs_dblCenterY;
ses.eyetracking.Eye_areasize = Rs_ecclipse_area;
ses.eyetracking.All_Blinks = BlinksResampledRaw;
ses.eyetracking.Manually_selected_Blinks = BlinksResampled;

%clearvars -except ses video_heel blinkdata eyetracking
