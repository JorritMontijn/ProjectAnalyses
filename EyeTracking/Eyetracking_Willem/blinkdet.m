strDir = 'D:\Data\Processed\video\';

%Load eyetrack file
if ~exist('eyetracking','var')
	eye_strFile = '20140430xyt01_vid_eyetracking.mat';
	load([strDir eye_strFile]);
end
%Load video file
if ~exist('video_heel','var')
	video_strFile = '20140430xyt01_video_heel.mat';
	load([strDir video_strFile]);
end

blink_strFile = [eye_strFile(1:(length(eye_strFile)-12)),'_blinkdata.mat'];

%clear b_saved
close all

%Initialize
ecclipse_area = eyetracking.ecclipse_area;
dblCenterX = eyetracking.dblCenterX;
dblCenterY = eyetracking.dblCenterY;
frame = eyetracking.frame;
eccback = ecclipse_area;
ecclipse_a=ecclipse_area;
ecclipse_a(ecclipse_area==0)=[];
frame=frame(1:length(ecclipse_a));
ecclipse_area=ecclipse_area(frame);
dblCenterX=dblCenterX(frame);
dblCenterY=dblCenterY(frame);

% 
% %DF over F of ecclipse area
% eccl_norm = calcdFoF(ecclipse_area,100,4);
% 
% nidx = find(eccl_norm==(mean(eccl_norm)-std(eccl_norm)));
% pidx = find(eccl_norm==(mean(eccl_norm)+std(eccl_norm)));
% if isempty(nidx)
% 	nidx = find(diff(sign(eccl_norm-(mean(eccl_norm)-std(eccl_norm)))));
% end
% if isempty(pidx)
% 	pidx = find(diff(sign(eccl_norm-(mean(eccl_norm)+std(eccl_norm)))));
% end
% 
% idx = unique([nidx,pidx]);
% 
% figure
% plot(frame,eccl_norm)
% xlabel('Frame #')
% ylabel('delta Ecclipse Area over Ecclipse Area')
% hold on
% plot(frame(idx),eccl_norm(idx),'r*');
% hold on
% plot(frame,mean(eccl_norm))
% plot(frame,mean(eccl_norm)-std(eccl_norm))
% plot(frame,mean(eccl_norm)+std(eccl_norm))
% threshold = 5;
% array = idx;
% sortedArray = sort(idx);
% nPerGroup = diff(find([1 (diff(sortedArray) > threshold) 1]));
% groupArray = mat2cell(sortedArray,1,nPerGroup);
% blink = zeros(1,length(groupArray));
% contin = 0;


% if isempty(mnidx)
% 	mnidx = find(diff(sign(ecclipse_area-(mean(ecclipse_area(ecclipse_area>2))-std(ecclipse_area(ecclipse_area>2))))));
% end
% if isempty(mpidx)
% 	mpidx = find(diff(sign(ecclipse_area-(mean(ecclipse_area(ecclipse_area>2))+std(ecclipse_area(ecclipse_area>2))))));
% end

%vanaf hier
mnidx = find(ecclipse_area<(mean(ecclipse_area(ecclipse_area>2))-2*std(ecclipse_area(ecclipse_area>2))));
mpidx = find(ecclipse_area>(mean(ecclipse_area(ecclipse_area>2))+2*std(ecclipse_area(ecclipse_area>2))));

midx = unique([mnidx,mpidx]);

figure
plot(frame,ecclipse_area)
xlabel('Frame #')
ylabel('Ecclipse Area')
hold on
plot(frame(midx),ecclipse_area(midx),'r*');
hold on
plot(frame,mean(ecclipse_area(ecclipse_area>2)))
plot(frame,mean(ecclipse_area(ecclipse_area>2))-2*std(ecclipse_area(ecclipse_area>2)))
plot(frame,mean(ecclipse_area(ecclipse_area>2))+2*std(ecclipse_area(ecclipse_area>2)))

threshold = 4;
array = midx; %idx
sortedArray = sort(midx); %idx
% array = idx; %idx
% sortedArray = sort(idx); %idx
nPerGroup = diff(find([1 (diff(sortedArray) > threshold) 1]));
groupArray = mat2cell(sortedArray,1,nPerGroup);
blink = zeros(1,length(groupArray));
contin = 0;

current = pwd;
cd('D:\Data\Processed\video\');
savefile = [eye_strFile(1:(length(eye_strFile)-12)),'_frame_blinks.mat'];
save(savefile, 'groupArray');
cd(current);

figure

uicontrol('Style', 'pushbutton', 'String', 'Blink',...
	'Position', [220 10 60 20],...
	'Callback', 'blink(b) = 1; contin = 1;');

uicontrol('Style', 'pushbutton', 'String', 'No Blink',...
	'Position', [290 10 60 20],...
	'Callback', 'contin = 1;');

uicontrol('Style', 'pushbutton', 'String', 'Save and Exit',...
	'Position', [360 10 90 20],...
	'Callback', 'contin = 2;');

uicontrol('Style', 'pushbutton', 'String', 'Load Session',...
	'Position', [120 10 90 20],...
	'Callback', 'contin = 3;');

frame_start=zeros(length(groupArray),1);
frame_end=zeros(length(groupArray),1);

for b = 1:length(groupArray)
	
	%Check if session is loaded, if so then resume at last blink
	if exist('b_saved') == 1
		if b < b_saved
			continue
		end
	end
	
	contin = 0;
	
	if groupArray{b}(1)<2
		fra = frame(groupArray{b}(1)):frame(groupArray{b}(end)+1);
		frame_end(b,1) = frame(groupArray{b}(end)+1);
		frame_start(b,1) = frame(groupArray{b}(1));
	elseif b == length(groupArray)
		fra = frame(groupArray{b}(1)-1):frame(groupArray{b}(end));
		frame_end(b,1) = frame(groupArray{b}(end));
		frame_start(b,1) = frame(groupArray{b}(1)-1);
	else
		fra = frame(groupArray{b}(1)-1):frame(groupArray{b}(end)+1);
		frame_end(b,1) = frame(groupArray{b}(end)+1);
		frame_start(b,1) = frame(groupArray{b}(1)-1);
	end
	
	if b>1
		if ishandle(texti)==1
			delete(texti)
		end
	end
	
	img1=cell(1,fra(end));
	
	for f = fra
		img1{f} = im2double(video_heel.frames(1,f).cdata);
	end
	
	while contin == 0
		
		for f = fra
			
			tic
			t = video_heel.times(1,f);
			
			if f>fra(1)
				delete(texti)
			else
			end
			
			image1 = image(img1{f});
			hold on
			
			texti = text(3*pi/4,-1,['time = ',int2str(t) ' frame = ' int2str(f) ' blink ' int2str(b) ' of ' int2str(length(groupArray))],'FontSize',16);
			drawnow
			
			if b < length(groupArray)
				if f == frame(groupArray{b}(end)+1)
					delete(texti)
				end
			end
			
			%pause(0.5)
			toc
			clearvars -except img1 eyetracking eccback video_heel contin groupArray b blink ecclipse_area frame t texti f fra dblCenterY dblCenterX eyeposxovertime eyeposyovertime frame_end frame_start strDir blink_strFile eye_strFile
			
			%Get out of for loop when button is pressed.
			if contin > 0
				break
			end
			
		end
		
	end
	
	%Save and stop
	if contin == 2
		break
	end
	
	if contin == 3
		
		clearvars -except video_heel strDir eye_strFile eyetracking frame texti contin blink_strFile
		load([strDir blink_strFile]);
		frame_start = blinkdata{:,2};
		frame_end = blinkdata{:,3};
		blink = blinkdata{:,4}';
		groupArray = blinkdata{:,5}';
		b_saved = find(frame_start>0,1,'last');
		contin = 0;
		
	end
	
end

%Create data table

blinkdata{:,1}= [1:length(groupArray)]';
blinkdata{:,2} = frame_start;
blinkdata{:,3} = frame_end;
blinkdata{:,4} = [blink]';
blinkdata{:,5} = groupArray';
current = pwd;
cd('D:\Data\Processed\video\');
savefile = [eye_strFile(1:(length(eye_strFile)-12)),'_blinkdata.mat'];
save(savefile, 'blinkdata');
cd(current);
close all

