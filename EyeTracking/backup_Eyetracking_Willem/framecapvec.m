function framecapvec(filename)

%Find all specified video files in mmread map to convert
fileList = dir('F:\WillemBruin\Processing\Toolboxes\mmread\20140129*');
nFiles = length(fileList);
	
%	close all
	
	%framecapvec imports video files at the maximum possible amount of frames which can be
	%proccesed at a time. A region of interrest is to be selected first.
	%Input requires 'filename' of the video file located in the mmread folder.
	%Output saves total video processed in .ses matlab file.
	%
	%Also requires to change current folder to this function!
	%
	%syntax: video_parts = framecapvec(filename)
	%	input:
	%	- 'filename': video file located in mmread folder
	%
	%	output:
	%	'filename_video_heel.mat'
	%
	%Dependencies:
	% - mmread
	% - exampleROIselect.jpg located in mmread folder.
	
	%Shows an example of how to select ROI.
	exampleROIselect = imread('exampleROIselect.jpg');
	figure;
	imshow(exampleROIselect)
	
	pilot = mmread(filename,1:2);
	%Runs mmread for first 2 frames of specified file to get total amount of
	%frames.
	a = abs(pilot.nrFramesTotal);
	%Divides total amount of frames by 1000 to get proccesable parts and saves
	%the determined frames in a vector.
	b = ceil((abs(pilot.nrFramesTotal))/1000);
	framecap = linspace(1,a,b);
	framecap = round(framecap);
	
	
	for i=1:(length(framecap)-1)
		tic
		
		video_parts = mmread(filename,framecap(i):framecap(i+1));
		
		for f = 1:length(video_parts.frames)
			
			%For first frame ROI is to be specified.
			if f == 1 && i == 1
				
				figure;
				hi_pick = image(video_parts.frames(1,round((length(video_parts.frames))/2)).cdata);
				colormap((video_parts.frames(1,(round(length(video_parts.frames)/2))).colormap));
				img1 = get(hi_pick,'CData');
				
				h_im = imshow(img1);
				
				e = imrect;
				BW = createMask(e,h_im);
				
				BW = uint8(BW);
				
				pos = getPosition(e);
				
				video_parts.frames(1,f).cdata = imcrop(video_parts.frames(1,f).cdata,pos);
				
			else
				
				video_parts.frames(1,f).cdata = imcrop(video_parts.frames(1,f).cdata,pos);
				
			end
			
		end
		
		if i < 10
			k = num2str(0);
		else k = '';
		end
		
		close all
		
		savefile = [ filename '_video_part_' k num2str(i) '.mat' ];
		save(savefile, 'video_parts');
		disp([num2str(i), ' of ', num2str(length(framecap)), ' processed'])
		clearvars -except f i framecap filename pos
		
		toc
		
	end
	
	%Find all video part files and load them in 1 variable loadedData
	
	fileList = dir([filename '_video_part*']);
	nFiles = length(fileList);
	
	%# loop through files and write results into loadedData
	
	for iFile = 1:nFiles
		
		tmp = load(fileList(iFile).name);
		loadedData(iFile) = tmp;
		
	end
	
	%Merge the structs in loaddata
	
	fields = fieldnames(loadedData(1).video_parts)';
	
	B = cat(3,loadedData.video_parts);
	out = reshape(B(1,1,:),1,[]);
	fields(2,:) = cellfun(@(f) [out(1,:).(f)], fields, 'unif', false);
	
	%fields(2,:) = cellfun(@(f) [loadedData(1,1).video_parts.(f) loadedData(1,2).video_parts.(f) loadedData(1,3).video_parts.(f) loadedData(1,4).video_parts.(f) loadedData(1,5).video_parts.(f) loadedData(1,6).video_parts.(f) loadedData(1,7).video_parts.(f) loadedData(1,8).video_parts.(f) loadedData(1,9).video_parts.(f) loadedData(1,10).video_parts.(f) loadedData(1,11).video_parts.(f) loadedData(1,12).video_parts.(f) loadedData(1,13).video_parts.(f) loadedData(1,14).video_parts.(f) loadedData(1,15).video_parts.(f) loadedData(1,16).video_parts.(f) loadedData(1,17).video_parts.(f) loadedData(1,18).video_parts.(f) loadedData(1,19).video_parts.(f) loadedData(1,20).video_parts.(f) loadedData(1,21).video_parts.(f) loadedData(1,22).video_parts.(f) loadedData(1,23).video_parts.(f) loadedData(1,24).video_parts.(f) loadedData(1,25).video_parts.(f) loadedData(1,26).video_parts.(f) loadedData(1,27).video_parts.(f)], fields, 'unif', false);
	%fields(2,:) = cellfun(@(f) cat(3,loadedData.(f).video_parts), fields, 'unif', false);
	
	%save file
	
	filename = regexprep(filename,'[^a-zA-Z-0-9]','');
	filename = filename(1:(length(filename)-3));
	video_heel = struct(fields{:});
	current = pwd;
	cd('D:\Data\Processed\video\');
	savefile = [ filename '_video_heel.mat' ];
	save(savefile, 'video_heel');
	cd(current);
	disp('Done!');
	