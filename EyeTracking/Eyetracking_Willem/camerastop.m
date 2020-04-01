function fStop = camerastop(video_heel)
	
	%calculate video stopping
	for c = length(video_heel.frames):-1:1
		
		if c>1
			
			compare_1 = mean(mean(mean(im2double(video_heel.frames(1,c-1).cdata))));
			compare_2 = mean(mean(mean(im2double(video_heel.frames(1,c).cdata))));	
			
			if ((compare_2 - compare_1)/compare_1) < -0.5
				fStop = c;
				%break
			end
			
		else fStop = length(video_heel.frames);
			
		end
		
	end