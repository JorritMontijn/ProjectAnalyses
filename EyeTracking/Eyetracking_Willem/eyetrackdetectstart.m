for c = 1:length(video_heel.frames)
	
	if c<length(video_heel.frames)
		compare_1 = mean(mean(mean(im2double(video_heel.frames(1,c).cdata))));
		compare_2 = mean(mean(mean(im2double(video_heel.frames(1,c+1).cdata))));
		
		if ((compare_2 - compare_1)/compare_1) > 1.5
			f1 = c+1;
			%break
		end
		
	end
end
