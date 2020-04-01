% This script detects different stimulus types from ses.structStim.
% Output: matrix 'trialMat'
% 			[   trial number ;
% 				trial type ;
% 				orientation of first presented stimuli, or, in situation 7+8, the orientation of stimulus with the highest contrast ;
% 				orientation of mask ;
% 				frame number where trial begins ;
% 				frame number where trial ends ;
%				duration of the trial in seconds ;
%				number of 3s bins in trial]
% Input: ses

function trialMat = stimDetection(ses)
	
	trialMat = nan(8, max(ses.structStim.TrialNumber));
	
	for i = 1 : length(ses.structStim.TrialNumber)-1
		
		% Situation 1 - 90 deg test, 0 deg mask
		
		if ses.structStim.Contrast(i) == 25
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Orientation(i+1) == 0
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 1;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 1 - 0 deg test, 90 deg mask
		
		if ses.structStim.Contrast(i) == 25
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Orientation(i+1) == 90
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 1;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 2 - 0 deg
		if ses.structStim.Orientation(i) == 0
			if ses.structStim.Orientation(i+1) == -1
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) <=100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 2;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 2;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
					
				end
				
			end
		end
		
		% Situation 2 - 90 deg
		if ses.structStim.Orientation(i) == 90
			if ses.structStim.Orientation(i+1) == -1
				
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) <= 100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 2;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 2;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
					
					
				end
				
			end
		end
		
		% Situation 3 - 0 deg
		if ses.structStim.Orientation(i) == 0
			if ses.structStim.Orientation(i+1) == -1
				if ses.structStim.Contrast(i) == 25
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) >= 100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 3;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 3;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
				
				
				
			end
		end
		
		
		% Situation 3 - 90 deg
		if ses.structStim.Orientation(i) == 90
			if ses.structStim.Orientation(i+1) == -1
				
				if ses.structStim.Contrast(i) == 25
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) >= 100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 3;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 3;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
						
					end
					
				end
				
			end
		end
		
		% Situation 4 - 0 deg test, 90 deg mask
		
		if ses.structStim.Contrast(i) == 50
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i+1) == 25
					if ses.structStim.Orientation(i+1) == 90
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 4;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 4 - 90 deg test, 0 deg mask
		
		if ses.structStim.Contrast(i) == 50
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i+1) == 25
					if ses.structStim.Orientation(i+1) == 0
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 4;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 5 - 0 deg
		if ses.structStim.Orientation(i) == 0
			if ses.structStim.Orientation(i+1) == -1
				
				
				if ses.structStim.Contrast(i) == 25
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) <=100
						if i == 1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 5;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
							
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 5;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
				
				
			end
		end
		
		% Situation 5 - 90 deg
		if ses.structStim.Orientation(i) == 90
			if ses.structStim.Orientation(i+1) == -1
				
				if ses.structStim.Contrast(i) == 25
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) <=100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 5;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 5;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
				
			end
		end
		
		% Situation 6 - 0 deg
		if ses.structStim.Orientation(i) == 0
			if ses.structStim.Orientation(i+1) == -1
				
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) >=100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 6;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 6;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
				
			end
		end
		
		% Situation 6 - 90 deg
		if ses.structStim.Orientation(i) == 90
			if ses.structStim.Orientation(i+1) == -1
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.FrameOff(i)-ses.structStim.FrameOn(i) >=100
						if i==1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 6;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						elseif ses.structStim.Orientation(i-1) == -1
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 6;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = NaN;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
			
		end
		
		% Situation 7 - 90 deg is high
		if ses.structStim.FrameOn(i) == ses.structStim.FrameOn(i+1)
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.Contrast(i+1) == 25
						if ses.structStim.FrameOff(i) <= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 7;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i+1);
						end
					end
					
				end
			elseif ses.structStim.Orientation(i+1) == 90
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Contrast(i) == 25
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 7;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
				
			end
		end
		
		
		% Situation 7 - 0 deg is high
		if ses.structStim.FrameOn(i) == ses.structStim.FrameOn(i+1)
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.Contrast(i+1) == 25
						if ses.structStim.FrameOff(i) <= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 7;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i+1);
						end
					end
				end
			elseif ses.structStim.Orientation(i+1) == 0
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Contrast(i) == 25
						if ses.structStim.FrameOff(i+1) <= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 7;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
					
				end
			end
		end
		
		% Situation 8 - 90 deg is high
		if ses.structStim.FrameOn(i) == ses.structStim.FrameOn(i+1)
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.Contrast(i+1) == 25
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 8;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
			elseif ses.structStim.Orientation(i+1) == 90
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Contrast(i) == 25
						if ses.structStim.FrameOff(i+1) >= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 8;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i+1);
						end
					end
				end
			end
		end
		
		% Situation 8 - 0 deg is high
		if ses.structStim.FrameOn(i) == ses.structStim.FrameOn(i+1)
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i) == 50
					if ses.structStim.Contrast(i+1) == 25
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 8;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
						end
					end
				end
			elseif ses.structStim.Orientation(i+1) == 0
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Contrast(i) == 25
						if ses.structStim.FrameOff(i+1) >= ses.structStim.FrameOff(i)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 8;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i+1);
						end
					end
				end
			end
		end
		
		% Situation 9 - 90 deg test, 0 deg mask
		
		if ses.structStim.Contrast(i) == 25
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i+1) == 25
					if ses.structStim.Orientation(i+1) == 0
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 9;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 9 - 0 deg test, 90 deg mask
		
		if ses.structStim.Contrast(i) == 25
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i+1) == 25
					if ses.structStim.Orientation(i+1) == 90
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 9;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 10 - 90 deg test, 0 deg mask
		
		if ses.structStim.Contrast(i) == 50
			if ses.structStim.Orientation(i) == 90
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Orientation(i+1) == 0
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 10;
							trialMat(3, trialNumber) = 90;
							trialMat(4, trialNumber) = 0;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		% Situation 9 - 0 deg test, 90 deg mask
		
		if ses.structStim.Contrast(i) == 50
			if ses.structStim.Orientation(i) == 0
				if ses.structStim.Contrast(i+1) == 50
					if ses.structStim.Orientation(i+1) == 90
						if ses.structStim.FrameOff(i) >= ses.structStim.FrameOff(i+1)
							trialNumber = ses.structStim.TrialNumber(i);
							
							trialMat(1, trialNumber) = trialNumber;
							trialMat(2, trialNumber) = 10;
							trialMat(3, trialNumber) = 0;
							trialMat(4, trialNumber) = 90;
							trialMat(5, trialNumber) = ses.structStim.FrameOn(i);
							trialMat(6, trialNumber) = ses.structStim.FrameOff(i);
							
						end
					end
				end
			end
		end
		
		%		trialMat(7,:) = ((trialMat(6,:)-trialMat(5,:))/ses.samplingFreq);
		
		for j = 1 : length(trialMat)
			trialMat(7,j) = round(((trialMat(6,j)-trialMat(5,j))/ses.samplingFreq));
			trialMat(8,j) = trialMat(7,j)/3;
		end
		
	end