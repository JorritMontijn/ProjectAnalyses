function vecOut = findvals(vecIn)
	%UNTITLED5 Summary of this function goes here
	%   Detailed explanation goes here
	
	vecOut = [];
	while ~isempty(vecIn)
		intMax = min(vecIn);
		vecOut(end+1) = intMax;
		vecIn = vecIn(vecIn~=intMax);
	end
end

