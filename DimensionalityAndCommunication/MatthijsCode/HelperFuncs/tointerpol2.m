function [inMA] = tointerpol2(Ma,freq,freqw)


tmpt = size(Ma,2);
inMA=zeros(length(freqw),tmpt);

for i=1:tmpt
	inMA(:,i) = interp1(freq,Ma(:,i),freqw);
end

