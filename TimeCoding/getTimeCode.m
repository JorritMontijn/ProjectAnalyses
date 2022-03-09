function dblTimeCode = getTimeCode(vecSpikeT,dblStartT,intSpikeNum,dblUseMaxDur)
	%getTimeCode Temporal code that converges to a rate code when all spikes are used over the whole period
	%   dblTimeCode = getTimeCode(vecSpikeT,dblStartT,intSpikeNum,dblUseMaxDur)
	
	vecISI = diff(cat(1,0,flat(vecSpikeT(vecSpikeT>dblStartT)-dblStartT),dblUseMaxDur));
	dblTimeCode = mean(vecISI(1:min(intSpikeNum,numel(vecISI))));%dblRate = (1/dblTimeCode) - 1/dblUseMaxDur; produces real spike rate
end

