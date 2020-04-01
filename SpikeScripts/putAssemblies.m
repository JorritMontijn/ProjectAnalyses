function matAssemblyActivity = putAssemblies(matAssemblyActivity,vecA,vecStarts,vecStops)
	%putAssemblies Transforms vector-based data to matrix data
	%   matAssemblyActivity = putAssemblies(matAssemblyActivity,vecA,vecStarts,vecStops)
	
	for intA=1:length(vecA)
		matAssemblyActivity(vecA(intA),vecStarts(intA):vecStops(intA)) = true;
	end
end
