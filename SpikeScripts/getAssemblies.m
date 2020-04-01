function matAssemblies = getAssemblies(matSpikingHiRes,vecStarts,vecStops)
	%getAssemblies Retrieves assemblies from spiking matrix
	%   getAssemblies(matSpikingHiRes,vecStarts,vecStops)
	
	global matSpiking;
	matSpiking = matSpikingHiRes;
	matAssemblies = cell2mat(arrayfun(@getAssemblyNeurons,vecStarts,vecStops,'UniformOutput',false));
	clear global matSpiking;
end
function vecAssembly = getAssemblyNeurons(intStart,intStop)
	global matSpiking;
	vecAssembly = sum(matSpiking(:,intStart:intStop),2);
end
