function matM = getBlockMeans(matIn,vecI)
	%getBlockMeans Returns mean values in [n x n] 2D matrix per identity
	%block given by [n x 1] vecI; between-block means are not returned
	%   matM = getBlockMeans(matIn,vecI)
	
	%globals
	global matV;
	global vecID;
	global vecStarts;
	global vecStops;
	
	%reorder by blocks
	if size(vecI,2) == 1,vecI = vecI';end
	[vecID,vecReorder] = sort(vecI,'ascend');
	matV = matIn(vecReorder,vecReorder);
	vecStarts = find([1 vecID(2:end)-vecID(1:(end-1))]);
	vecStops = find([vecID(2:end)-vecID(1:(end-1)) 1]);
	
	%get matrix
	matM = cell2mat(arrayfun(@getSingleBlock,unique(vecID),'UniformOutput',false));
	
	%clear globals
	clear global matV;
	clear global vecID;
	clear global vecStarts;
	clear global vecStops;
end
function vecM = getSingleBlock(intID)
	%globals
	global matV;
	global vecID;
	global vecStarts;
	global vecStops;

	
	matSingleBlock = matV(:,vecID==intID);
	vecM = nan(length(unique(vecID)),1);
	for intID2=unique(vecID)
		intStart = vecStarts(intID2);
		intStop = vecStops(intID2);
		vecM(intID2) = nanmean(nanmean(matSingleBlock(intStart:intStop,:)));
	end
end
