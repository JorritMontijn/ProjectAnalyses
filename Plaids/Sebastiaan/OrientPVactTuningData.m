function [orientPVact] = OrientPVactTuningData(orientationdata,contrastdata)
%% concatenate
%add contrast field to gratings
sesG.structStim.Contrast = 100*ones(size(orientationdata.ses.structStim.Orientation));
sSesAggregate = buildMultiSesAggregate(orientationdata.ses, contrastdata.ses);

for intCellType=1:numel(sSesAggregate.cellType)
	strCellType = sSesAggregate.cellType{intCellType};
	if isfield(orientationdata.ses,strCellType)
		indKeepG = ~cellfun(@strcmpi,cellfill('absent',[1 numel(orientationdata.ses.(strCellType))]),{orientationdata.ses.(strCellType).strPresence});
		indKeepP = ~cellfun(@strcmpi,cellfill('absent',[1 numel(contrastdata.ses.(strCellType))]),{contrastdata.ses.(strCellType).strPresence});
		
		%check which to keep
		vecKeepAgg = zeros([1 numel(sSesAggregate.(strCellType))]);
		vecKeepAgg(indKeepG) = 1;
		vecKeepAgg(indKeepP) = vecKeepAgg(indKeepP) + 1;
		indKeepAgg = vecKeepAgg==2;
		%remove entries
		sSesAggregate.(strCellType) = sSesAggregate.(strCellType)(indKeepAgg);
	end
end

orientPVact.PV   = sSesAggregate.PV;
orientPVact.structStim   = sSesAggregate.structStim;
orientPVact.n_neuronen = length(vertcat(orientPVact.PV.id));
orientPVact.trials = orientPVact.structStim.Orientation;
orientPVact.n_trials = length(orientPVact.trials);

%preallocation for variables

orientPVact.activ = zeros(orientPVact.n_neuronen, orientPVact.n_trials);


for i = 1:orientPVact.n_neuronen
    for k = 1:orientPVact.n_trials
        orientPVact.activ(i,k) = mean(orientPVact.PV(i).dFoF(orientPVact.structStim.FrameOn((k)):orientPVact.structStim.FrameOff((k))));
    end
end


indGratings = [241:320];

orientPVact.activG = orientPVact.activ(:,indGratings);
orientPVact.activZG = zscore(orientPVact.activG, [], 2);

end


