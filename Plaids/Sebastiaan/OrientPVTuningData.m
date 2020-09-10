function [orientPV] = OrientPVTuningData(orientationdata,contrastdata, datalocation)
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

orientPVact(datalocation).PV   = sSesAggregate.PV;
orientPVact(datalocation).structStim   = sSesAggregate.structStim;

orientPVact(datalocation).h = unique(orientPVact(datalocation).structStim.Orientation); %vindt alle bestaande unieke oriëntaties uit de data
orientPVact(datalocation).n_PV = length(vertcat(orientPVact(datalocation).PV.id));
orientPVact(datalocation).n_orientaties = length(orientPVact(datalocation).h);
orientPVact(datalocation).trials = orientPVact(datalocation).structStim.Orientation;
orientPVact(datalocation).n_trials = length(orientPVact(datalocation).trials);


%preallocation for variables


orientPVact(datalocation).stimindex = struct("indices", cell(1, orientPVact(datalocation).n_orientaties), 'orientation', cell(1,orientPVact(datalocation).n_orientaties));
orientPVact(datalocation).PVcell = struct('meansvect', cell(1, orientPVact(datalocation).n_orientaties));
orientPVact(datalocation).PVcellvect = zeros(1,orientPVact(datalocation).n_PV);
orientPVact(datalocation).max_val = zeros(1,orientPVact(datalocation).n_trials);
orientPVact(datalocation).activ = zeros(orientPVact(datalocation).n_PV, orientPVact(datalocation).n_trials);
orientPVact(datalocation).activity = struct('tuning', cell(1, orientPVact(datalocation).n_orientaties));

%finding the different stimulus orientations and their respective indexes.

    for i = 1:orientPVact(datalocation).n_orientaties
        orientPVact(datalocation).stimindex(i).indices=find(orientPVact(datalocation).structStim.Orientation==orientPVact(datalocation).h(i));
        orientPVact(datalocation).stimindex(i).orientation=orientPVact(datalocation).h(i);
    end


%determining the mean response of every neuron to the different stimulus
%presentations, corresponding to the orientation of that stimulus
    for k = 1:orientPVact(datalocation).n_PV
        for i = 1:orientPVact(datalocation).n_orientaties
            for j = 1:length(orientPVact(datalocation).stimindex(i).indices)
                orientPVact(datalocation).stimindex(i).PVcell(k).meanresponse(j) = mean(orientPVact(datalocation).PV(k).dFoF(orientPVact(datalocation).structStim.FrameOn(orientPVact(datalocation).stimindex(i).indices(j)):orientPVact(datalocation).structStim.FrameOff(orientPVact(datalocation).stimindex(i).indices(j))));
            end
        end
    end
  % voor elke oriëntatie worden voor alle neuronen de gemiddelde activiteiten tijdens stimuluspresentatie bepaald en als .meanresponse opgeslagen.


    for k = 1:orientPVact(datalocation).n_PV
        for i = 1:orientPVact(datalocation).n_orientaties
            orientPVact(datalocation).stimindex(i).PVcell(k).grandmean = mean(orientPVact(datalocation).stimindex(i).PVcell(k).meanresponse);%voor alle oriëntaties wordt voor alle neuronen het gemiddelde berekend van de gemiddelde activiteiten tijdens stimuluspresentatie.
            orientPVact(datalocation).PVcell(k).meansvect(i) = orientPVact(datalocation).stimindex(i).PVcell(k).grandmean; %maakt een nieuwe struct aan voor de gemiddeldes. Maakt vorige stap arbitrair?
            orientPVact(datalocation).PVcellvect(i,k) = orientPVact(datalocation).PVcell(k).meansvect(i); %maakt een matrix met de gemiddeldes, waarbij de rows de orientaties zijn en de collumns de neuronen
        end
    end



    for i = 1:orientPVact(datalocation).n_PV
        for k = 1:orientPVact(datalocation).n_trials
            orientPVact(datalocation).activ(i,k) = mean(orientPVact(datalocation).PV(i).dFoF(orientPVact(datalocation).structStim.FrameOn((k)):orientPVact(datalocation).structStim.FrameOff((k))));
        end
    end
orientPVact(datalocation).activZ = zscore(orientPVact(datalocation).activ,[],1); %z-scored at each trial to remove mean population activity



indGratings = [241:320];
orientPVact(datalocation).activZG = orientPVact(datalocation).activZ(:,indGratings);
orientPVact(datalocation).activG = orientPVact(datalocation).activ(:,indGratings);


[sOut] = getTuningCurves(orientPVact(datalocation).activZG, orientationdata.ses.structStim.Orientation);

    orientPV(datalocation).significantindex = find(sOut.vecOriTtest <=1);
    orientPV(datalocation).degrees = (rad2deg(sOut.matFittedParams(:,1)))/2;
    orientPV(datalocation).significant = orientPV(datalocation).degrees(orientPV(datalocation).significantindex);


end


