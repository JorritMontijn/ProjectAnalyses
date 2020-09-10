function [CorrMatrixCOI] = COIfunc(dat, plaidsPVactall, plaidsact, plaidsactall, orientall, conditions, dataset)
%% Variable preparation
OrientVect = [2 4 5 7 9 11 14 15]; %these are all datasets that vary on the stimulus orientation, with a constant contrast
ContrastVect = [1 3 6 8 10 12 13 16]; %these are the datasets that vary on stimulus contrast, with a constant orientation
OrientPVVect = [4 5 7 14 15]; %These are the datasets that vary on stimulus orientation, with PV cell recordings
ContrastPVVect = [3 6 8 13 16]; %These are the datasets that vary on stimulus contrast, with PV cell recordings.


%for every plaid dataset's activity matrix I select the rows and collums
%corresponding to the significant neurons, and the trials belonging to the
%different conditions.

for i = 1:length(orientall(dataset).significant)
    for k = 1:length(conditions(dataset).condition1)
        activity.single_0_25.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition1(k)));
        activity.single_25_0.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition6(k)));
    end
end

%% Selecting the relevant conditions out of every PV cell activity matrix
hx = find(ContrastVect == dataset);
for i = 1:length(dat(ContrastVect(dataset)).ses.PV)
    for k = 1:length(conditions(dataset).condition3)
        PVactivity.all_0_25.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition1(k)));
        PVactivity.all_25_0.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition6(k)));
    end
end
%concatenating the matrices

%taking the mean of the trials for every neuron, for every condition.
for j = 1:length(dat(ContrastVect(dataset)).ses.PV)
    PVactivity.f_0_25(j) = mean(PVactivity.all_0_25.matr(j,:));
    PVactivity.f_25_0(j) = mean(PVactivity.all_25_0.matr(j,:));
end
%making a correlation matrix of the vector of activity of every PV cell's activity
%over the trials in a condition, and the vector of the cross orientation
%suppression per pyramidal neuron over the different trials in a condition
%and selecting only the non-diagonal correlation.

for i = 1:length(orientall(dataset).significant)
    for j = 1:length(dat(ContrastVect(dataset)).ses.PV)
        [CorrMatrixCOI.vals.val_0_25(j,i), CorrMatrixCOI.vals.p_0_25(j,i)] = corr(PVactivity.all_0_25.matr(j,:)', activity.single_0_25.matr(i,:)');
        CorrMatrixCOI.final.COImean_0_25 = mean(mean(CorrMatrixCOI.vals.val_0_25));%calculating the mean correlation.
        CorrMatrixCOI.final.Pmean_0_25 = mean(mean(CorrMatrixCOI.vals.p_0_25));

        [CorrMatrixCOI.vals.val_25_0(j,i), CorrMatrixCOI.vals.p_25_0(j,i)] = corr(PVactivity.all_25_0.matr(j,:)', COImatrix.m_25_0(i,:)');
        CorrMatrixCOI.final.COImean_25_0 = mean(mean(CorrMatrixCOI.vals.val_25_0));
        CorrMatrixCOI.final.Pmean_25_0 = mean(mean(CorrMatrixCOI.vals.p_25_0));
        
        
        
    end
end

