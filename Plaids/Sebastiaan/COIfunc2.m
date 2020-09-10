
function [CorrMatrixCOI] = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, dataset)
%% Variable preparation
%new variable for the activity of every neuron during presentation of the 10 different plaid stimuli

OrientVect = [2 4 5 7 9 11 14 15]; %these are all datasets that vary on the stimulus orientation, with a constant contrast
ContrastVect = [1 3 6 8 10 12 13 16]; %these are the datasets that vary on stimulus contrast, with a constant orientation
OrientPVVect = [4 5 7 14 15]; %These are the datasets that vary on stimulus orientation, with PV cell recordings
ContrastPVVect = [3 6 8 13 16]; %These are the datasets that vary on stimulus contrast, with PV cell recordings.
for i = 1:10
    plaidsref(i).ref = plaidsact.all(i).condition; 
end
%selecting the reference conditions.
ref_0_25 = plaidsref(1).ref;
ref_0_50 = plaidsref(2).ref;
ref_25_0 = plaidsref(6).ref;
ref_50_0 = plaidsref(7).ref;
%for every plaid dataset's activity matrix I select the rows and collums
%corresponding to the significant neurons, and the trials belonging to the
%different conditions.

for i = 1:length(orientall(dataset).significant)
    for k = 1:length(conditions(dataset).condition3)
        activity.single_25_25.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition3(k)));
        activity.single_25_50.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition4(k)));
        activity.single_50_25.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition8(k)));
        activity.single_50_50.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition5(k)));
    end
end

for i = 1:length(orientall(dataset).significant)
    for k = 1:length(conditions(dataset).condition1)
        activity.single_0_25.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition1(k)));
        activity.single_25_0.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition6(k)));
        activity.single_0_50.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition2(k)));
        activity.single_50_0.matr(i,k) = plaidsactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition7(k)));
    end
end

for i = 1:length(orientall(dataset).significant)
    for k = 1:length(conditions(dataset).condition10)
        activity.single_0_100.matr(i,k) = orientactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition10(k)));
        activity.single_100_0.matr(i,k) = orientactall(dataset).activ(orientall(dataset).significantindex(i),horzcat(conditions(dataset).condition9(k)));
    end
end
%% Calculating the COI for every neuron and trial within a condition, with the (A-R1-R2)/(R1+R2) formula
COImatrix.m_25_25 = COIcalc(activity.single_25_25.matr, ref_25_0, ref_0_25);
COImatrix.m_25_50 = COIcalc(activity.single_25_50.matr, ref_25_0, ref_0_50);
COImatrix.m_50_25 = COIcalc(activity.single_50_25.matr, ref_50_0, ref_0_25);
COImatrix.m_50_50 = COIcalc(activity.single_50_50.matr, ref_50_0, ref_0_50);
%calculating the mean COI for every neuron in the different conditions.
for j = 1:length(orientall(dataset).significant)
    COImatrix.f_25_25(j) = mean(COImatrix.m_25_25(j,:));
    COImatrix.f_25_50(j) = mean(COImatrix.m_25_50(j,:));
    COImatrix.f_50_25(j) = mean(COImatrix.m_50_25(j,:));
    COImatrix.f_50_50(j) = mean(COImatrix.m_50_50(j,:));
end
%% Selecting the relevant conditions out of every PV cell activity matrix
hx = find(ContrastVect == dataset);
[nRows, nCols] = size(plaidsPVactall(dataset).activ);
for i = 1:nRows
    for k = 1:length(conditions(dataset).condition3)
        PVactivity.all_25_25.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition3(k)));
        PVactivity.all_25_50.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition4(k)));
        PVactivity.all_50_25.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition8(k)));
        PVactivity.all_50_50.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition5(k)));
    end
end

for i = 1:nRows
    for k = 1:length(conditions(dataset).condition1)
        for h = 1:length(conditions(dataset).condition10)
            PVactivity.all_0_25.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition1(k)));
            PVactivity.all_25_0.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition6(k)));
            PVactivity.all_0_50.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition2(k)));
            PVactivity.all_50_0.matr(i,k) = plaidsPVactall(dataset).activ(i, horzcat(conditions(dataset).condition7(k)));
            PVactivity.all_0_100.matr(i,h) = orientpvactall(dataset).activ(i, horzcat(conditions(dataset).condition10(h)));
            PVactivity.all_100_0.matr(i,h) = orientpvactall(dataset).activ(i, horzcat(conditions(dataset).condition9(h)));
        end
    end
end


%taking the mean of the neurons for every trial, for every condition.
for j = 1:length(conditions(dataset).condition3)
    PVactivity.f_25_25(j) = mean(PVactivity.all_25_25.matr(:,j));
    PVactivity.f_25_50(j) = mean(PVactivity.all_25_50.matr(:,j));
    PVactivity.f_50_25(j) = mean(PVactivity.all_50_25.matr(:,j));
    PVactivity.f_50_50(j) = mean(PVactivity.all_50_50.matr(:,j));
end

for j = 1:length(conditions(1).condition1)
    PVactivity.f_0_25(j) = mean(PVactivity.all_0_25.matr(:,j));
    PVactivity.f_25_0(j) = mean(PVactivity.all_25_0.matr(:,j));
    PVactivity.f_0_50(j) = mean(PVactivity.all_0_50.matr(:,j));
    PVactivity.f_50_0(j) = mean(PVactivity.all_50_0.matr(:,j));
end

for j = 1:length(conditions(1).condition10)
    PVactivity.f_0_100(j) = mean(PVactivity.all_0_100.matr(:,j));
    PVactivity.f_100_0(j) = mean(PVactivity.all_100_0.matr(:,j));
end
%making a correlation matrix of the vector of activity of every PV cell's activity
%over the trials in a condition, and the vector of the cross orientation
%suppression per pyramidal neuron over the different trials in a condition
%and selecting only the non-diagonal correlation.

for i = 1:length(orientall(dataset).significant)

        [CorrMatrixCOI.vals.val_25_25(i), CorrMatrixCOI.vals.p_25_25(i)] = corr(PVactivity.f_25_25', COImatrix.m_25_25(i,:)');
        CorrMatrixCOI.final.COImean_25_25 = mean(CorrMatrixCOI.vals.val_25_25);%calculating the mean correlation.
        CorrMatrixCOI.final.Pmean_50_50 = mean(CorrMatrixCOI.vals.p_25_25);

        [CorrMatrixCOI.vals.val_25_50(i), CorrMatrixCOI.vals.p_25_50(i)] = corr(PVactivity.f_25_50', COImatrix.m_25_50(i,:)');
        CorrMatrixCOI.final.COImean_25_50 = mean(CorrMatrixCOI.vals.val_25_50);
        CorrMatrixCOI.final.Pmean_25_50 = mean(CorrMatrixCOI.vals.p_25_50);
        
        [CorrMatrixCOI.vals.val_50_25(i), CorrMatrixCOI.vals.p_50_25(i)] = corr(PVactivity.f_50_25', COImatrix.m_50_25(i,:)');
        CorrMatrixCOI.final.COImean_50_25 = mean(CorrMatrixCOI.vals.val_50_25);
        CorrMatrixCOI.final.Pmean_50_25 = mean(CorrMatrixCOI.vals.p_50_25);

        [CorrMatrixCOI.vals.val_50_50(i), CorrMatrixCOI.vals.p_50_50(i)] = corr(PVactivity.f_50_50', COImatrix.m_50_50(i,:)');
        CorrMatrixCOI.final.COImean_50_50 = mean(CorrMatrixCOI.vals.val_50_50);
        CorrMatrixCOI.final.Pmean_50_50 = mean(CorrMatrixCOI.vals.p_50_50);
        
        [CorrMatrixCOI.vals.val_0_25(i), CorrMatrixCOI.vals.p_0_25(i)] = corr(PVactivity.f_0_25', activity.single_0_25.matr(i,:)');
        CorrMatrixCOI.final.COImean_0_25 = mean(mean(CorrMatrixCOI.vals.val_0_25));%calculating the mean correlation.
        CorrMatrixCOI.final.Pmean_0_25 = mean(mean(CorrMatrixCOI.vals.p_0_25));

        [CorrMatrixCOI.vals.val_25_0(i), CorrMatrixCOI.vals.p_25_0(i)] = corr(PVactivity.f_25_0', activity.single_25_0.matr(i,:)');
        CorrMatrixCOI.final.COImean_25_0 = mean(mean(CorrMatrixCOI.vals.val_25_0));
        CorrMatrixCOI.final.Pmean_25_0 = mean(mean(CorrMatrixCOI.vals.p_25_0));
        
        [CorrMatrixCOI.vals.val_0_50(i), CorrMatrixCOI.vals.p_0_50(i)] = corr(PVactivity.f_0_50', activity.single_0_50.matr(i,:)');
        CorrMatrixCOI.final.COImean_0_25 = mean(mean(CorrMatrixCOI.vals.val_0_25));%calculating the mean correlation.
        CorrMatrixCOI.final.Pmean_0_25 = mean(mean(CorrMatrixCOI.vals.p_0_25));

        [CorrMatrixCOI.vals.val_50_0(i), CorrMatrixCOI.vals.p_50_0(i)] = corr(PVactivity.f_50_0', activity.single_50_0.matr(i,:)');
        CorrMatrixCOI.final.COImean_50_0 = mean(mean(CorrMatrixCOI.vals.val_50_0));
        CorrMatrixCOI.final.Pmean_50_0 = mean(mean(CorrMatrixCOI.vals.p_50_0));
        
        [CorrMatrixCOI.vals.val_0_100(i), CorrMatrixCOI.vals.p_0_100(i)] = corr(PVactivity.f_0_100', activity.single_0_100.matr(i,:)');
        CorrMatrixCOI.final.COImean_0_100 = mean(mean(CorrMatrixCOI.vals.val_0_100));%calculating the mean correlation.
        CorrMatrixCOI.final.Pmean_0_100 = mean(mean(CorrMatrixCOI.vals.p_0_100));

        [CorrMatrixCOI.vals.val_100_0(i), CorrMatrixCOI.vals.p_100_0(i)] = corr(PVactivity.f_100_0', activity.single_100_0.matr(i,:)');
        CorrMatrixCOI.final.COImean_100_0 = mean(mean(CorrMatrixCOI.vals.val_100_0));
        CorrMatrixCOI.final.Pmean_100_0 = mean(mean(CorrMatrixCOI.vals.p_100_0));
    end
end


