%every dataset is loaded onto the same structure for easier handling with
%loops
dat(1) = load("20130612xyt01_ses.mat");
dat(2) = load("20130612xyt02_ses.mat");
dat(3) = load("20130625xyt01_ses.mat");
dat(4) = load("20130625xyt02_ses.mat");
dat(5) = load("20131016xyt01_ses.mat");
dat(6) = load("20131016xyt02_ses.mat");
dat(7) = load("20131022xyt01_ses.mat");
dat(8) = load("20131022xyt02_ses.mat");
dat(9) = load("20140129xyt01_ses.mat");
dat(10) = load("20140129xyt02_ses.mat");
dat(11) = load("20140314xyt08_ses.mat");
dat(12) = load("20140314xyt09_ses.mat");
dat(13) = load("20140423xyt01_ses.mat");
dat(14) = load("20140423xyt02_ses.mat");
dat(15) = load("20140425xyt09_ses.mat");
dat(16) = load("20140425xyt10_ses.mat");

OrientVect = [2 4 5 7 9 11 14 15]; %these are all datasets that vary on the stimulus orientation, with a constant contrast
ContrastVect = [1 3 6 8 10 12 13 16]; %these are the datasets that vary on stimulus contrast, with a constant orientation
OrientPVVect = [4 5 7 14 15]; %These are the datasets that vary on stimulus orientation, with PV cell recordings
ContrastPVVect = [3 6 8 13 16]; %These are the datasets that vary on stimulus contrast, with PV cell recordings.\

%% Preallocation
orientall = struct('significantindex', cell(1,8), 'degrees', cell(1,8), 'significant', cell(1,8));
orientactall = struct('neuron', cell(1,8), 'structStim', cell(1,8), 'n_neuronen', cell(1,8), 'trials', cell(1,8), 'n_trials', cell(1,8), 'activ', cell(1,8), 'activZG', cell(1,8));
orientpvactall = struct('activ', cell(1,8), 'activZ', cell(1,8));
contrastactall = struct('neuron', cell(1,8), 'structStim', cell(1,8), 'UC', cell(1,8), 'n_neuronen', cell(1,8), 'n_contrasts', cell(1,8), 'trials', cell(1,8), 'n_trials', cell(1,8), 'stimindex', cell(1,8), 'neur', cell(1,8), 'neurvect', cell(1,8), 'max_val', cell(1,8), 'activ', cell(1,8), 'activZ', cell(1,8));
plaidsactall = struct('neuron', cell(1,8), 'structStim', cell(1,8), 'UC', cell(1,8), 'n_neuronen', cell(1,8), 'n_contrasts', cell(1,8), 'trials', cell(1,8), 'n_trials', cell(1,8), 'activ', cell(1,8), 'activZ', cell(1,8));
plaidsPVactall = struct('activ', cell(1,8), 'activZ', cell(1,8), 'structStim', cell(1,8));
plaid = struct('trial', cell(1,8), 'contrindex', cell(1,8));
plaids = struct('comborientation', cell(1,8), 'combcontrast', cell(1,8));
conditions = struct('condition1', cell(1,9), 'condition2', cell(1,9), 'condition3', cell(1,9), 'condition4', cell(1,9), 'condition5', cell(1,9), 'condition6', cell(1,9), 'condition7', cell(1,9), 'condition8', cell(1,9), 'condition9', cell(1,9), 'condition10', cell(1,9));
plaidconditions = struct('indices', cell(1,8));
Neuron = struct('tuning', cell(1,10));
plaidneurons = struct('tuning', cell(1,10)); 
%%
for j = 1:8
    orientall(j) = OrientTuningData(dat(OrientVect(j)), dat(ContrastVect(j))); %creates vectors of the preferred orientations of every neuron in every dataset, and the indices of the significant tuning curves.
    orientactall(j) = OrientactTuningData(dat(OrientVect(j)), dat(ContrastVect(j))); %creates the response matrices for the orientation datasets.
end



for j = 1:length(OrientPVVect)
    orientpvactall(ismember(OrientVect, OrientPVVect(j))).activ = OrientPVactTuningData(dat(OrientPVVect(j)), dat(ContrastPVVect(j))).activG;
    orientpvactall(ismember(OrientVect, OrientPVVect(j))).activZ = OrientPVactTuningData(dat(OrientPVVect(j)), dat(ContrastPVVect(j))).activZG;
end

for j = 1:8
    contrastactall(j) = ContrastactTuningData(dat(OrientVect(j)),dat(ContrastVect(j)));
end
for j = 1:8
    plaidsactall(j) = PlaidsactTuningData(dat(OrientVect(j)), dat(ContrastVect(j)));
end

xt = find(ismember(OrientVect, OrientPVVect))

for i = 1:5
    plaidsPVactall(xt(i)).activ = PlaidsPVactTuningData(dat(OrientPVVect(i)), dat(ContrastPVVect(i))).activ;
    plaidsPVactall(xt(i)).activZ = PlaidsPVactTuningData(dat(OrientPVVect(i)), dat(ContrastPVVect(i))).activZ;
end
clear xt
%still need to optimize this part.
orientPV1 = OrientPVTuningData(dat(4), dat(3), 2);
orientPV2 = OrientPVTuningData(dat(5), dat(6), 3);
orientPV3 = OrientPVTuningData(dat(7), dat(8), 4);
orientPV4 = OrientPVTuningData(dat(14), dat(13), 7);
orientPV5 = OrientPVTuningData(dat(15), dat(16), 8);
orientPVall(1).significant = orientPV1(2).significant;
orientPVall(2).significant = orientPV2(3).significant;
orientPVall(3).significant = orientPV3(4).significant;
orientPVall(4).significant = orientPV4(7).significant;
orientPVall(5).significant = orientPV5(8).significant;
orientPVall(1).degrees = orientPV1(2).degrees;
orientPVall(2).degrees = orientPV2(3).degrees;
orientPVall(3).degrees = orientPV3(4).degrees;
orientPVall(4).degrees = orientPV4(7).degrees;
orientPVall(5).degrees = orientPV5(8).degrees;
orientPVall(1).significantindex = orientPV1(2).significantindex;
orientPVall(2).significantindex = orientPV2(3).significantindex;
orientPVall(3).significantindex = orientPV3(4).significantindex;
orientPVall(4).significantindex = orientPV4(7).significantindex;
orientPVall(5).significantindex = orientPV5(8).significantindex;
clear orientPV1 orientPV2 orientPV3 orientPV4 orientPV5

%% Initial parameter creation
for j = 1:length(ContrastVect)
    plaid(j).trial = unique(dat(ContrastVect(j)).ses.structStim.TrialNumber);
    for i = 1:length(plaid(j).trial)
        plaid(j).contrindex(i).indices = find(dat(ContrastVect(j)).ses.structStim.TrialNumber == plaid(j).trial(i));
        plaid(j).contrindex(i).orientation1 = dat(ContrastVect(j)).ses.structStim.Orientation(plaid(j).contrindex(i).indices(1));
        if length(plaid(j).contrindex(i).indices) == 2
            plaid(j).contrindex(i).orientation2 = dat(ContrastVect(j)).ses.structStim.Orientation(plaid(j).contrindex(i).indices(2));
        end
        plaid(j).contrindex(i).contrast1 = dat(ContrastVect(j)).ses.structStim.Contrast(plaid(j).contrindex(i).indices(1));
        if length(plaid(j).contrindex(i).indices) == 2
            plaid(j).contrindex(i).contrast2 = dat(ContrastVect(j)).ses.structStim.Contrast(plaid(j).contrindex(i).indices(2));
        end
    end
    for i = 1:length(plaid(j).contrindex)
      if (plaid(j).contrindex(i).orientation1 == 90)
        plaid(j).contrindex(i).orientation2 = 0;
      end
      if (plaid(j).contrindex(i).orientation1 == 0)
        plaid(j).contrindex(i).orientation2 = 90;
      end
      if isempty(plaid(j).contrindex(i).contrast2)
        plaid(j).contrindex(i).contrast2 = 0;
      end
    end
end
for j = 1:length(ContrastVect)
        plaids(j).comborientation.comb90_0 = find(vertcat(plaid(j).contrindex.orientation1) == 90 & vertcat(plaid(j).contrindex.orientation2) == 0);
        plaids(j).comborientation.comb0_90 = find(vertcat(plaid(j).contrindex.orientation1) == 0 & vertcat(plaid(j).contrindex.orientation2) == 90);
        plaids(j).combcontrast.comb0_25 = find(vertcat(plaid(j).contrindex.contrast1) == 0 & vertcat(plaid(j).contrindex.contrast2) == 25);
        plaids(j).combcontrast.comb25_0 = find(vertcat(plaid(j).contrindex.contrast1) == 25 & vertcat(plaid(j).contrindex.contrast2) == 0);
        plaids(j).combcontrast.comb0_50 = find(vertcat(plaid(j).contrindex.contrast1) == 0 & vertcat(plaid(j).contrindex.contrast2) == 50);
        plaids(j).combcontrast.comb50_0 = find(vertcat(plaid(j).contrindex.contrast1) == 50 & vertcat(plaid(j).contrindex.contrast2) == 0);
        plaids(j).combcontrast.comb25_25 = find(vertcat(plaid(j).contrindex.contrast1) == 25 & vertcat(plaid(j).contrindex.contrast2) == 25);
        plaids(j).combcontrast.comb25_50 = find(vertcat(plaid(j).contrindex.contrast1) == 25 & vertcat(plaid(j).contrindex.contrast2) == 50);
        plaids(j).combcontrast.comb50_25 = find(vertcat(plaid(j).contrindex.contrast1) == 50 & vertcat(plaid(j).contrindex.contrast2) == 25);
        plaids(j).combcontrast.comb50_50 = find(vertcat(plaid(j).contrindex.contrast1) == 50 & vertcat(plaid(j).contrindex.contrast2) == 50);
     
        conditions(j).condition1 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 0) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 25)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 0))) ;
        conditions(j).condition2 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 0) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 50)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 0))) ;
        conditions(j).condition3 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 25)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 25))) ;
        conditions(j).condition4 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 50)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 25))) ;
        conditions(j).condition5 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 50)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 50))) ;
        conditions(j).condition6 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 0)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 0) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 25))) ;
        conditions(j).condition7 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 0)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 0) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 50))) ;
        conditions(j).condition8 = find(((vertcat((plaid(j).contrindex.orientation1)) == 0) & (vertcat((plaid(j).contrindex.contrast1)) == 50) & (vertcat((plaid(j).contrindex.orientation2)) == 90) & (vertcat((plaid(j).contrindex.contrast2)) == 25)) | ((vertcat((plaid(j).contrindex.orientation1)) == 90) & (vertcat((plaid(j).contrindex.contrast1)) == 25) & (vertcat((plaid(j).contrindex.orientation2)) == 0) & (vertcat((plaid(j).contrindex.contrast2)) == 50))) ;    
end
for j = 1:length(ContrastVect)
    for i = 1:12
        plaidconditions(j).indices(i).vector1 = plaid(1).contrindex([conditions(j).condition1(i)]).indices(1);
        plaidconditions(j).indices(i).vector2 = plaid(1).contrindex([conditions(j).condition2(i)]).indices(1);
        plaidconditions(j).indices(i).vector6 = plaid(1).contrindex([conditions(j).condition6(i)]).indices(1);
        plaidconditions(j).indices(i).vector7 = plaid(1).contrindex([conditions(j).condition7(i)]).indices(1);
    end
end
for j = 1:length(ContrastVect)
    for i = 1:24
        plaidconditions(j).indices(i).vector3 = plaid(1).contrindex([conditions(j).condition3(i)]).indices(1);
        plaidconditions(j).indices(i).vector4 = plaid(1).contrindex([conditions(j).condition4(i)]).indices(1);
        plaidconditions(j).indices(i).vector5 = plaid(1).contrindex([conditions(j).condition5(i)]).indices(1);
        plaidconditions(j).indices(i).vector8 = plaid(1).contrindex([conditions(j).condition8(i)]).indices(1);
    end
end
for j = 1:8
    for i = 1:length(plaid(j).contrindex)
        if length(plaid(j).contrindex(i).indices) == 2
            if dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(1)) > dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(2))
                plaid(j).contrindex(i).FrameOn = dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(1));
                plaid(j).contrindex(i).TrialfixOn = plaid(j).contrindex(i).indices(1);
            elseif dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(2)) >= dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(1))
                plaid(j).contrindex(i).FrameOn = dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices(2));
                plaid(j).contrindex(i).TrialfixOn = plaid(j).contrindex(i).indices(2);
            end
        elseif length(plaid(j).contrindex(i).indices) == 1
            plaid(j).contrindex(i).FrameOn = dat(ContrastVect(j)).ses.structStim.FrameOn(plaid(j).contrindex(i).indices);
            plaid(j).contrindex(i).TrialfixOn = plaid(j).contrindex(i).indices;
        end
            if length(plaid(j).contrindex(i).indices) == 2
            if dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(1)) < dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(2))
                plaid(j).contrindex(i).FrameOff = dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(1));
                plaid(j).contrindex(i).TrialfixOff = plaid(j).contrindex(i).indices(1);
            elseif dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(2)) <= dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(1))
                plaid(j).contrindex(i).FrameOff = dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices(2));
                plaid(j).contrindex(i).TrialfixOff = plaid(j).contrindex(i).indices(2);
            end
        elseif length(plaid(j).contrindex(i).indices) == 1
            plaid(j).contrindex(i).FrameOff = dat(ContrastVect(j)).ses.structStim.FrameOff(plaid(j).contrindex(i).indices);
            plaid(j).contrindex(i).TrialfixOff = plaid(j).contrindex(i).indices;
        end
        
    end
end
for j = 1:8
    for i = 1:length(orientall(j).significant)
        plaidsact.vect(j).condition1(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition1)));
        plaidsact.vect(j).condition2(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition2)));
        plaidsact.vect(j).condition3(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition3)));
        plaidsact.vect(j).condition4(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition4)));
        plaidsact.vect(j).condition5(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition5)));
        plaidsact.vect(j).condition6(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition6)));
        plaidsact.vect(j).condition7(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition7)));
        plaidsact.vect(j).condition8(i)           =                  mean(plaidsactall(j).activ(orientall(j).significantindex(i),horzcat(conditions(j).condition8)));
    end
end
for j = 1:length(OrientVect)
    conditions(j).condition9 = find(dat(OrientVect(j)).ses.structStim.Orientation == 0);
    conditions(j).condition10 = find(dat(OrientVect(j)).ses.structStim.Orientation == 90);    
end
for j = 1:8
    for i = 1:length(orientall(j).significant)
        plaidsact.vect(j).condition9(i) = mean(orientactall(j).activ(orientall(j).significantindex(i),transpose(conditions(j).condition9)));
        plaidsact.vect(j).condition10(i) = mean(orientactall(j).activ(orientall(j).significantindex(i),transpose(conditions(j).condition10)));
    end
end
plaidsact.all(1).condition = horzcat(plaidsact.vect(1:8).condition1);
plaidsact.all(2).condition = horzcat(plaidsact.vect(1:8).condition2);
plaidsact.all(3).condition = horzcat(plaidsact.vect(1:8).condition3);
plaidsact.all(4).condition = horzcat(plaidsact.vect(1:8).condition4);
plaidsact.all(5).condition = horzcat(plaidsact.vect(1:8).condition5);
plaidsact.all(6).condition = horzcat(plaidsact.vect(1:8).condition6);
plaidsact.all(7).condition = horzcat(plaidsact.vect(1:8).condition7);
plaidsact.all(8).condition = horzcat(plaidsact.vect(1:8).condition8);
plaidsact.all(9).condition = horzcat(plaidsact.vect(1:8).condition9);
plaidsact.all(10).condition = horzcat(plaidsact.vect(1:8).condition10);
for i = 1:10
    Neuron(i).tuning = transpose(vertcat(orientall(1:8).significant));
end
n_tunings = 15;
tuningspace= linspace(0, 180, (n_tunings+1));
for j = 1:10
    for i = 1:n_tunings
        plaidneurons(j).tuning(i).indices = find(Neuron(j).tuning >= tuningspace(i) & Neuron(j).tuning < tuningspace(i+1));
        plaidneurons(j).tuning(i).AvrgOri = (tuningspace(i) + tuningspace (i+1))/2;
        plaidneurons(j).tuning(i).activity = mean(plaidsact.all(j).condition(plaidneurons(j).tuning(i).indices));
    end
end

for j =1:10
    for i=1:n_tunings
        plaidsact.vectstd.condition(j).tuning(i) = std(plaidsact.all(j).condition(plaidneurons(j).tuning(i).indices));
        plaidsact.vectse.condition(j).tuning(i) = plaidsact.vectstd.condition(j).tuning(i)/sqrt(length(plaidneurons(j).tuning(i).indices));
    end
end

conditions(9).condition1 = "0&25";
conditions(9).condition2 = "0&50";
conditions(9).condition3 = "25&25";
conditions(9).condition4 = "25&50";
conditions(9).condition5 = "50&50";
conditions(9).condition6 = "25&0";
conditions(9).condition7 = "50&0";
conditions(9).condition8 = "50&25";
conditions(9).condition9 = "100&0";
conditions(9).condition10 = "0&100";

for i = 1:8
    condition.condition1(i,:) = conditions(i).condition1;
    condition.condition2(i,:) = conditions(i).condition2;
    condition.condition3(i,:) = conditions(i).condition3;
    condition.condition4(i,:) = conditions(i).condition4;
    condition.condition5(i,:) = conditions(i).condition5;
    condition.condition6(i,:) = conditions(i).condition6;
    condition.condition7(i,:) = conditions(i).condition7;
    condition.condition8(i,:) = conditions(i).condition8;
    condition.condition9(i,:) = conditions(i).condition9;
    condition.condition10(i,:) = conditions(i).condition10;
end
