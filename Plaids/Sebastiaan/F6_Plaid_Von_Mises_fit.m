
%% contrast 0-0 code
for i = 1:8
    datco(ContrastVect(i)).ses.neuron = dat(ContrastVect(i)).ses.neuron(orientall(i).significantindex);
    datco(OrientVect(i)).ses.neuron = dat(OrientVect(i)).ses.neuron(orientall(i).significantindex);
end

activity90.peak = max(vertcat(plaidneurons(10).tuning.activity));
activity0.peak = max(vertcat(plaidneurons(9).tuning.activity));

for j = 1:8
    for i = 1:length(datco(ContrastVect(j)).ses.neuron)
        for k = 1:length(plaid(j).contrindex)
            for h = 1:length(dat(OrientVect(j)).ses.structStim.FrameOn)
            baseline.means(ContrastVect(j)).neuron(i).timeframe(k) = mean(datco(ContrastVect(j)).ses.neuron(i).dFoF((plaid(j).contrindex(k).FrameOn-25):plaid(j).contrindex(k).FrameOn));
            baseline.means(OrientVect(j)).neuron(i).timeframe(h) = mean(datco(OrientVect(j)).ses.neuron(i).dFoF((dat(OrientVect(j)).ses.structStim.FrameOn(h)-25):dat(OrientVect(j)).ses.structStim.FrameOn(h)));
            end
        end
    end
end

for j = 1:8
    for i = 1:length(datco(ContrastVect(j)).ses.neuron)
        baseline.combined(j).neuron(i).timeframe = horzcat(baseline.means(ContrastVect(j)).neuron(i).timeframe, baseline.means(OrientVect(j)).neuron(i).timeframe);
    end
end

for j = 1:8
    for i = 1:length(baseline.combined(j).neuron)
        baseline.combined(j).neuron(i).gm = mean(baseline.combined(j).neuron(i).timeframe);
        baselinegms(j).neuron = vertcat(baseline.combined(j).neuron.gm);
    end
end

COmeans = transpose(vertcat(baselinegms(1:8).neuron));

%% CONTRAST 0-0 PLOTS
n_tunings = 15;
tuningspace= linspace(0, 180, (n_tunings+1));
for i = 1:n_tunings
    COneurons.combined(i).indices = find(Neuron(1).tuning >= tuningspace(i) & Neuron(1).tuning < tuningspace(i+1));
    COneurons.combined(i).AvrgOri = (tuningspace(i) + tuningspace (i+1))/2;
    COneurons.combined(i).activity = mean(COmeans(COneurons.combined(i).indices));
end

    for i=1:n_tunings
        CO.combinedstd.tuning(i) = std(COmeans(COneurons.combined(i).indices));
        CO.combinedse.tuning(i) = CO.combinedstd.tuning(i)/sqrt(length(COneurons.combined(i).indices));
    end

%% Plot functions
PAB_1 = figure;
subplot(4,4,1)
PlaidsVMplot(COneurons.combined, CO.combinedse.tuning, " ");

subplot(4,4,2)
PlaidsVMplot(plaidneurons(1).tuning, plaidsact.vectse.condition(1).tuning, " ");

subplot(4,4,3)
PlaidsVMplot(plaidneurons(2).tuning, plaidsact.vectse.condition(2).tuning, " ");

subplot(4,4,4)
PlaidsVMplot(plaidneurons(10).tuning, plaidsact.vectse.condition(10).tuning, " ");

subplot(4,4,5)
PlaidsVMplot(plaidneurons(6).tuning, plaidsact.vectse.condition(6).tuning, " ");

subplot(4,4,6)
PlaidsVMplot(plaidneurons(3).tuning, plaidsact.vectse.condition(3).tuning, " ");

subplot(4,4,7)
PlaidsVMplot(plaidneurons(4).tuning, plaidsact.vectse.condition(4).tuning, " ");

subplot(4,4,9)
PlaidsVMplot(plaidneurons(7).tuning, plaidsact.vectse.condition(7).tuning, " ");

subplot(4,4,10)
PlaidsVMplot(plaidneurons(8).tuning, plaidsact.vectse.condition(8).tuning, " ");

subplot(4,4,11)
PlaidsVMplot(plaidneurons(5).tuning, plaidsact.vectse.condition(5).tuning, " ");

subplot(4,4,13)
PlaidsVMplot(plaidneurons(9).tuning, plaidsact.vectse.condition(9).tuning, " ");

print(gcf, '-bestfit','PBA_1','-dpdf')
%%

