lh.Color = [0.7,0.7,0.7];
lh.Width = 0.5;
lh.AverageWidth = 1.5;

y = deg2rad(1:360);
for i = [1, 2, 4, 5, 6, 7, 8] %3 does not have significant neurons.
    sOut(i) = getTuningCurves(orientactall(i).activZG, dat(OrientVect(i)).ses.structStim.Orientation);
    plot_VM(i).significantparams = sOut(i).matFittedParams(sOut(i).vecOriTtest<=0.05,:);
    
    for j = 1:length(plot_VM(i).significantparams)
        plot_VM(i).neuron(j).neuronfit = vonMisesSingleFitPX(plot_VM(i).significantparams(j,:),y); %werkt niet meer. Kan het kloppen dat vonMisesSingleFitPX is aangepast? anders moet ik even kijken wat er anders kan zijn.
    end
    for j = 1:length(plot_VM(i).significantparams)
        for k = 1:length(y)
            plot_VM(i).neuron(j).neuronfitpct(k) = (plot_VM(i).neuron(j).neuronfit(k)/max(plot_VM(i).neuron(j).neuronfit))*100;
            [ymax ind] = max(plot_VM(i).neuron(j).neuronfitpct);
            plot_VM(i).neuron(j).shiftneuronfitpct = circshift(plot_VM(i).neuron(j).neuronfitpct, 90-(ind));
        end
    end
    
end

for j = 1:8
    for i = 1:length(plot_VM(j).neuron)
        allneurons(j).matshiftpct(i,:) = plot_VM(j).neuron(i).shiftneuronfitpct;
        allneurons(j).matshift(i,:) = plot_VM(j).neuron(i).neuronfitpct;
    end
end
allneurons_gm = mean(vertcat(allneurons(1:8).matshiftpct));

