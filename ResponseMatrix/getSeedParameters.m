function vecParameters = getSeedParameters(matResponse,cellSelect,vecOrientations)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	%pre-allocate output vector
	vecParameters = nan(1,8);
	%1 = MeanPopResp
	%2 = StDevPopResp
	%3 = MeanPopBandwidth
	%4 = StDevPopBandwidth
	%5 = NeuronMeanStDevResp
	%6 = NeuronStDevStDevResp
	%7 = NeuronMeanFano
	%8 = NeuronStDevFano
	
	%get response params
	sTuning = calcTuningRespMat(matResponse,cellSelect,vecOrientations);
	
	%get relevant parameters
	vecMeanPrefResp = sTuning.vecMeanPrefResp;
	vecSDPrefResp = sTuning.vecSDPrefResp;
	vecFano = sTuning.vecFano;
	
	%calc mean/sd's
	dblMeanPopResp = mean(vecMeanPrefResp);
	dblStDevPopResp = std(vecMeanPrefResp);
	dblNeuronMeanStDevResp = mean(vecSDPrefResp);
	dblNeuronStDevStDevResp = std(vecSDPrefResp);
	dblNeuronMeanFano = mean(vecFano);
	dblNeuronStDevFano = std(vecFano);
	
	%get bandwidth
	dblMeanPopBandwidth = 10;
	dblStDevPopBandwidth = 2;
	
	%assign to output
	vecParameters(1) = dblMeanPopResp;
	vecParameters(2) = dblStDevPopResp;
	vecParameters(3) = dblMeanPopBandwidth;
	vecParameters(4) = dblStDevPopBandwidth;
	vecParameters(5) = dblNeuronMeanStDevResp;
	vecParameters(6) = dblNeuronStDevStDevResp;
	vecParameters(7) = dblNeuronMeanFano;
	vecParameters(8) = dblNeuronStDevFano;
end