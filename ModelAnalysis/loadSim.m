
if intLoadSim == -1 && boolLoad
	strSimulation = 'xAreaExperiment_106r001p26-t_2018-02-21';
elseif intLoadSim == -2 && boolLoad
	strSimulation = 'xAreaExperiment_107l003p143-t_2018-02-21';
	
elseif intLoadSim == -3
	%monyetp31: 124 neurons (recorded on 144 channels) in awake V1 in a monkey doing passive fixation (1 deg window). 2 gratings with 5 deg offset in orientation. Animal performed 831 trials with 3 stim presentations per trial, 200 ms each  with 150 ms ISI. Effectively ~2500 stimulus presentations divided across the two stimuli.
	%vecOris = [0 5];
	matCompareTypes = [1 2; 2 1];
	strSimulation = 'xAreaExperiment_monyetp31_alex_2018-06-14';
elseif intLoadSim == -4
	%monyetp33: 106 neurons in awake V1 with 4 gratings (0,5 and 90,95). 845 trials or~600 presentations of each stimulus.
	%vecOris = [0 5 90 95];
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_monyetp33_alex_2018-06-14';

elseif intLoadSim == -6
	%vecOris = [-5 -5 -5 0 0 0 5 5 5 85 85 85 90 90 90 95 95 95];
	%vecSFs = repmat([1 2 0.5],[1 6]);
	matCompareTypes = [1 4; 4 1; 4 7; 7 4;...
		10 13; 13 10; 13 16; 16 13];
	%matCompareTypes = [matCompareTypes;matCompareTypes+1;matCompareTypes+2];
	strSimulation = 'xAreaExperiment_139l001p113_alex_2018-06-20';
elseif intLoadSim == -7
	%vecOris = [-5 -5 -5 0 0 0 5 5 5 85 85 85 90 90 90 95 95 95];
	%vecSFs = repmat([1 2 0.5],[1 6]);
	matCompareTypes = [1 4; 4 1; 4 7; 7 4;...
		10 13; 13 10; 13 16; 16 13];
	%matCompareTypes = [matCompareTypes;matCompareTypes+1;matCompareTypes+2];
	%matCompareTypes = matCompareTypes(1:2,:);
	strSimulation = 'xAreaExperiment_139l001p115_alex_2018-06-20';
	
elseif intLoadSim == -8 && boolLoad
	matCompareTypes = [3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp13_2019-02-20';
elseif intLoadSim == -9 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp70_2019-02-20';
elseif intLoadSim == -10 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp149_2019-02-20';
elseif intLoadSim == -11 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp185_2019-02-20';
elseif intLoadSim == -12 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp428_2019-02-20';
elseif intLoadSim == -13 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_cadetp429_2019-02-20';
elseif intLoadSim == -14 && boolLoad
	matCompareTypes = [1 2; 2 1; 3 4; 4 3];
	strSimulation = 'xAreaExperiment_alcaponep42_alex_2019-10-08';

		

elseif intLoadSim == -21 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuS008_2018-07-19';
elseif intLoadSim == -22 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuS006_2018-07-19';
elseif intLoadSim == -23 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuS010_2018-07-19';
elseif intLoadSim == -24 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuP011_2018-07-19';
elseif intLoadSim == -25 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuP012_2018-07-19';
elseif intLoadSim == -26 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuP016_2018-07-19';
elseif intLoadSim == -27 && boolLoad
	strSimulation = 'xAreaExperiment_JoGuP014_2018-07-19';
	
elseif intLoadSim == 00 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise0RP_2018-04-30';
elseif intLoadSim == 01 && boolLoad
	strSimulation = 'xAreaDistributed_IndInpOri5Noise063_2018-11-26';
elseif intLoadSim == 02 && boolLoad
	%strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise063_2018-06-28';
	strSimulation = 'xAreaDistributed_IndInpOri5Noise088_2018-11-26';
elseif intLoadSim == 03 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise1RP_2018-04-30';
elseif intLoadSim == 04 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise5RP_2018-04-30';
elseif intLoadSim == 05 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise5_2018-04-20';


elseif intLoadSim == 11 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise0RP_2018-04-30'; %to do
elseif intLoadSim == 12 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise063_2018-06-28';
elseif intLoadSim == 13 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise1RP_2018-04-30'; %to do
elseif intLoadSim == 14 && boolLoad
	strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise5RP_2018-04-30'; %to do
	


elseif intLoadSim == 100 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise00_2020-10-08';
elseif intLoadSim == 102 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise02_2020-09-29';
elseif intLoadSim == 104 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise04_2020-09-29';
elseif intLoadSim == 106 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise06_2020-09-29';
elseif intLoadSim == 108 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise08_2020-09-29';
elseif intLoadSim == 110 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise10_2020-09-29';
	
	elseif intLoadSim == 112 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise12_2020-09-29';
	elseif intLoadSim == 114 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise14_2020-09-29';
	elseif intLoadSim == 116 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise16_2020-09-29';
	elseif intLoadSim == 118 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise18_2020-09-29';
	elseif intLoadSim == 120 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise20_2020-09-29';
	
	elseif intLoadSim == 122 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise22_2020-09-29';
	elseif intLoadSim == 124 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise24_2020-09-29';
	elseif intLoadSim == 126 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise26_2020-09-29';
	elseif intLoadSim == 128 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise28_2020-09-29';
	elseif intLoadSim == 130 && boolLoad
	strSimulation = 'xAreaDistributed_LargeRetNoise30_2020-09-29';
	
elseif intLoadSim == 150 && boolLoad
	strSimulation = 'xAreaDistributed_IndInpOri5Noise5_2018-11-26';
	
	
elseif intLoadSim == 200 && boolLoad
	strSimulation = 'xAreaDistributed_OnlyFF_Noise00_2020-10-06';
end
