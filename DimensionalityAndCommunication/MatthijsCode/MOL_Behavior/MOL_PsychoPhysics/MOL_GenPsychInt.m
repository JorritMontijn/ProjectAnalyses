
function [psyvalues] = MOL_GenPsychInt(mu,stddev)
%generate 5 values for best sampling of psychometric curve:

stddev = abs(stddev);
stds = [-2 -1 0 1 2];
psyvalues = mu + stds*stddev;

end


