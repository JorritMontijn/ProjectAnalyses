function [matL, vecMu, vecPsi, vecLLH] = factorAnal(X, zDim, dblTolerance)
	% Perform EM algorithm for factor analysis model
	% Input:
	%   X: n x d data matrix
	%   zDim: dimension of target space
	%	dblTolerance: tolerance for fitting convergence
	% Output:
	%   matL: d x m weight matrix
	%   vecMu: d x 1 mean vector
	%   vecPsi: d x 1 variance vector
	%   vecLLH: loglikelihood
	% Written by Jorrit Montijn, based on code by Byron Yu, based on code by Zoubin Ghahramani.
	
	%% input parameters
	minVarFrac = -inf;
	if ~exist('dblTolerance','var')
		tol = 1e-8;
	else
		tol = dblTolerance;
	end
	intIterMax = round(1/tol);
	verbose = false;
	
	%reformat data
	X = X';
	[xDim, N] = size(X);
	vecMu = mean(X,2);
	X = bsxfun(@minus,X,vecMu);
	
	%% based on fastfa
	% Initialization of parameters
	cX    = cov(X', 1);
	if rank(cX) == xDim
		scale = exp(2*sum(log(diag(chol(cX))))/xDim);
	else
		% cX may not be full rank because N < xDim
		fprintf('WARNING in fastfa.m: Data matrix is not full rank.\n');
		r     = rank(cX);
		e     = sort(eig(cX), 'descend');
		scale = geomean(e(1:r));
	end
	matL     = randn(xDim,zDim)*sqrt(scale/zDim);
	vecPsi    = diag(cX);
	vecMu     = mean(X, 2);
	
	varFloor = minVarFrac * diag(cX);
	
	I     = eye(zDim);
	const = -xDim/2*log(2*pi);
	LLi   = 0;
	vecLLH    = nan(1,intIterMax);
	
	for i = 1:intIterMax
		% =======
		% E-step
		% =======
		iPh  = diag(1./vecPsi);
		iPhL = iPh * matL;
		MM   = iPh - iPhL / (I + matL' * iPhL) * iPhL';
		beta = matL' * MM; % zDim x xDim
		
		cX_beta = cX * beta'; % xDim x zDim
		EZZ     = I - beta * matL + beta * cX_beta;
		
		% Compute log likelihood
		LLold = LLi;
		ldM   = sum(log(diag(chol(MM))));
		LLi   = N*const + N*ldM - 0.5*N*sum(sum(MM .* cX));
		if verbose
			fprintf('EM iteration %5i lik %8.1f \r', i, LLi);
		end
		vecLLH(i) = LLi;
		
		% =======
		% M-step
		% =======
		matL  = cX_beta / EZZ;
		vecPsi = diag(cX) - sum(cX_beta .* matL, 2);
		
		% Set minimum private variance
		vecPsi = max(varFloor, vecPsi);
		
		if i<=2
			LLbase = LLi;
		elseif (LLi < LLold)
			disp('VIOLATION');
		elseif ((LLi-LLbase) < (1+tol)*(LLold-LLbase))
			break;
		end
	end
	
	if verbose
		fprintf('\n');
	end
	
	if any(vecPsi == varFloor)
		fprintf('Warning: Private variance floor used for one or more observed dimensions in FA.\n');
	end
	
	%remove over-allocated entries
	vecLLH(isnan(vecLLH)) = [];
	
	%%
	%{
llh = -inf(1,maxiter);

I = eye(m);
r = dot(X,X,2);

W = randn(d,m);
lambda = 1./rand(d,1);
for iter = 2:maxiter
    T = bsxfun(@times,W,sqrt(lambda));
    M = T'*T+I;                     % M = W'*inv(Psi)*W+I
    U = chol(M);
    WInvPsiX = bsxfun(@times,W,lambda)'*X;       % WInvPsiX = W'*inv(Psi)*X
    
    % likelihood
    logdetC = 2*sum(log(diag(U)))-sum(log(lambda));              % log(det(C))
    Q = U'\WInvPsiX;
    trInvCS = (r'*lambda-dot(Q(:),Q(:)))/n;  % trace(inv(C)*S)
    llh(iter) = -n*(d*log(2*pi)+logdetC+trInvCS)/2;
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter-1)); break; end   % check likelihood for convergence
    
    % E step
    Ez = M\WInvPsiX;                                         % 12.66
    V = inv(U);
    Ezz = n*(V*V')+Ez*Ez';                                        % 12.67
    
    % M step
    U = chol(Ezz);
    XEz = X*Ez';
    W = (XEz/U)/U';                                         % 12.69
    lambda = n./(r-dot(W,XEz,2));                           % 12.70
end
llh = llh(2:iter);
psi = 1./lambda;

if iter == maxiter
	warning([mfilename ':NoConvergence'],'Factor analysis did not converge; maximum number of steps used');
end
	%}
