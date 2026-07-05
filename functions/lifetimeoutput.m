function [meanLE, varLE] = lifetimeoutput(N, I_ngroupstage, selectionvect)
%**********************************************************************
% Function that computes and returns the mean and variance of life expectancy. 
%INPUT
%N: The fundamental matrix  
%I_ngroupstage: An identity matrix of dimension ngroupstagexngroupstage
%selectionvect: Vector containing the stage from which we compute expected longevity
%OUTPUT
%meanLE: mean longevity for each heterogeneity group.
%varLE: within-group variance in longevity for each heterogeneity group.
%*************************************************************************

eta1tilde = sum(N, 1)';
eta2tilde = (eta1tilde'*(2*N - I_ngroupstage))';
varetatilde = eta2tilde - (eta1tilde.*eta1tilde);
meanLE = eta1tilde(selectionvect);
varLE = varetatilde(selectionvect);

end
