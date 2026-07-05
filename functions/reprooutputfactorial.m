function [meanLRO, varLRO] = reprooutputfactorial(U, F, N, selectionvect, s1, t0, Z)
%**********************************************************************
% Function that computes and returns the mean and variance of lifetime reproductive output. See. For details about the equations, see Van Daalen and Caswell 2017, doi: 10.1007/s12080-017-0335-2.
%INPUT
%U: The transition matrix
%F: The fertility matrix
%N: The fundamental matrix  
%s1, t0, Z: utility verctors and matrices to perform calculations.
%selectionvect: Vector containing the stage from which we compute expected longevity
%OUTPUT
%meanLRO: mean lifetime reproductive output for each heterogeneity group.
%varLRO: within-group variance in lifetime reproductive output for each heterogeneity group.
%*************************************************************************

% theta = sj1 sa1 g1 f1 sj2 sa2 g2 f2 r
% We consider a situation in which the only absorbing state is death

f = sum(F, 1); % Horizontal fertility vect
m = 1 - sum(U, 1); % Horizontal mortality vect
R1 = s1*[f 0]; % Rewards correspond to state-based reproductive output
R1tilde = Z*R1*Z';
P = [U t0; m 1]; % Markov chain with absorbing state

% We suppose a poisson distribution of offspring
R2f = (R1.^2);
% R3f = (R1.^3);

R2p = R1+R2f;
% R2ptilde = Z*R2p*Z';
% R3p = R1+3*R2f+R3f;
% R3ptilde = Z*R3p*Z';

% calculation of moments vectors rho for the lifetime reproductive output
% Formulas given in Van Daalen and Caswell 2017 (Theoretical Ecology) and
% Van Daalen and Caswell 2020 (Ecological modelling)
rhotilde1 = N'*Z*(P.*R1)'*s1;
rhotilde2 = N'*((Z*(P.*R2p)'*s1)+(2*(U.*R1tilde)'*rhotilde1));

rhotildemean = rhotilde1;
rhotildevar = rhotilde2-(rhotilde1.*rhotilde1);

% We did not gather juveniles and adults together in the model structure,
% so equations 26 and 27 from Van Daalen and Caswell 2020 cannot be applied
% directly, but the principle here is simple: we decide that the juvenile
% population is the starting point of our cohort
rhomean = rhotildemean(selectionvect);
rhovar = rhotildevar(selectionvect);

% partitioning of variance
meanLRO = rhomean;
varLRO = rhovar;
end
