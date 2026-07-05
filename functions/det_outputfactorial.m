function [lambda1, SSD, GT]= det_outputfactorial(Umat, Fmat, Amat)
%**********************************************************************
% Function that returns the population's growth rate, it's stable stage distribution.
%INPUT
%Umat: transition matrix
%Fmat: fertility matrix
%Amat: population matrix
%OUTPUT
%Lambda1: population's growth rate
%SSD: stable stage distribution
%GT: generation time (in years)
%*************************************************************************
[wmat, dmat, vmat]=eig(Amat);

% Calculate the rate of population change (lambda1 = the dominant eigenvalue) 
lambda = diag(dmat);
lambda1 = max(lambda);

imax=find(lambda==lambda1);
v = vmat(:,imax);
w = wmat(:,imax);

GT= 1 + ((v'*Umat*w)/(v'*Fmat*w)); 

RawAge=wmat(:,imax);
SSD=RawAge./sum(RawAge);

return
