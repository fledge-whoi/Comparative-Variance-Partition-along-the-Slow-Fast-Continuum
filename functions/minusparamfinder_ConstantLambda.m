function [diff] = minusparamfinder_ConstantLambda(thetax, Y, one_ngroup_nparam, Py, Rcomb)
%**********************************************************************
% Function that computes the difference between the population growth rate of any baseline population and the growth rate of a corresponding heterogeneous population (for which we know the phenotype distribution at birth and the vital rates of high quality individuals). 
%INPUT
%thetax: Vector containing the vital rate values for high-quality individuals.
%Y: Multiplicative factor that defines the vital rates of low quality individuals.
%one_ngroup_nparam: Matrix of dimension ngroupxnparam filled with ones. 
%Py: Template that maps low vital rates values in each heterogeneity group. 
%Rcomb: Phenotype distribution at birth
%OUTPUT
%diff: difference between genration time in the baseline population and in the 
%	heterogeneous population.
%*************************************************************************

Ymat = Y * one_ngroup_nparam;
thetaxy = thetax .* (Ymat .^ Py);

sjvect = thetaxy(:, 1);
savect = thetaxy(:, 2);
gvect = thetaxy(:, 3);
fvect = thetaxy(:, 4);

rzerovect = (fvect.*sjvect.*gvect)./((1 - savect).*(1 - (sjvect.*(1 - gvect))));

diff = (Rcomb'*rzerovect) - 1; % - 1 means we are looking for a stable population
	
end
