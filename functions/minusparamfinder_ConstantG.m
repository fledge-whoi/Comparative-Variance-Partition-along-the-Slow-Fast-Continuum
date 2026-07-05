function [diff] = minusparamfinder_ConstantG(thetax, Y, one_ngroup_nparam, Py, Rcomb, G, ngroup, nstage)
%**********************************************************************
% Function that computes the difference between the generation time of any baseline population and the generation time of a corresponding heterogeneous population (for which we know the phenotype distribution at birth and the vital rates of high quality individuals). 
%INPUT
%thetax: Vector containing the vital rate values for high-quality individuals.
%Y: Multiplicative factor that defines the vital rates of low quality individuals.
%one_ngroup_nparam: Matrix of dimension ngroupxnparam filled with ones. 
%Py: Template that maps low vital rates values in each heterogeneity group. 
%Rcomb: Phenotype distribution at birth
%G: Generation time
%ngroup: number of heterogeneity groups 
%nstage: number of stages in the matrix model
%OUTPUT
%diff: difference between genration time in the baseline population and in the 
%	heterogeneous population.
%*************************************************************************

Ymat = Y * one_ngroup_nparam;
thetaxy = thetax .* (Ymat .^ Py);

[Uhet, Fhet, Ahet] = popHeteromatfact(thetaxy, Rcomb, ngroup, nstage);
[wmat, dmat, vmat]=eig(Ahet);
lambda=diag(dmat);
imax=find(lambda==max(lambda));
v = vmat(:,imax);
w = wmat(:,imax);
g = 1 + ((v'*Uhet*w)./(v'*Fhet*w));
diff = G - g;

end
