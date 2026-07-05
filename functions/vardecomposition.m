function [sortedvarfactoriallist] = vardecomposition(M, PI, Px, Py, ngroup, excluded)
%**********************************************************************
% This is the code for factorial variance partition. The code is abstract in
% order to apply to all heterogeneity structures in the populations 
% (i.e. all possible structures of the Px matrices). 
 
%INPUT
%M: Matrix with mean values of any life history outcome. 
%PI: Matrix with relative distribution of juvenile stages at birth.
%Px, Py: Template matrices that allow us to differenciate between heterogeneity groups 
%ngropu: number of groups
%excluded: Vector stores the vital rates for which there is or is not heterogeneity.
%OUTPUT
%VarX: Between group variance.
%*************************************************************************

% Here comb (for combination) refers to Py (not Px) because we sum on dimensions that are
% NOT part of the combination of heterogeneity groups we use. 
% M and PI should be in comptact form (not a vectorized form)

% Here we permute because reshape flipped the dimensions
% for the parameters (e.g. sj corresponds to the 4th
% dimension instead of the 1st and conversely f corresponds
% to the 1st dimension instead of the 4th.
comb = flip(Py, 2);
nparam = length(comb(1, :));
paramind = 1:nparam;
varfactoriallist = zeros([1, ngroup]);
MPI = M .* PI;

for i=1:ngroup
    contractionvect = paramind(logical(comb(i, :)));
    if ~isempty(contractionvect)
        PIpartial = sum(PI, contractionvect);
        Mpartial = sum(MPI, contractionvect) ./ PIpartial;
    else
        Mpartial = M;
        PIpartial = PI;
    end
    varfactoriallist(i) = varcalc(Mpartial(:), PIpartial(:));
end    

varfactoriallist(isnan(varfactoriallist))=0;

% We remove groups that are identical so we can automatize variance
% partition calculations
[sortedPx, hetindex, ~] = unique(Px, "rows", "first");
sortedPx = flip(sortedPx, 1);
hetindex = flip(hetindex, 1);
sortedvarfactoriallist = varfactoriallist(hetindex);

[~, newngroup] = size(sortedvarfactoriallist);
degreesvect = sum(sortedPx, 2);
one_ngroup = ones([newngroup 1]);
indexes = 1:newngroup;

% We calculate variance for interactions (more than one factor)s
for deg=2:max(degreesvect)
    degreepos = logical(degreesvect == deg); 
    degreeindex = indexes(degreepos);
    % We get all subcombinations of factors (lower degrees)
    indexvarsubdegree = sum(degreesvect == 0:(deg - 1), 2);

    for ind=degreeindex
        % We all factor combinations with a common factor with heterogeneity
        combi = one_ngroup * sortedPx(ind, :); % Now we need sortedPx
        indexvarcommonfact = sum(abs(sortedPx - combi) - combi, 2) + degreesvect == 0;
        
        % We combine the two pieces of information
        indexvartotake = indexvarcommonfact .* indexvarsubdegree;
        vartotake = sortedvarfactoriallist(logical(indexvartotake));
        sumvartotake = sum(vartotake);
        
        % Be aware that when variances are low, are very close to one
        % another the substraction sometimes yield results 
        % that are very very close to zero, but negative.
        sortedvarfactoriallist(ind) = sortedvarfactoriallist(ind) - sumvartotake;
    end
end

% We reinsert sortedvarfactoriallist into a 16x1 vector to ease code
% implementation down the road.

% This code does not work yet for more than 1 excluded parameter
initialpowers = 2.^(4 - find(excluded == 0));
powers = [];
supplement = NaN([1, ngroup]);
queue = 1:ngroup;

for i=1:length(initialpowers)
    dummy = nchoosek(initialpowers, i);
    sumdummy = sum(dummy, 2);
    powers = [powers, sumdummy'];
end

for i=1:newngroup
    k = queue(1);
    supplement(k) = sortedvarfactoriallist(i); 
    queue = queue(queue ~= k);
    for p=powers
        % p
        kbis = k + p;
        queue = queue(queue ~= kbis);
        supplement(kbis) = sortedvarfactoriallist(i);
    end
end

% supplement
supplement(isnan(supplement)) = 0;
sortedvarfactoriallist = supplement;
