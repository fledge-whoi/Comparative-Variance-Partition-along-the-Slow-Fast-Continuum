function [VarX] = varcalc(X, PI)
%**********************************************************************
% Function that returns the population's growth rate, it's stable stage distribution.
%INPUT
%X: withingroup mean of any life history outcome
%PI: relative distribution of juvenile stages at birth.
%OUTPUT
%VarX: Between group variance.
%*************************************************************************

VarX = PI'*(X.*X) - (PI'*X).^2;
end
