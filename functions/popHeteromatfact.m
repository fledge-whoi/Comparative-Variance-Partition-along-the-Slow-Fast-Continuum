function [Utilde, Ftilde, Atilde] = popHeteromatfact(thetamat, Rcomb, g, s)
%**********************************************************************
% Function that computes the transition matrix (Utilde), fertility matrix (Ftilde), and returns Utilde,Ftilde and Atilde (Atilde = Ftilde + Utilde).
% To construct the matrix population models we follow formulas presented in Caswell et al. 2018 (doi: 10.1002/ecm.1306). The code we implemented reflect that omega (the number of age classes) = 1

%INPUT
%thetamat: Matrix containing vital rates for each heterogeneity group. 
%Py: Template that maps low vital rates values in each heterogeneity group. 
%Rcomb: Phenotype distribution at birth
%g: number of groups
%s: number of stages
%OUTPUT
%...
%*************************************************************************

Is=eye(s);
Ig=eye(g);

% vec permutation matrix
p = zeros(s*g);
a = zeros(s,g);
for i = 1:s
    for j = 1:g
        e = a;
        e(i,j) = 1;
        p = p + kron(e,e');
    end
end
K=p;

H = repmat(Rcomb, [1, g]);
bbH=kron(Is(:,1)*Is(1,:),H);
for i=2:s
    bbH=bbH+kron(Is(:,i)*Is(i,:),H);
end

% 1: sj, 2: sa, 3: g, 4: F
bbU=kron(Ig(:,1)*Ig(1,:),[thetamat(1, 1)*(1-thetamat(1, 3)) 0; thetamat(1, 1)*thetamat(1, 3) thetamat(1, 2)]);
for i=2:g
    bbU=bbU+kron(Ig(:,i)*Ig(i,:),[thetamat(i, 1)*(1-thetamat(i, 3)) 0; thetamat(i, 1)*thetamat(i, 3) thetamat(i, 2)]);
end

% create block fertilty matrix
bbF=kron(Ig(:,1)*Ig(1,:),[0 thetamat(1, 4); 0 0]);
for i=2:g
    bbF=bbF+kron(Ig(:,i)*Ig(i,:),[0 thetamat(i, 4); 0 0]);
end

Utilde = bbU; % Fixed heterogeneity, so no transitions between groups and no age classes
Ftilde = K'*bbH*K*bbF;
Atilde = Utilde + Ftilde;
end
