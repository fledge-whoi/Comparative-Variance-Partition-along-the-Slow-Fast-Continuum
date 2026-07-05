%% Choosing vital rates, and determining where heterogeneity is introduced
% Default code is for QUALITY heterogeneity. See lines 40, 173, 290, for
% adapting the code for TRADEOFF heterogeneity.

evalin('base','clear')
close all
addpath(genpath("./functions"))

% Number of demo parameters considered: sj sa g f
% Population demographic rates
thetamatpop = ... %sj sa g f
    [
     0.01  0.05  0.9    105.4500;
     0.2   0.1   0.6    6.9000;
     0.5   0.5   0.4    1.7500;
     0.7   0.75  0.1    1.3214;
     0.8   0.9   0.075  0.4333;
     0.85  0.95  0.05   0.2265
    ]; 

[npop, nparam] = size(thetamatpop);

nstage = 2; % Two stages: Juveniles and Adults; within each class
nclass = 2; % Number of classes per phenotypecombinations
heterogeneitymat = [1 1 1 1]; % We select demographic parameters on which to introduce individual heterogeneity
ngroup = nclass^nparam; % Total number of heterogeneity groups
heterogeneityfilter = repmat(heterogeneitymat, [ngroup 1]);

% Combinations of frequencies for each phenotype at birth that we want to
% try. We could choose different frequency vects for different parameters
Rvect = [0.01 0.1 0.25 0.5 0.75 0.9 0.99];
Rallcomb = combinations(Rvect, Rvect, Rvect, Rvect).Variables;
[ncomb, ~] = size(Rallcomb);

% Px: Combinations of parameter values in the groups (1 means greater parameter
% value than base model, 0 means lower parameter value than base model)
% QUALITY
% Pxunfiltered = combinations([1 0], [1 0], [1 0], [1 0]).Variables; % Sj Sa g f
% Pyunfiltered = 1 - Pxunfiltered; % Py: Complementary matrix of Px

%#############################################################################################################
% Note that there are other possible configurations for the combinations of parameter values. For instance:
% Only best and worst possible individuals
% Pxunfiltered = [repmat([1 1 1 1], [8 1]); repmat([0 0 0 0], [8 1])];

% Or the matrices to be used when introducing trade-off heterogeneity:
% TRADEOFF sa F
% Pxunfiltered = [repmat([0], [8 3]), repmat([1], [8 1]); repmat([0], [8 4])]; % Tradeoff sa F
% Pyunfiltered = [repmat([0], [8 1]), repmat([1], [8 1]), repmat([0], [8 2]); repmat([0], [8 4])]; % Tradeoff sa F

% TRADEOFF sj sa
Pxunfiltered = [repmat([0], [8 1]), repmat([1], [8 1]), repmat([0], [8 2]); repmat([0], [8 4])]; % Tradeoff sj sa
Pyunfiltered = [repmat([1], [8 1]), repmat([0], [8 3]); repmat([0], [8 4])]; % Tradeoff sj sa
%#############################################################################################################

% We modify Px and Py to only account for parameters on which we want to
% introduce heterogeneity
% The final Px matrix should remind you of figure 1 (from the method
% section of the article)
Px = Pxunfiltered .* heterogeneityfilter;
Py = Pyunfiltered .* heterogeneityfilter;

% For all combinations of sj,sa,g,f values, we create a Rmat that shows the
% frequencies of for each group at birth (we consider independance between
% phenotypes)
Rmat = zeros([ngroup nparam ncomb]);
for comb=1:ncomb
    Rj = Rallcomb(comb, 1);
    Ra = Rallcomb(comb, 2);
    Rg = Rallcomb(comb, 3);
    Rf = Rallcomb(comb, 4);
    % Rmatcombunfiltered = combinations([Rj 1-Rj], [Ra 1-Ra], [Rg 1-Rg], [Rf 1-Rf]).Variables;
    % Rmatcombfiltered = Rmatcombunfiltered.*heterogeneityfilter + (1 - heterogeneityfilter); % We remove groups we do not need
    % Rmat(:, :, comb) = Rmatcombfiltered;
    Rmat(:, :, comb) = combinations([Rj 1-Rj], [Ra 1-Ra], [Rg 1-Rg], [Rf 1-Rf]).Variables;
end

%% Amplitude of heterogeneity
% Here we choose the demographic parameter value in the (+) group first,
% then, latter in the process, we calculate the value of parameters in the (-) group
% according to specific criteria (either constante growth rate or
% generation time)
% demoratetocheck is a variable that simply allows us to determine the
% vital rates to account for when calculating the maximum increase in vital
% rate. Here we constructed it in a way that includes all parameters from
% heterogeneitymat and excludes parameters you don't want (set their value
% to 0, for instance here we excluded fertility because it is not a transition rate)
transitionrate = [1 1 1 0];
demoratetocheck = transitionrate.*heterogeneitymat;
demoratetocheck = logical(demoratetocheck);
xmaxmat = 1 ./ thetamatpop(:, demoratetocheck) - 1; % We try to get the maximum increase in vital rates
% for each species (since vital rates cannot exceed one).

% While xmat is the maximum threshold for the vital rates of high-quality individuals
% nx is the number of vital rates we want to test for between baseline parameters and 
% the threshold. Change this parameter according to the resolution you need for exploring the parameter space.  
xmat = min(xmaxmat, [], 2);
nx = 50;

%% Basic utilities and population statistics

% Useful vectors and matrices to conduct matrix calculations
one_nclass = ones([nclass 1]);
one_nparam = ones([nparam 1]);
one_ngroup = ones([ngroup 1]);
one_ngroup_nparam = ones([ngroup nparam]);
I_ngroupstage = eye(ngroup*nstage);
I_nstage = eye(nstage);

selectionvect = repmat(logical([1 0]), [1 ngroup]); % We choose a stage for which 
% to compute the expected longevity and lifetime reproductive output (here juveniles)

% Useful vectors and matrices to conduct calculations for the rewards
t = nstage*ngroup;
s = t+1;
s1 = ones([s 1]);
t1 = ones([t 1]);
t0 = zeros([t 1]);
It = diag(t1);
Z = [It t0];

% Options to pass to fzero function
% options = optimset('Display','iter'); % show iterations
options = optimset();

% Getting basic information on the homogeneous populations with baseline vital rates. 
GTbasevect = NaN([1, npop]);
Rcomb = 1;
initialstats = NaN([4 npop]);
for pop=1:npop
    if (~anynan(thetamatpop(pop, :)))
        thetamat = thetamatpop(pop, :);
        [Uhet, Fhet, Ahet] = popHeteromatfact(thetamat, Rcomb, 1, nstage);
        [wmat, dmat, vmat]=eig(Ahet);
        lambda=diag(dmat);
        imax=find(lambda==max(lambda));
        v = vmat(:,imax);
        w = wmat(:,imax);
        gt= 1 + ((v'*Uhet*w)./(v'*Fhet*w));
        GTbasevect(1, pop) = round(gt, 2);
        
	    %Ntilde is the fundamental matrix (see Caswell et al. 2018, doi:10.1002/ecm.1306)	
        Ntilde = (eye(2)-Uhet)\eye(2);
        [meanLE, varLE] = lifetimeoutput(Ntilde, eye(2), logical([0 1]));
        [meanLRO, varLRO] = reprooutputfactorial(Uhet, Fhet, Ntilde, logical([1 0]), ones([3 1]), zeros([2 1]), [diag(ones([2 1])) zeros([2 1])]);
    
        initialstats(1, pop) = meanLE; initialstats(2, pop) = varLE; 
        initialstats(3, pop) = meanLRO; initialstats(4, pop) = varLRO;
    end
end

% Results storage
[sortedPx, hetindex, ~] = unique(Px, "rows", "first"); % We get all unique combinations of heterogeneity groups
[nhet, ~] = size(sortedPx);

gapvalresults = NaN([ncomb, nx, npop]); %1x1
lambda1valsresults = NaN([ncomb, nx, npop]); %1x1
SSDresults = NaN([ngroup*nstage, 1, ncomb, nx, npop]); %ngroupnstagex1
PIresults = NaN([1, ngroup, ncomb, nx, npop]); %1xngroup
GTresults = NaN([ncomb, nx, npop]); %1x1 

meanLEresults = NaN([ncomb, nx, npop]); %1x1
varLEresults = NaN([ncomb, nx, npop]); %1x1
IHeteroLEresults = NaN([ncomb, nx, npop]); %1x1
LEbetweenvarlistresults = NaN([ngroup, 1, ncomb, nx, npop]); %1xvariable

meanLROresults = NaN([ncomb, nx, npop]); %1x1
varLROresults = NaN([ncomb, nx, npop]); %1x1
IHeteroLROresults = NaN([ncomb, nx, npop]); %1x1
LRObetweenvarlistresults = NaN([ngroup, 1, ncomb, nx, npop]); %1xvariable


%% Core loops
% Following code is written to reproduce results for QUALITY heterogeneity
% If producing results for TRADE-OFFS, change:
% "minusparamfinder_ConstantG" by "minusparamfinder_ConstantLambda" at line 207. 
% Please notice the arguments of these two functions are different.
% You may also want to change the name of the save file at line 290.


tic % setting a timer
for pop=1:npop
    theta = thetamatpop(pop, :)
    thetamat = one_ngroup * theta;
    
    % We set the amount of heterogeneity amplitude we need to apply
    xvect = linspace(0, xmaxmat(pop), nx); 

    if (~anynan(theta))
        for xind=1:nx
            % Define an X value for the phenotype with better quality
            X = 1 + xvect(xind);
            Xmat = X * one_ngroup_nparam;
	        % thetax contains the vital rate values for high-quality individuals.
            thetax = thetamat .* (Xmat .^ Px);

            % We check that probabilities / rates are not greater than 1
            % or lower than 0. We ignore the 4th parameter since it is
            % fertility
            if (~any([1 1 1 1] - [thetax(1, 1:3) 0] <= 1e-15)) && (~any(thetax(1, 1:4) < 0))
    
                for comb=1:ncomb
                    Rmatsmall = Rmat(:, :, comb);
                    % Find the corresponding Y value to keep the population
                    % stable or with the same generation time
		            % Rcomb is the phenotype distribution at birth.
                    Rcomb = Rmatsmall(:, 1, 1) .* Rmatsmall(:, 2, 1) .* Rmatsmall(:, 3, 1) .* Rmatsmall(:, 4, 1);
                    % f = @(Y) minusparamfinder_ConstantG(thetax, Y, one_ngroup_nparam, Py, Rcomb, GTbasevect(pop), ngroup, nstage);
                    f = @(Y) minusparamfinder_ConstantLambda(thetax, Y, one_ngroup_nparam, Py, Rcomb);
                    
                    try 
                        Y = fzero(f, [0 1], options); % When X > 1
                    catch ME
                        % Checks if a solution does not exist between 0 and 1
                        if (strcmp(ME.identifier, 'MATLAB:fzero:ValuesAtEndPtsSameSign'))
                            Y = NaN;
                        else
                            disp(ME.identifier)
                        end
                    end    

                    Ymat = Y * one_ngroup_nparam;
                    thetaxy = thetax .* (Ymat .^ Py);
        
                    if ((~anynan(Y)) && (~any(thetax(16, 1:4) < 0)))
                        % We construct the matrix population model
                        [Uhet, Fhet, Ahet] = popHeteromatfact(thetaxy, Rcomb, ngroup, nstage);
                        [lambda1, SSD, GT] = det_outputfactorial(Uhet, Fhet, Ahet);
                        
                        % Generate the fundamental matrix
                        Ntilde = (I_ngroupstage-Uhet)\I_ngroupstage;
        
                        % Calculate stable distribution of the groups for a
                        % particular stages
                        PI = SSD(selectionvect);
                        PI = PI / sum(PI);
                        reshapedPI = reshape(PI, nclass*one_nparam');
        
                        % Calculate mean and variance in longevity
                        [meanLE, varLE] = lifetimeoutput(Ntilde, I_ngroupstage, selectionvect);
                        
                        % Variance decomposition in longevity
                        LEwithinvar = PI'*varLE;
                        LEbetweenvartot = varcalc(meanLE, PI);
                        IHeteroLE = LEbetweenvartot./(LEwithinvar + LEbetweenvartot);
                             
                        reshapedmeanLE = reshape(meanLE, nclass*one_nparam');
                        LEbetweenvarlist = vardecomposition(reshapedmeanLE, reshapedPI, Px, Py, ngroup, heterogeneitymat);
                        LEbetweenvarlist = LEbetweenvarlist';
        
                        % Calculate mean and variance in lifetime reproductive
                        % success
                        [meanLRO, varLRO] = reprooutputfactorial(Uhet, Fhet, Ntilde, selectionvect, s1, t0, Z);
                        
                        % Variance decomposition in lifetime reproductive
                        % success
                        LROwithinvar = PI'*varLRO;
                        LRObetweenvartot = varcalc(meanLRO, PI);
                        IHeteroLRO = LRObetweenvartot./(LROwithinvar + LRObetweenvartot);
                        
                        reshapedmeanLRO = reshape(meanLRO, nclass*one_nparam');
                        LRObetweenvarlist = vardecomposition(reshapedmeanLRO, reshapedPI, Px, Py, ngroup, heterogeneitymat);
                        LRObetweenvarlist = LRObetweenvarlist';
                       
                        % Saving results
                        %array_to_save_results_to %nxn = dimension of the data
                        %that is saved
                        gapvalresults(comb, xind, pop) = X - Y; %1x1
                        lambda1valsresults(comb, xind, pop) = lambda1; %1x1
                        SSDresults(:, :, comb, xind, pop) = SSD; %ngroupnstagex1
                        PIresults(:, :, comb, xind, pop) = PI; %1xngroup
                        GTresults(comb, xind, pop) = GT; %1x1 
        
                        meanLEresults(comb, xind, pop) = PI'*meanLE; %1x1
                        varLEresults(comb, xind, pop) = LEbetweenvartot + LEwithinvar; %1x1
                        IHeteroLEresults(comb, xind, pop) = IHeteroLE; %1x1
                        LEbetweenvarlistresults(:, :, comb, xind, pop) = LEbetweenvarlist; %variablex1
        
                        meanLROresults(comb, xind, pop) = PI'*meanLRO; %1x1
                        varLROresults(comb, xind, pop) = LRObetweenvartot + LROwithinvar; %1x1
                        IHeteroLROresults(comb, xind, pop) = IHeteroLRO; %1x1
                        LRObetweenvarlistresults(:, :, comb, xind, pop) = LRObetweenvarlist; %variablex1
                    end
                end
            end
        end
    end
end
toc

%% EXECUTE THIS SECTION TO SAVE YOUR RESULTS

% save("./data/simulation_results_Quality_Gconstant.mat", ...
%     "xvect", "Rvect", "gapvalresults", "varLEresults", "meanLROresults",...
%     "meanLEresults", "varLEresults", "varLROresults", "IHeteroLEresults", "IHeteroLROresults", "LEbetweenvarlistresults", ...
%     "LRObetweenvarlistresults", "Px", "Py", "thetamatpop", "lambda1valsresults",...
%     "GTresults", "SSDresults", "PIresults", "GTbasevect", "initialstats", "heterogeneitymat")

% save("./data/simulation_results_SaFTradeoff_LAMBDAconstant.mat", ...
%     "xvect", "Rvect", "gapvalresults", "varLEresults", "meanLROresults",...
%     "meanLEresults", "varLEresults", "varLROresults", "IHeteroLEresults", "IHeteroLROresults", "LEbetweenvarlistresults", ...
%     "LRObetweenvarlistresults", "Px", "Py", "thetamatpop", "lambda1valsresults",...
%     "GTresults", "SSDresults", "PIresults", "GTbasevect", "initialstats", "heterogeneitymat")

save("./data/simulation_results_SaSjTradeoff_LAMBDAconstant.mat", ...
    "xvect", "Rvect", "gapvalresults", "varLEresults", "meanLROresults",...
    "meanLEresults", "varLEresults", "varLROresults", "IHeteroLEresults", "IHeteroLROresults", "LEbetweenvarlistresults", ...
    "LRObetweenvarlistresults", "Px", "Py", "thetamatpop", "lambda1valsresults",...
    "GTresults", "SSDresults", "PIresults", "GTbasevect", "initialstats", "heterogeneitymat")

