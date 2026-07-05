% SETTING UP GENERAL VARIABLES
clear
close all
load("./data/simulation_results_Quality_Gconstant.mat")
addpath(genpath("./functions"))

[nr, nx, npop] = size(IHeteroLEresults);

q10line = 'magenta-';
q25line = 'yellow-';
meanline = 'b--';
medianline = 'b:';

cmap = flipud(hot); % reverses the colorbar
nc = length(cmap);

piepalette = [
  252, 255, 164;
  249, 142, 9;
  188, 55, 84;
  87, 16, 110;
  0, 0, 4
];
piepalette = piepalette / 256;

xresults = repmat(xvect, [length(Rvect)^4 1 npop]);
xresults = xresults + 1;
yresults = xresults - gapvalresults;
xyratio = xresults ./ yresults;
ncomb = length(Rvect)^(4);

titlefontsize = 40;
labelfontsize = 26;
tickfontsize = 14;
linewidthvalue = 0.75;
ticklengthrelativevalue = 0.05;
ticklengthvect = [ticklengthrelativevalue 0.025];

%% LE PIE
%..........................................................................
f = figure('Position',[10 10 2000 350], 'Visible','on');

indvect = 1:16; 
indextokeep = [8 12 14 15];
remain = [1 2 3 4 5 6 7 9 10 11 13];

for j=1:npop
    pos = j;
    subplot(1, npop, pos);

    newmat = LEbetweenvarlistresults(:,:,:,:,j);
    newmat = reshape(newmat, [16, ncomb * nx]);
    newmat(newmat < 0) = 0;
    size(newmat);

    Kmean = round(mean(IHeteroLEresults(:, :, j), "all", "omitnan"), 1);
    Kmean = Kmean / 1;

    sumvect = repelem(1, 16)' * sum(newmat, 1);
    newmat = newmat ./ sumvect;
    newmat(newmat == Inf) = NaN;

    sjalone = newmat(8, :);
    saalone = newmat(12, :);
    galone = newmat(14, :);
    falone = newmat(15, :);

    newmat = newmat(remain, :);
    interactions = sum(newmat, 1, "omitmissing") / 2^(4 - sum(heterogeneitymat));

    X = [mean(sjalone, "omitnan"), mean(saalone, "omitnan"),...
        mean(galone, "omitnan"), mean(falone, "omitnan"), mean(interactions, "omitnan")];

    if ((sum(X) ~= 0) && ~isnan(sum(X)))
        sum(X)
        X = X / sum(X)
        
        pc = piechart(X);
        pc.Labels = ["", "", "", "", ""];
        pc.FaceAlpha = 1;
        pc.EdgeColor = [0 0 0];
        pc.LineWidth = 2;
        colororder(viridis(5)) 
    end
end

%% LRO PIE
%..........................................................................
f = figure('Position',[10 10 2000 350], 'Visible','on');

indvect = 1:16;
indextokeep = [8 12 14 15];
remain = [1 2 3 4 5 6 7 9 10 11 13];

for j=1:npop
    pos = j;
    subplot(1, npop, pos);

    newmat = LRObetweenvarlistresults(:,:,:,:,j);
    newmat = reshape(newmat, [16, ncomb * nx]);
    newmat(newmat < 0) = 0;
    size(newmat);

    Kmean = round(mean(IHeteroLROresults(:, :, j), "all", "omitnan"), 1);
    Kmean = Kmean / 1;

    sumvect = repelem(1, 16)' * sum(newmat, 1);
    newmat = newmat ./ sumvect;
    newmat(newmat == Inf) = NaN;

    sjalone = newmat(8, :);
    saalone = newmat(12, :);
    galone = newmat(14, :);
    falone = newmat(15, :);

    newmat = newmat(remain, :);
    interactions = sum(newmat, 1, "omitmissing") / 2^(4 - sum(heterogeneitymat));

    X = [mean(sjalone, "omitnan"), mean(saalone, "omitnan"),...
        mean(galone, "omitnan"), mean(falone, "omitnan"), mean(interactions, "omitnan")];

    if ((sum(X) ~= 0) && ~isnan(sum(X)))
        sum(X)
        X = X / sum(X)
        
        pc = piechart(X);
        pc.Labels = ["", "", "", "", ""];
        pc.FaceAlpha = 1;
        pc.EdgeColor = [0 0 0];
        pc.LineWidth = 2;
        colororder(viridis(5))
    end
end

%% CODE for exploring the relantionship between Vb(LRO) and K(LRO)
% Code used to represent the "moving mean" for all six populations, in a
% single graph.
f = figure('Position',[10 10 350 350], 'Visible','on');
colorvect = viridis(6);

testmat = varLROresults .* IHeteroLROresults;

% Creating plot grid
nexttile(1)
hold on

for j=1:npop 
    % Processing data
    xvalues = testmat(:,:,j);
    yvalues = IHeteroLROresults(:,:,j);
    colorvalues = lambda1valsresults(:,:,j);
    xvalues = xvalues(:); yvalues = yvalues(:); colorvalues = colorvalues(:);

    % plotting the moving average:
    [sortedx, sortedxindexes] = sort(xvalues);
    sortedy = yvalues(sortedxindexes);
    ymovmean = 100*movmean(sortedy, 600); % computing the 600 point average

    pa = plot(sortedx, ymovmean, 'LineStyle', '-', 'Color', colorvect(j, :), 'Marker', '.', 'LineWidth', 2.5);

    % Formatting axes
    ax = gca;
    ax.LineWidth = 2; % Increases linewidth for all axes
    ax.TickLength = ticklengthvect * 0.5; % Increases tick length for all axes

    xlim([min(testmat, [], "all") max(testmat, [], "all")])
    xscale log
    ylim([0 max(IHeteroLROresults, [], "all")*1.05*100])

    XTickmin = round(log10(min(testmat, [], "all"))) - 1;
    XTickmax = round(log10(max(testmat, [], "all"))) + 1;
    XTick = XTickmin:XTickmax;
    XTickTen = 10.^XTick;
    XTickLabels = cellstr(num2str(round(log10(XTickTen(:))), '10^{%d}'));
    
    % Formatting labels
    xlabel('$V_b(LRO)$', 'Interpreter', 'latex', 'Rotation', 0, "FontSize", 18) 
    xtickangle(0)
    ylabel('$K(LRO)$', 'Interpreter', 'latex', 'Rotation', 90, "FontSize", 18)  

    % Setting fonts
    ax.XAxis.FontSize = tickfontsize * 0.9; % changes fontsize for label and ticks
    ax.XLabel.FontSize = labelfontsize; % changes fontsize for label only
    ax.YAxis.FontSize = tickfontsize; 
    ax.YLabel.FontSize = labelfontsize; 

    xticks(XTickTen);
    xticklabels(XTickLabels);
    xtickangle(45);
 
    grid off
end
legend('G = 2.05', 'G = 2.20', 'G = 3.43', 'G = 6.70', 'G = 13.85', 'G = 25.19', "Interpreter","latex", 'Fontsize', 0.8*labelfontsize)
hold off

%% CODE for exploring the relantionship between Vb(LE) and K(LE)
% Code used to represent the "moving mean" for all six populations, in a
% single graph.
f = figure('Position',[10 10 350 350], 'Visible','on');

testmat = varLEresults .* IHeteroLEresults;

% Creating plot grid
nexttile(1)
hold on

for j=1:npop 
    % Processing data
    xvalues = testmat(:,:,j);
    yvalues = IHeteroLEresults(:,:,j);
    colorvalues = lambda1valsresults(:,:,j);
    xvalues = xvalues(:); yvalues = yvalues(:); colorvalues = colorvalues(:);

    % plotting the moving average:
    [sortedx, sortedxindexes] = sort(xvalues);
    sortedy = yvalues(sortedxindexes);
    ymovmean = 100*movmean(sortedy, 600); % computing the 600 point average

    pa = plot(sortedx, ymovmean, 'LineStyle', '-', 'Marker', '.', 'LineWidth', 2.5);

    % Formatting axes
    ax = gca;
    ax.LineWidth = 2; % Increases linewidth for all axes
    ax.TickLength = ticklengthvect * 0.5; % Increases tick length for all axes

    xlim([min(testmat, [], "all") max(testmat, [], "all")])
    xscale log
    ylim([0 max(IHeteroLROresults, [], "all")*1.05*100])

    XTickmin = round(log10(min(testmat, [], "all"))) - 1;
    XTickmax = round(log10(max(testmat, [], "all"))) + 1;
    XTick = XTickmin:XTickmax;
    XTickTen = 10.^XTick;
    XTickLabels = cellstr(num2str(round(log10(XTickTen(:))), '10^{%d}'));
    
    % Formatting labels
    xlabel('$V_b(L)$', 'Interpreter', 'latex', 'Rotation', 0, "FontSize", 18) 
    xtickangle(0)
    ylabel('$K(L)$', 'Interpreter', 'latex', 'Rotation', 90, "FontSize", 18)  

    % Setting fonts
    ax.XAxis.FontSize = tickfontsize * 0.9; % changes fontsize for label and ticks
    ax.XLabel.FontSize = labelfontsize; % changes fontsize for label only
    ax.YAxis.FontSize = tickfontsize; 
    ax.YLabel.FontSize = labelfontsize; 

    xticks(XTickTen);
    xticklabels(XTickLabels);
    xtickangle(45);
     
    grid off
end
legend('G = 2.05', 'G = 2.20', 'G = 3.43', 'G = 6.70', 'G = 13.85', 'G = 25.19', "Interpreter","latex", 'Fontsize', 0.8*labelfontsize)
hold off
%% CODE for getting the quantiles, extrema and mean of K(L/LRO) for each population (QUALITY)
% Code used to represent the "moving mean" for all six populations, in a
% single graph.

dummymat = 100*permute(IHeteroLEresults, [3, 2, 1]);
maxvec = max(dummymat, [], [2, 3], "omitmissing");
minvec = min(dummymat, [], [2, 3], "omitmissing");
meanvec = mean(dummymat, [2, 3], "omitmissing");
medianvec = median(dummymat, [2, 3], "omitmissing");
D1vec = quantile(dummymat, 0.25, [2, 3]);
D9vec = quantile(dummymat, 0.75, [2, 3]);
dispvec = (D9vec - D1vec)./(maxvec - minvec);
disp("D1 median D9 max dispersion -- LE QUALITY")
disp([D1vec medianvec D9vec maxvec dispvec]) % populations are stored "vertically"

dummymat = 100*permute(IHeteroLROresults, [3, 2, 1]);
maxvec = max(dummymat, [], [2, 3], "omitmissing");
minvec = min(dummymat, [], [2, 3], "omitmissing");
meanvec = mean(dummymat, [2, 3], "omitmissing");
medianvec = median(dummymat, [2, 3], "omitmissing");
D1vec = quantile(dummymat, 0.25, [2, 3]);
D9vec = quantile(dummymat, 0.75, [2, 3]);
dispvec = (D9vec - D1vec)./(maxvec - minvec);
disp("D1 median D9 max dispersion -- LRO QUALITY")
disp([D1vec medianvec D9vec maxvec dispvec])

%% CORRELATION BETWEEN K(L) AND K(LRO)
A_vec = IHeteroLEresults(:);
B_vec = IHeteroLROresults(:);

mask = ~isnan(A_vec) & ~isnan(B_vec);   % keep only valid pairs

disp("Pearson correlation coefficient")
r = corr(A_vec(mask), B_vec(mask))

%% CODE for getting the quantiles, extrema and mean of K(L) for each population (TRADEOFFS)
% Code used to represent the "moving mean" for all six populations, in a
% single graph.

load("./data/simulation_results_SaSjTradeoff_LAMBDAconstant.mat")
dummymat = 100*permute(IHeteroLEresults, [3, 2, 1]);
maxvec = max(dummymat, [], [2, 3], "omitmissing");
minvec = min(dummymat, [], [2, 3], "omitmissing");
meanvec = mean(dummymat, [2, 3], "omitmissing");
medianvec = median(dummymat, [2, 3], "omitmissing");
D1vec = quantile(dummymat, 0.25, [2, 3]);
D9vec = quantile(dummymat, 0.75, [2, 3]);
dispvec = (D9vec - D1vec)./(maxvec - minvec);
disp("D1 median D9 max dispersion -- TRADEOF SA/SJ")
disp([D1vec medianvec D9vec maxvec dispvec]) % populations are stored "vertically"

load("./data/simulation_results_SaFTradeoff_LAMBDAconstant.mat")
dummymat = 100*permute(IHeteroLEresults, [3, 2, 1]);
maxvec = max(dummymat, [], [2, 3], "omitmissing");
minvec = min(dummymat, [], [2, 3], "omitmissing");
meanvec = mean(dummymat, [2, 3], "omitmissing");
medianvec = median(dummymat, [2, 3], "omitmissing");
D1vec = quantile(dummymat, 0.25, [2, 3]);
D9vec = quantile(dummymat, 0.75, [2, 3]);
dispvec = (D9vec - D1vec)./(maxvec - minvec);
disp("D1 median D9 max dispersion -- TRADEOF SA/F")
disp([D1vec medianvec D9vec maxvec dispvec])


