%==========================================================================
% SCRIPT DESCRIPTION
%==========================================================================
%
% PURPOSE:
% Processes building-level seismic analysis data ('your_data_file.csv') and 
% generates lognormal seismic fragility curves.
%
% METHOD:
% Uses the Hueste & Bai (2004) lognormal CDF formulation with standard 
% uncertainty components (β_CL = 0.30, β_M = 0.30). The drift–PGV relationship 
% for each building is obtained via a power-law regression:
%       ln(D) = a + b·ln(PGV)
%
% MAIN OPERATIONS:
% • Performs building-by-building regression and produces associated plots.
% • Computes fragility parameters (C, k, β_fs) and limit-state capacities.
% • Exports key results to:
%       - 'Fragility_Summary.csv'
%       - 'Median_PGV_Table.csv'
%       - 'Regression_Stats.csv'
%
% GROUP ANALYSIS:
% • Classifies buildings into Groups A, B, and C.
% • Computes group-average model and capacity parameters.
% • Generates representative fragility curves for each group.
%
%==========================================================================

%% --- 1. Initialization ---
clc; close all;
warning('off','all');

% Load the data
data = readtable('your_data_file.csv'); % Adjust the file name accordingly

% Extract unique building names from folder names
buildingNames = unique(extractBetween(data.FolderName, 'Data_build_', '_EQ')); 
buildingNums = str2double(buildingNames);
[~, sortIdx] = sort(buildingNums);
buildingNames = buildingNames(sortIdx);
numBuildings = length(buildingNames);

% Define uncertainty parameters
beta_CL = 0.30;   % Capacity uncertainty
beta_M = 0.30;    % Modeling uncertainty

%% --- 2. Function Definitions and Folder Setup ---

% Define functions used in the Hueste & Bai (2004) methodology
local_power_fit = @(PGV,DR) powerlaw_fit_impl(PGV,DR);
lambdaD = @(a,b,pgv) a + b.*log(pgv); % ln(predicted drift)
beta_total_fun = @(beta_fs,bL,bM) sqrt(beta_fs.^2 + bL.^2 + bM.^2);
Pfrag = @(pgv,a,b,thetaL,beta_fs,bL,bM) 1 - normcdf( (log(thetaL) - lambdaD(a,b,pgv)) ./ beta_total_fun(beta_fs,bL,bM));
pgv50 = @(a,b,thetaL) exp( (log(thetaL) - a) / b ); % Solution for P = 0.5 (median)

% Create folders for figures
if ~exist('fragility_figs','dir'); mkdir fragility_figs; end
if ~exist('driftvspgv_figs','dir'); mkdir driftvspgv_figs; end
if ~exist('powerlawfit_fits','dir'); mkdir powerlawfit_fits; end 

%% --- 3. Individual Building Analyses ---

Summary = table('Size',[0 10], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Building','Scenario','C','k','R2','s2','beta_fs','LS1_pct','LS2_pct','LS3_pct'});

RegStats = table('Size',[0 6],'VariableTypes', {'string','string','double','double','double','double'},'VariableNames', {'Building','Scenario','s2','beta_fs','R2_log','n_points'});
MedianSummary = table('Size',[0 4],'VariableTypes', {'string','string','double','double'},'VariableNames', {'Building','LimitState','PGV50_EQ','PGV50_Seq'});

for i = 1:numBuildings
    % Identify building name pattern
    buildingPattern = sprintf('Data_build_%s_EQ', buildingNames{i});
    bLabel = sprintf('B%02d', i); % Name B01, B02...
    
    % Filter EQ and Seq data
    EQData = data(startsWith(data.FolderName, buildingPattern) & ~contains(data.FolderName, '_Seq'), :);
    SeqData = data(startsWith(data.FolderName, buildingPattern) & contains(data.FolderName, '_Seq'), :);
    hasSeq = ~isempty(SeqData);

    if isempty(EQData)
        warning('Bina %s için EQData bulunamadı, atlanıyor.', bLabel);
        continue;
    end
    
    % Get height and LS values for this building (assumed constant)
    Height = EQData.Height(1);
    LS1_pct = (EQData.LS1(1) / Height) * 100; % Convert to %
    LS2_pct = (EQData.LS2(1) / Height) * 100; % Convert to %
    LS3_pct = (EQData.LS3(1) / Height) * 100; % Convert to %
    
    % Calculate drift ratios in percentage
    Drift_EQ = (EQData.MaxTopDisplacement ./ Height) * 100;
    Drift_Seq = (SeqData.MaxTopDisplacement ./ Height) * 100;
    
    % --- PGV vs Drift Graph ---
    f_drift = figure('Visible','off'); 
    ax_drift = axes('Parent', f_drift); 
    
    scatter(ax_drift, EQData.PGV, Drift_EQ, 'kx'); 
    hold(ax_drift, 'on');
    if hasSeq
        scatter(ax_drift, SeqData.PGV, Drift_Seq, 'ks'); 
    end
    yline(ax_drift, LS1_pct, 'k--', 'LineWidth', 2, 'Label', 'LS1');
    yline(ax_drift, LS2_pct, 'k--', 'LineWidth', 2, 'Label', 'LS2');
    yline(ax_drift, LS3_pct, 'k--', 'LineWidth', 2, 'Label', 'LS3');
    set(ax_drift, 'YScale', 'log'); 
    xlabel(ax_drift, 'PGV (cm/s)'); ylabel(ax_drift, 'Drift Ratio (%) [log scale]');
    
    title(ax_drift, sprintf('PGV vs Drift for Building %s', bLabel)); 
    
    grid(ax_drift, 'on');  
    
    % Y-Limit Alignment (Step 1: Retrieve Limits)
    y_limits_log = ylim(ax_drift); 
    
    saveas(f_drift, fullfile('driftvspgv_figs', sprintf('Building_%s_PGV_vs_Drift.png', buildingNames{i})));
    close(f_drift);
    
    % Log-Log Fit
    [a_EQ,b_EQ,beta_fs_EQ,R2_EQ,n_EQ,s2_EQ] = local_power_fit(EQData.PGV, Drift_EQ);
    C_EQ = exp(a_EQ);  k_EQ = b_EQ; 
    
    if hasSeq
        [a_SQ,b_SQ,beta_fs_SQ,R2_SQ,n_SQ,s2_SQ] = local_power_fit(SeqData.PGV, Drift_Seq);
        C_SQ = exp(a_SQ);  k_SQ = b_SQ;
    else
        a_SQ=NaN; b_SQ=NaN; beta_fs_SQ=NaN; R2_SQ=NaN; n_SQ=0;
        C_SQ=NaN; k_SQ=NaN;
    end

    PGV_all = [EQData.PGV; SeqData.PGV]; 
    PGV_all = PGV_all(PGV_all > 0 & isfinite(PGV_all));
    if isempty(PGV_all)
        pmin_data = 1e-3; pmax_data = 150; 
    else
        pmax_data = max(PGV_all); 
        pmin_data = min(PGV_all);
    end
    PGV_grid_fit = linspace(0.1, 150, 400); 
    pmin_frag = max(pmin_data * 0.8, 1e-3);
    pmax_frag = min(pmax_data * 1.2, 150); 
    PGV_grid_frag = linspace(pmin_frag, pmax_frag, 400);
    f_fit = figure('Visible','off'); 
    ax_fit = axes('Parent', f_fit); 
    hold(ax_fit, 'on'); grid(ax_fit, 'on'); 

    % Plot the fitted curve
    Drift_fit_EQ = exp(a_EQ + b_EQ .* log(PGV_grid_fit));
    % Legend
    eq_display_name = sprintf('EQ: y = %.4fx^{%.4f}\nR%c = %.3f', C_EQ, k_EQ, 178, R2_EQ);
    p1 = plot(ax_fit, PGV_grid_fit, Drift_fit_EQ, 'r-', 'LineWidth', 2, 'DisplayName', eq_display_name);
    
    % Plot the raw data for the EQ dataset
    scatter(ax_fit, EQData.PGV, Drift_EQ, 'kx');
    
    legend_items = [p1]; 
    
    if hasSeq
        % Plot the fitted curve for the Seq dataset
        Drift_fit_Seq = exp(a_SQ + b_SQ .* log(PGV_grid_fit));
        seq_display_name = sprintf('Seq: y = %.4fx^{%.4f}\nR%c = %.3f', C_SQ, k_SQ, 178, R2_SQ);
        p2 = plot(ax_fit, PGV_grid_fit, Drift_fit_Seq, 'b--', 'LineWidth', 2, 'DisplayName', seq_display_name);
        
        % Plot the raw data for the EQ dataset
        scatter(ax_fit, SeqData.PGV, Drift_Seq, 'ks');
        legend_items = [legend_items, p2]; 
    end
    
    % Y axis (Log Scale)
    set(ax_fit, 'YScale', 'log'); 
    xlabel(ax_fit, 'PGV (cm/s)'); 
    ylabel(ax_fit, 'Drift Ratio (%) [log scale]'); 
    legend(ax_fit, legend_items, 'Location', 'southeast', 'FontSize', 8);
    xlim(ax_fit, [0 150]); 
    ylim(ax_fit, y_limits_log); 
    box(ax_fit, 'off'); 
    
    % Building Labels Added
    text(ax_fit, 1, y_limits_log(2)*0.9, bLabel, ... 
         'Units','data', 'FontWeight','bold','HorizontalAlignment','left', ...
         'VerticalAlignment','top','FontSize',11, 'FontName','Times New Roman');

    saveas(f_fit, fullfile('powerlawfit_fits', sprintf('Fit_Plot_%s.png', buildingNames{i})));
    close(f_fit);
    
    % Individual Fragility Plot
    f2 = figure('Visible','off'); ax2 = axes('Parent', f2); 
    hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
    plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_EQ, b_EQ, LS1_pct, beta_fs_EQ, beta_CL, beta_M), '-',  'LineWidth',1.8, 'DisplayName','EQ–LS1');
    plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_EQ, b_EQ, LS2_pct, beta_fs_EQ, beta_CL, beta_M), '-',  'LineWidth',1.8, 'DisplayName','EQ–LS2');
    plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_EQ, b_EQ, LS3_pct, beta_fs_EQ, beta_CL, beta_M), '-',  'LineWidth',1.8, 'DisplayName','EQ–LS3');
    if hasSeq
        plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_SQ, b_SQ, LS1_pct, beta_fs_SQ, beta_CL, beta_M), '--', 'LineWidth',1.8, 'DisplayName','Seq–LS1');
        plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_SQ, b_SQ, LS2_pct, beta_fs_SQ, beta_CL, beta_M), '--', 'LineWidth',1.8, 'DisplayName','Seq–LS2');
        plot(ax2, PGV_grid_frag, Pfrag(PGV_grid_frag, a_SQ, b_SQ, LS3_pct, beta_fs_SQ, beta_CL, beta_M), '--', 'LineWidth',1.8, 'DisplayName','Seq–LS3');
    end
    xlabel(ax2, 'PGV (cm/s)'); ylabel(ax2, 'P(Exceedance)');   
    xlim(ax2, [0 pmax_frag]); 
    yl = ylim(ax2); xl = xlim(ax2);
    text(ax2, xl(1) + (xl(2)-xl(1))*0.02, yl(2)*0.98, bLabel,'Units','data', 'FontWeight','bold','HorizontalAlignment','left', 'VerticalAlignment','top','FontSize',11, 'FontName','Times New Roman');
    legend(ax2, 'Location','southeast');
    saveas(f2, fullfile('fragility_figs',sprintf('Fragility_B%s_EQ_vs_Seq.png', buildingNames{i})));
    close(f2);
    
    % Median PGV Calculations
    LS1_PGV50_EQ = pgv50(a_EQ,b_EQ,LS1_pct);
    LS2_PGV50_EQ = pgv50(a_EQ,b_EQ,LS2_pct);
    LS3_PGV50_EQ = pgv50(a_EQ,b_EQ,LS3_pct);
    if hasSeq
        LS1_PGV50_SQ = pgv50(a_SQ,b_SQ,LS1_pct);
        LS2_PGV50_SQ = pgv50(a_SQ,b_SQ,LS2_pct);
        LS3_PGV50_SQ = pgv50(a_SQ,b_SQ,LS3_pct);
    else
        LS1_PGV50_SQ = NaN; LS2_PGV50_SQ = NaN; LS3_PGV50_SQ = NaN;
    end
    
    % Save Results
    Summary = [Summary; {string(bLabel),'EQ',C_EQ,k_EQ,R2_EQ,s2_EQ,beta_fs_EQ, LS1_pct, LS2_pct, LS3_pct}];
    RegStats = [RegStats; {string(bLabel),'EQ',round(s2_EQ,4),round(beta_fs_EQ,4),round(R2_EQ,4),round(double(n_EQ),4)}];
    if hasSeq
        Summary = [Summary; {string(bLabel),'Seq',C_SQ,k_SQ,R2_SQ,s2_SQ,beta_fs_SQ, LS1_pct, LS2_pct, LS3_pct}];
        RegStats = [RegStats; {string(bLabel),'Seq',round(s2_SQ,4),round(beta_fs_SQ,4),round(R2_SQ,4),round(double(n_SQ),4)}];
    end
    
    LS_names = ["LS1","LS2","LS3"];
    PGV_EQ_values = [LS1_PGV50_EQ, LS2_PGV50_EQ, LS3_PGV50_EQ];
    PGV_SQ_values = [LS1_PGV50_SQ, LS2_PGV50_SQ, LS3_PGV50_SQ];
    for k = 1:3
        MedianSummary = [MedianSummary; {string(bLabel), LS_names(k), PGV_EQ_values(k), PGV_SQ_values(k)}];
    end
end

%% --- 4. Data Export ---

% Round the Values
MedianSummary.PGV50_EQ  = round(MedianSummary.PGV50_EQ, 4);
MedianSummary.PGV50_Seq = round(MedianSummary.PGV50_Seq, 4);
Summary{:,3:end} = round(Summary{:,3:end}, 4);
RegStats{:, 3:end} = round(RegStats{:, 3:end}, 4);

% Write to CSV files
writetable(Summary, 'Fragility_Summary_MAIN.csv');
disp('OK: Full summary table (params + LS) written to Fragility_Summary.csv.');
writetable(MedianSummary, 'Median_PGV_Table.csv');
disp('Median PGV values saved to Median_PGV_Table.csv');
writetable(RegStats,'Regression_Stats.csv');
disp('Regression statistics (s², β_fs, R²) saved to Regression_Stats.csv');

%% --- 5. Group Statistics (R², beta_fs vb.) ---

% Define building groups according to the report
GroupA = compose("%02d", [14 15 16 17 18 19 20]);
GroupB = compose("%02d", [3 6 7 11 12 13]);
GroupC = compose("%02d", [1 2 4 5 8 9 10]);

Summary.BuildNum = str2double(extractAfter(Summary.Building, 'B'));

% Group Statistics Table
GroupSummary = table({'A';'B';'C'}, NaN(3,1), NaN(3,1), NaN(3,1), ...
    'VariableNames', {'Group','Mean_R2','Mean_s2','Mean_beta_fs'});
for g = 1:3
    switch g
        case 1; groupList_nums = str2double(GroupA);
        case 2; groupList_nums = str2double(GroupB);
        case 3; groupList_nums = str2double(GroupC);
    end
    idx = ismember(Summary.BuildNum, groupList_nums);
    GroupSummary.Mean_R2(g) = mean(Summary.R2(idx), 'omitnan');
    GroupSummary.Mean_s2(g) = mean(Summary.s2(idx), 'omitnan');
    GroupSummary.Mean_beta_fs(g) = mean(Summary.beta_fs(idx), 'omitnan');
end

GroupSummary{:,2:end} = round(GroupSummary{:,2:end}, 4);
writetable(GroupSummary, 'Group_Average_Stats.csv');
disp('Group-based average statistics saved to Group_Average_Stats.csv');

%% --- 6. Representative Group Fragility Curves (Parameter Averaging Method) ---

GroupList = {GroupA, GroupB, GroupC};
GroupNames = {'A','B','C'};
PGV_grid = linspace(1, 150, 400); 

for g = 1:3
    groupList = GroupList{g};
    groupName = GroupNames{g};
    
    % Identify EQ and Seq indices for the buildings in the group
    idx_EQ = ismember(extractAfter(Summary.Building,'B'), groupList) & Summary.Scenario == "EQ";
    idx_SQ = ismember(extractAfter(Summary.Building,'B'), groupList) & Summary.Scenario == "Seq";
    
    % Compute the average model parameters (a, b, beta_fs)
    a_EQ_mean = mean(log(Summary.C(idx_EQ)), 'omitnan');  % C = exp(a)
    b_EQ_mean = mean(Summary.k(idx_EQ), 'omitnan');
    beta_fs_EQ_mean = mean(Summary.beta_fs(idx_EQ), 'omitnan');
    
    a_SQ_mean = mean(log(Summary.C(idx_SQ)), 'omitnan');
    b_SQ_mean = mean(Summary.k(idx_SQ), 'omitnan');
    beta_fs_SQ_mean = mean(Summary.beta_fs(idx_SQ), 'omitnan');
    
    % Compute the average capacity (LS) parameters
    LS1_EQ_mean = mean(Summary.LS1_pct(idx_EQ), 'omitnan');
    LS2_EQ_mean = mean(Summary.LS2_pct(idx_EQ), 'omitnan');
    LS3_EQ_mean = mean(Summary.LS3_pct(idx_EQ), 'omitnan');
    
    LS1_SQ_mean = mean(Summary.LS1_pct(idx_SQ), 'omitnan');
    LS2_SQ_mean = mean(Summary.LS2_pct(idx_SQ), 'omitnan');
    LS3_SQ_mean = mean(Summary.LS3_pct(idx_SQ), 'omitnan');
    
    % Total uncertainty based on the averaged parameters
    beta_total_EQ = sqrt(beta_fs_EQ_mean^2 + beta_CL^2 + beta_M^2);
    beta_total_SQ = sqrt(beta_fs_SQ_mean^2 + beta_CL^2 + beta_M^2);
    
    % Compute exceedance probabilities using the Parameter Averaging method
    P_EQ_LS1 = Pfrag(PGV_grid, a_EQ_mean, b_EQ_mean, LS1_EQ_mean, beta_fs_EQ_mean, beta_CL, beta_M);
    P_EQ_LS2 = Pfrag(PGV_grid, a_EQ_mean, b_EQ_mean, LS2_EQ_mean, beta_fs_EQ_mean, beta_CL, beta_M);
    P_EQ_LS3 = Pfrag(PGV_grid, a_EQ_mean, b_EQ_mean, LS3_EQ_mean, beta_fs_EQ_mean, beta_CL, beta_M);
    
    P_SQ_LS1 = Pfrag(PGV_grid, a_SQ_mean, b_SQ_mean, LS1_SQ_mean, beta_fs_SQ_mean, beta_CL, beta_M);
    P_SQ_LS2 = Pfrag(PGV_grid, a_SQ_mean, b_SQ_mean, LS2_SQ_mean, beta_fs_SQ_mean, beta_CL, beta_M);
    P_SQ_LS3 = Pfrag(PGV_grid, a_SQ_mean, b_SQ_mean, LS3_SQ_mean, beta_fs_SQ_mean, beta_CL, beta_M);
    
    % Graphs
    f = figure('Visible','on'); ax_rep = axes('Parent', f); 
    hold(ax_rep, 'on'); box(ax_rep, 'on'); grid(ax_rep, 'on');
    plot(ax_rep, PGV_grid, P_EQ_LS1, '-',  'LineWidth',1.8, 'DisplayName','EQ–LS1');
    plot(ax_rep, PGV_grid, P_EQ_LS2, '-',  'LineWidth',1.8, 'DisplayName','EQ–LS2');
    plot(ax_rep, PGV_grid, P_EQ_LS3, '-',  'LineWidth',1.8, 'DisplayName','EQ–LS3');
    plot(ax_rep, PGV_grid, P_SQ_LS1, '--', 'LineWidth',1.8, 'DisplayName','Seq–LS1');
    plot(ax_rep, PGV_grid, P_SQ_LS2, '--', 'LineWidth',1.8, 'DisplayName','Seq–LS2');
    plot(ax_rep, PGV_grid, P_SQ_LS3, '--', 'LineWidth',1.8, 'DisplayName','Seq–LS3');
    xlabel(ax_rep, 'PGV (cm/s)');
    ylabel(ax_rep, 'P(Exceedance)');
    legend(ax_rep, 'Location','southeast');
    xlim(ax_rep, [0 150]); 
    yl = ylim(ax_rep); xl = xlim(ax_rep);
    text(ax_rep, 5, 0.95, sprintf('Grup %s', groupName), ...
        'Units','data', 'FontWeight','bold','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontSize',11, 'FontName','Times New Roman');
        
    saveas(f, sprintf('fragility_figs/Representative_Fragility_Group_%s_AvgParams.png', groupName));
    close(f);
end

%% --- 7. Local Function: powerlaw_fit_impl ---

function [a,b,beta_fs,R2_log,n,s2] = powerlaw_fit_impl(PGV,DRIFT)
    % Select clean and valid data (positive and finite)
    ok = PGV>0 & DRIFT>0 & isfinite(PGV) & isfinite(DRIFT);
    PGV = PGV(ok); DRIFT = DRIFT(ok);
    n   = numel(DRIFT);
    
    % Check whether there are sufficient points for fitting (minimum of 3)
    if n < 3
        warning('Insufficient data points for log–log fitting (n=%d). Returning NaN.', n);
        a=NaN; b=NaN; beta_fs=NaN; R2_log=NaN; s2=NaN; return;
    end
    
    % Transform data into logarithmic domain
    x = log(PGV(:));
    y = log(DRIFT(:));
    
    % Perform linear regression: y = p(1) + p(2)*x  (i.e., ln(D) = a + b·ln(PGV))
    X = [ones(size(x)) x];
    p = X\y; 
    a = p(1); 
    b = p(2);
    
    % Compute regression statistics
    yhat = X*p;                 % Model predictions
    res  = y - yhat;            % Residuals
    SSE = sum(res.^2);          % Sum os squared errors
    ybar = mean(y);             % mean of y
    SST  = sum((y - ybar).^2);  % Total sum of squares
    
    R2_log = 1 - SSE/SST; % R-squared (in log domain)
    
    % Mean squared error (s2) and beta_fs
    s2 = SSE / max(n-2, 1); % Degrees of freedom = n - 2
    beta_fs = sqrt(log(1 + s2)); 
end