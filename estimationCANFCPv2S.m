%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cascade Affine Term Structure Model Estimation
%  Based on: Calvet, Fisher, Wu (2018) - "Staying on Top of the Curve"
%
%  Scale in kappa_Q, constant sigma, add the same market price g0+g1X to all risk sources.
%  Original code: Liuren Wu, liurenwu@gmail.com, April 2009
%  Adapted for French sovereign debt analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; format compact; format short;

%% =========================== CONFIGURATION ===========================
% All paths are relative to the current working directory.
% Users should place their data in a 'data/' subfolder.

% --- Path Configuration ---
DATA_DIR = 'data';                    % Directory containing input data
OUTPUT_DIR = 'output';                % Directory for output files
FIGURES_DIR = 'figures';              % Directory for generated figures

% Create output directories if they don't exist
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end
if ~exist(FIGURES_DIR, 'dir'), mkdir(FIGURES_DIR); end

% --- Data Configuration ---
% Set to 1 for French sovereign debt data, 0 for original US LIBOR/swap data
dataDette = 1;

% Data file names (should be placed in DATA_DIR)
DATA_FILE_DETTE = 'nusrates_dette_cleaned.mat';   % French debt data
DATA_FILE_US = 'nusrates.mat';                     % US LIBOR/swap data

% --- Estimation Settings ---
estimation = 1;         % 1 = run optimization, 0 = use existing parameters
unc = 1;                % 1 = unconstrained optimization (fminunc), 0 = fminsearch
stderror = 1;           % 1 = compute standard errors
AttemptNumber = '31';   % Run identifier for saving outputs

% --- Visualization Settings ---
gammavplot = 0;         % 1 = plot gamma coefficients
draw = 0;               % 1 = draw yield curve plots

% --- Prediction Settings ---
prediction_after = 1;   % 1 = forecast future yields
prediction_before = 1;  % 1 = backtest model predictions
info = 4;               % Recalibration frequency (weeks)
Tpred_before = 357;     % Backtesting horizon (weeks)
Tpred_after = 357;      % Forecasting horizon (weeks)
plot_maturity = 7;      % Maturity index for prediction plots

% --- Optimization Parameters ---
tol = 1e-4;
nit = 50000;
fopt = optimset('Display', 'iter', 'MaxIter', nit, 'MaxFunEvals', nit, 'TolX', tol, 'TolFun', tol);

% --- Model Settings ---
filter = 'ukf_lfnlh';
likefun = 'ratelikefunlf';
nx = 10;  % Number of latent factors

%% =========================== DATA LOADING ===========================

if dataDette
    % Load French sovereign debt data
    datafile = fullfile(DATA_DIR, DATA_FILE_DETTE);
    if ~exist(datafile, 'file')
        error('Data file not found: %s\nPlease place your data file in the ''%s'' directory.', datafile, DATA_DIR);
    end
    load(datafile, 'rates', 'swapmat', 'mdate', '-mat');

    cdate = [mdate(end):mdate(1)]';
    wdate = cdate(weekday(cdate) == 4);
    dt = 1/52;
    rates = interp1(mdate, rates, wdate);
    mat = [swapmat];
    [T, ny] = size(rates);
    datevec([wdate(1); wdate(end)]);
    lastdate = datestr(wdate(end), 1);

    % Model configuration
    termModel = 'CANFCPv2';
    hfun = 'liborswap';
    hfunpar.dt = dt;
    hfunpar.ny = ny;
    hfunpar.swapmat = swapmat;
    hfunpar.libormat = [];
    hfunpar.nx = nx;
    modelflag = [termModel, '_FS', num2str(nx)];
    hfunpar.modelflag = modelflag;
    hfunpar.output_dir = OUTPUT_DIR;  % Pass output directory to likelihood function

    % Load or initialize parameters
    parfile = fullfile(OUTPUT_DIR, ['par_', modelflag, '_', AttemptNumber, '.txt']);
    if exist(parfile, 'file')
        par = load(parfile);
    else
        % Default initial parameters for French debt
        par = [-2.5556001063807963e+00
            -2.6384878563863501e-01
            -8.5619137343435925e-01
            -1.0498554832404401e+00
            -1.2585177801722682e+00
            -7.8172780580145824e+00
            -6.2807525243175635e-02
            1.5485760897156067e-01
            4.0269433282983880e-02
            -3.0223589931092056e-02
            -9.0114087099664231e-02
            -5.5738327012341410e-03
            -4.9786674741841633e-02
            1.4841023869127912e-02
            -6.9679602544003147e-03
            -1.1378875655991154e-02]';
    end
else
    % Load US LIBOR/swap data
    datafile = fullfile(DATA_DIR, DATA_FILE_US);
    if ~exist(datafile, 'file')
        error('Data file not found: %s\nPlease place your data file in the ''%s'' directory.', datafile, DATA_DIR);
    end
    load(datafile, 'rates', 'mat', 'swapmat', 'libormat', 'mdate', '-mat');

    cdate = [mdate(1):mdate(end)]';
    wdate = cdate(weekday(cdate) == 4);
    dt = 1/52;
    rates = interp1(mdate, rates(:, [4, 7:end]), wdate);
    libormat = 6;
    mat = [6/12; swapmat];
    [T, ny] = size(rates);
    datevec([wdate(1); wdate(end)]);
    lastdate = datestr(wdate(end), 1);

    % Model configuration
    termModel = 'CANFCPv2';
    hfun = 'liborswap';
    hfunpar.dt = dt;
    hfunpar.ny = ny;
    hfunpar.swapmat = swapmat;
    hfunpar.libormat = libormat' / 12;
    hfunpar.nx = nx;
    modelflag = [termModel, '_FS', num2str(nx)];
    hfunpar.modelflag = modelflag;
    hfunpar.output_dir = OUTPUT_DIR;

    % Default initial parameters for US data
    par = [-2.9702, -4.2022, -9.9750, -0.5565, -0.1706, -9.6040, -0.2710, 0.1458, 0.0338, 0.0892, 3.7794, -0.0607, -0.2011, -0.9491, -1.6781, 0.0099]';
end

fprintf('Data loaded: %d observations, %d maturities\n', T, ny);
fprintf('Sample period: %s to %s\n', datestr(wdate(1)), datestr(wdate(end)));

%% =========================== PARAMETER EXTRACTION ===========================

epar = exp(par);
kappar = epar(1);
sigmar = epar(2);
thetarp = epar(3);
b = exp(epar(4));
gamma0 = par(5);
R = epar(6) * eye(ny);
gamma1 = par(7:6+nx);
gamma0v = gamma0 * sigmar;

%% =========================== ESTIMATION ===========================

t0 = clock;
[loglike, likeliv, predErr, mu_dd, y_dd] = feval(likefun, par, rates, hfun, filter, termModel, hfunpar);
fprintf('Initial log-likelihood: %.4f\n', loglike);
runtime = etime(clock, t0);

if estimation
    fprintf('Starting optimization...\n');
    if unc
        par = fminunc(likefun, par, fopt, rates, hfun, filter, termModel, hfunpar);
    else
        par = fminsearch(likefun, par, fopt, rates, hfun, filter, termModel, hfunpar);
    end

    [loglike, likeliv, predErr, mu_dd, y_dd] = feval(likefun, par, rates, hfun, filter, termModel, hfunpar);
    fprintf('Final log-likelihood: %.4f\n', loglike);

    % Save results
    save(fullfile(OUTPUT_DIR, ['par_', modelflag, '_', AttemptNumber, '.txt']), 'par', '-ascii', '-double');
    save(fullfile(OUTPUT_DIR, ['mu_dd_', modelflag, '_', AttemptNumber, '.txt']), 'mu_dd', '-ascii', '-double');
    save(fullfile(OUTPUT_DIR, ['y_dd_', modelflag, '_', AttemptNumber, '.txt']), 'y_dd', '-ascii', '-double');
    fprintf('Results saved to %s/\n', OUTPUT_DIR);
end

%% =========================== STANDARD ERRORS ===========================

if stderror
    fprintf('Computing standard errors...\n');
    npar = length(par);
    h = zeros(T, npar);
    stepsizei = 1e-4;

    for i = 1:npar
        bb1 = par; bb2 = par;
        bb1(i) = par(i) - stepsizei/2;
        bb2(i) = par(i) + stepsizei/2;
        [~, x1] = feval(likefun, bb1, rates, hfun, filter, termModel, hfunpar);
        [~, x2] = feval(likefun, bb2, rates, hfun, filter, termModel, hfunpar);
        h(:, i) = (x2 - x1) ./ stepsizei;
    end

    [Q, R_qr] = qr(h, 0);
    Rinv = pinv(R_qr);
    avar = Rinv * Rinv';
    stdpar = sqrt(diag(avar));

    % Extract parameters for display
    epar = exp(par);
    kappar = epar(1);
    sigmar = epar(2);
    thetarp = epar(3);
    b = exp(epar(4));
    gamma0 = par(5);
    R = epar(6);
    gamma1 = par(7:6+nx);

    % Format for display
    if size(par, 1) == 1
        parpr = exp(par)';
    else
        parpr = exp(par);
    end

    stdpp = parpr .* stdpar;
    parpr(4) = b;
    stdpp(4) = b * epar(4) * stdpar(4);
    ind = [5, 7:6+nx];
    parpr(ind) = par(ind);
    stdpp(ind) = stdpar(ind);

    table = [parpr, stdpp];
    ind = [1:5, 7:6+nx]';
    fprintf('\n=== Parameter Estimates ===\n');
    fprintf('   &  %7.4f  &  ( %7.4f )   \\\\ \n', table(ind, :)');
end

%% =========================== GAMMA PLOTS ===========================

if gammavplot
    epar = exp(par);
    kappar = epar(1);
    b = exp(epar(4));
    gamma1 = par(7:6+nx);

    kappav = zeros(nx, 1);
    kappav(nx) = kappar;
    for n = nx-1:-1:1
        kappav(n) = kappav(n+1) * b;
    end
    fv = [nx:-1:1];

    figure(1); clf
    plot(fv, gamma1, 'o-', 'LineWidth', 2, 'MarkerSize', 10)
    xlabel('Frequency, j', 'FontSize', 16)
    ylabel('\lambda_j', 'FontSize', 16)
    grid on
    set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    saveas(gcf, fullfile(FIGURES_DIR, ['figgammav_', modelflag, '_', AttemptNumber, '.png']));

    lk = -log(kappav);
    figure(2); clf
    plot(lk, gamma1, 'o-', 'LineWidth', 2, 'MarkerSize', 10)
    ylabel('\lambda_j', 'FontSize', 16)
    xlabel('ln \kappa_j', 'FontSize', 16)
    grid on
    set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    saveas(gcf, fullfile(FIGURES_DIR, ['figlnkappavgammav_', modelflag, '_', AttemptNumber, '.png']));
end

%% =========================== YIELD CURVE PLOTS ===========================

if draw
    drawperiod = 1:T;
    [loglike, likeliv, predErr, mu_dd, y_dd] = feval(likefun, par, rates, hfun, filter, termModel, hfunpar);

    % Plot latent factors
    figure(6); clf
    columns = [1, 3, 10];
    colors = {[0 0 0.5], [0 0.5 0], [0.5 0 0]};
    for i = 1:length(columns)
        subplot(3, 1, i)
        plot(wdate(drawperiod), mu_dd(drawperiod, columns(i)), 'Color', colors{i}, 'LineWidth', 2)
        xlabel('Time', 'FontSize', 16)
        ylabel(['Factor ', num2str(columns(i))], 'FontSize', 16)
        datetick('x', 'mmmyy')
        grid on
        set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    end
    saveas(gcf, fullfile(FIGURES_DIR, ['rowstotime_', modelflag, '_', AttemptNumber, '.png']));

    % Plot yield curves
    figure(3); clf
    maturities = mat;
    if dataDette
        columns = [1, 3, 9, 10];
    else
        columns = [1, 3, 10];
    end
    colors = {[0 0 0.5], [0 0.5 0], [0.5 0 0], [0.5 0 0]};

    for i = 1:length(columns)
        subplot(length(columns), 1, i)
        plot(wdate(drawperiod), rates(drawperiod, columns(i)), 'LineWidth', 1, 'Color', colors{i})
        hold on
        plot(wdate(drawperiod), y_dd(drawperiod, columns(i)), 'r--', 'LineWidth', 1)
        hold off
        datetick('x', 'mmmyy')
        grid on
        legendLabels = sprintf('Maturity %.1f years', maturities(columns(i)));
        legend({legendLabels, ['Model ', legendLabels]}, 'Location', 'Best')
        ylabel('Yield (%)', 'FontSize', 16)
        set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    end
    xlabel('Date', 'FontSize', 16)
    saveas(gcf, fullfile(FIGURES_DIR, ['figyieldcurve_', modelflag, '_', AttemptNumber, '.png']));

    % Plot pricing errors
    figure(4); clf
    for i = 1:length(columns)
        subplot(length(columns), 1, i)
        plot(wdate(drawperiod), rates(drawperiod, columns(i)) - y_dd(drawperiod, columns(i)), 'LineWidth', 2)
        datetick('x', 'mmmyy')
        grid on
        legendLabels = sprintf('Error: Maturity %.1f years', maturities(columns(i)));
        legend(legendLabels, 'Location', 'Best')
        set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    end
    saveas(gcf, fullfile(FIGURES_DIR, ['figyieldcurveerror_', modelflag, '_', AttemptNumber, '.png']));

    % Display RMSE
    fprintf('\n=== Root-Mean-Squared Pricing Errors ===\n');
    for i = 1:min(10, ny)
        rmse = sqrt(mean((rates(drawperiod, i) - y_dd(drawperiod, i)).^2));
        fprintf('Maturity %.1f years: %.4f\n', maturities(i), rmse);
    end
end

%% =========================== FORWARD PREDICTION ===========================

T = Tpred_after;

if prediction_after
    [ffunpar, hfunpar, xEst, PEst, Q, R] = feval(termModel, par, hfunpar);

    PredY = zeros(T, ny);
    PredX = zeros(T, nx);
    xPred = mu_dd(end, :)';

    npar = length(par(:));
    par = reshape(par, npar, 1);
    epar = exp(par);
    R = epar(6) * eye(ny);
    Q = ffunpar.Q;
    PEst = Q;
    F = ffunpar.Phi;

    for i = 1:T
        A = ffunpar.A;
        phi = ffunpar.Phi;

        xPred = A + phi * xPred;
        xPredVar = F * PEst * F' + Q;
        [xPredSigmaPts, wSigmaPts, wSigmaPtsc, nsp] = SigmaPoints(xPred, xPredVar);
        yPredSigmaPts = feval(hfun, xPredSigmaPts, 0, 0, hfunpar);
        yPred = yPredSigmaPts * wSigmaPts';

        PredY(i, :) = yPred;
        PredX(i, :) = xPred;

        wSigmaPts_ymat = repmat(wSigmaPtsc, ny, 1);
        exSigmaPt = xPredSigmaPts - repmat(xPred, 1, nsp);
        eySigmaPt = yPredSigmaPts - repmat(yPred, 1, nsp);
        yPredVar = (wSigmaPts_ymat .* eySigmaPt) * eySigmaPt' + R;
        xyPredVar = exSigmaPt * (wSigmaPts_ymat .* eySigmaPt)';

        K = xyPredVar / yPredVar;
        PEst = xPredVar - K * yPredVar * K';
    end

    figure(5); clf
    plot(wdate(end-T:end), rates(end-T:end, plot_maturity), 'LineWidth', 2, 'DisplayName', 'Actual Yield')
    hold on
    plot(wdate(end-T:end), y_dd(end-T:end, plot_maturity), 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted Yield')
    future_dates = wdate(end) + (1:T)' * 7;
    plot(future_dates, PredY(:, plot_maturity), 'g-', 'LineWidth', 2, 'DisplayName', 'Forecast')
    hold off
    datetick('x', 'mmmyy')
    grid on
    legend('show', 'Location', 'Best')
    ylabel('Yield (%)', 'FontSize', 16)
    xlabel('Date', 'FontSize', 16)
    title(['Yield Forecast (', num2str(T), ' weeks ahead)'], 'FontSize', 16)
    set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    saveas(gcf, fullfile(FIGURES_DIR, ['prediction_forward_', modelflag, '_', AttemptNumber, '.png']));
end

%% =========================== BACKTESTING ===========================

T = Tpred_before;

if prediction_before
    [ffunpar, hfunpar, xEst, PEst, Q, R] = feval(termModel, par, hfunpar);

    PredY = zeros(T, ny);
    PredX = zeros(T, nx);
    xEst = mu_dd(end-T, :)';

    npar = length(par(:));
    par = reshape(par, npar, 1);
    epar = exp(par);
    R = epar(6) * eye(ny);
    A = ffunpar.A;
    F = ffunpar.Phi;
    Q = ffunpar.Q;

    for i = 1:T
        xPred = A + F * xEst;
        xPredVar = F * PEst * F' + Q;

        [xPredSigmaPts, wSigmaPts, wSigmaPtsc, nsp] = SigmaPoints(xPred, xPredVar);
        yPredSigmaPts = feval(hfun, xPredSigmaPts, 0, 0, hfunpar);
        yPred = yPredSigmaPts * wSigmaPts';

        wSigmaPts_ymat = repmat(wSigmaPtsc, ny, 1);
        exSigmaPt = xPredSigmaPts - repmat(xPred, 1, nsp);
        eySigmaPt = yPredSigmaPts - repmat(yPred, 1, nsp);
        yPredVar = (wSigmaPts_ymat .* eySigmaPt) * eySigmaPt' + R;
        xyPredVar = exSigmaPt * (wSigmaPts_ymat .* eySigmaPt)';

        PredY(i, :) = yPred';
        PredX(i, :) = xPred;

        if mod(i, info) == 0
            K = xyPredVar / yPredVar;
            innovation = rates(size(rates, 1) - T + i, :)' - yPred;
            xEst = xPred + K * innovation;
            PEst = xPredVar - K * yPredVar * K';
        else
            K = xyPredVar / yPredVar;
            xEst = xPred;
            PEst = xPredVar - K * yPredVar * K';
        end
    end

    % Calculate MAE
    fprintf('\n=== Mean Absolute Errors (Backtesting) ===\n');
    for i = 1:min(10, ny)
        MAE = mean(abs(PredY(:, i) - rates(end-T+1:end, i)));
        fprintf('Maturity %d: %.4f\n', i, MAE);
    end

    figure(6); clf
    subplot(2, 1, 1)
    plot(wdate(end-T:end), rates(end-T:end, plot_maturity), 'LineWidth', 2, 'DisplayName', 'Actual')
    hold on
    plot(wdate(end-T:end), y_dd(end-T:end, plot_maturity), 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted')
    hold off
    datetick('x', 'mmmyy')
    grid on
    legend('show', 'Location', 'Best')
    ylabel('Yield (%)', 'FontSize', 16)
    title(['In-Sample Fit (Last ', num2str(T), ' Weeks)'], 'FontSize', 16)
    set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)

    subplot(2, 1, 2)
    plot(wdate(end-T+1:end), PredY(:, plot_maturity), 'LineWidth', 2, 'DisplayName', 'Predicted')
    hold on
    plot(wdate(end-T+1:end), rates(end-T+1:end, plot_maturity), 'r--', 'LineWidth', 2, 'DisplayName', 'Actual')
    hold off
    datetick('x', 'mmmyy')
    grid on
    legend('show', 'Location', 'Best')
    ylabel('Yield (%)', 'FontSize', 16)
    xlabel('Date', 'FontSize', 16)
    title(['Backtest: Recalibrated Every ', num2str(info), ' Weeks'], 'FontSize', 16)
    set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 16)
    saveas(gcf, fullfile(FIGURES_DIR, ['backtest_', modelflag, '_', AttemptNumber, '.png']));
end

fprintf('\n=== Estimation Complete ===\n');
