%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nlnL, lnLv, predErrv, xEstv, yFitv] = ratelikefunlf(par, data, hfun, filter, termModel, hfunpar)
% RATELIKEFUNLF - Log-likelihood function for term structure models
%
% Inputs:
%   par       - Parameters to be optimized
%   data      - Observations (nobs x ny matrix of yields)
%   hfun      - Measurement equation function handle
%   filter    - Filtering technique (e.g., 'ukf_lfnlh')
%   termModel - Term structure model (e.g., 'CANFCPv2')
%   hfunpar   - Structure with model parameters including:
%               .modelflag  - Model identifier string
%               .output_dir - Output directory (optional, default: 'output')
%
% Outputs:
%   nlnL     - Negative log-likelihood (to be minimized)
%   lnLv     - Log-likelihood values by time period
%   predErrv - Prediction errors
%   xEstv    - Extracted latent states
%   yFitv    - Fitted observations
%
% Original: Liuren Wu, July 2002, revised April 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get output directory from hfunpar or use default
if isfield(hfunpar, 'output_dir')
    output_dir = hfunpar.output_dir;
else
    output_dir = 'output';
end

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Initialize model
[ffunpar, hfunpar, xEst, PEst, Q, R] = feval(termModel, par, hfunpar);

A = ffunpar.A;
F = ffunpar.Phi;

if isfinite(det(PEst))
    nx = length(xEst(:));
    [nobs, ny] = size(data);
    
    lnLv = zeros(nobs, 1);
    predErrv = zeros(nobs, ny);
    yFitv = predErrv;
    xEstv = zeros(nobs, nx);
    
    for t = 1:nobs
        y = data(t, :)';
        ind = find(isfinite(y));
        R2 = R(ind, ind);  % Handle missing data
        
        [xEst, PEst, xPred, yPred, predErr, yVar] = feval(filter, y(ind), xEst, PEst, Q, R2, t, ind, A, F, hfun, hfunpar);
        dp = det(yVar);
        
        if isfinite(dp)
            % Compute likelihood
            lnLv(t) = -(log(dp) + predErr' * pinv(yVar) * predErr) / 2;
            
            % Store outputs
            predErrv(t, ind) = predErr';
            xEstv(t, :) = xEst';
            yFitv(t, ind) = feval(hfun, xEst, t, ind, hfunpar)';
        else
            lnLv = -1e10;
            disp('Bad parameter: non-finite determinant');
            break;
        end
    end
    
    nlnL = -mean(lnLv(4:end));
    modelflag = hfunpar.modelflag;
    
    % Save best results
    nlnfile = fullfile(output_dir, ['nln_', modelflag, '.txt']);
    if exist(nlnfile, 'file')
        nln0 = load(nlnfile);
        if nlnL < nln0
            save(nlnfile, 'nlnL', '-ascii', '-double');
            save(fullfile(output_dir, ['par_', modelflag, '.txt']), 'par', '-ascii', '-double');
            save(fullfile(output_dir, ['X_', modelflag, '.txt']), 'xEstv', '-ascii', '-double');
            save(fullfile(output_dir, ['yFit_', modelflag, '.txt']), 'yFitv', '-ascii', '-double');
        end
    else
        save(nlnfile, 'nlnL', '-ascii', '-double');
        save(fullfile(output_dir, ['par_', modelflag, '.txt']), 'par', '-ascii', '-double');
        save(fullfile(output_dir, ['X_', modelflag, '.txt']), 'xEstv', '-ascii', '-double');
        save(fullfile(output_dir, ['yFit_', modelflag, '.txt']), 'yFitv', '-ascii', '-double');
    end
else
    disp('Bad parameter: non-finite covariance');
    nlnL = 1e10;
    lnLv = [];
    predErrv = [];
    xEstv = [];
    yFitv = [];
end

end
