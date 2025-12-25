%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [xEst,PEst,xPred,yPred,innovation,yPredVar,K] = ukf_lfnlh(y,xEst,PEst,Q,R,t,ind,A,F, hfun, hfunpar)
% TITLE    :  UNSCENTED KALMAN FILTER for additive system and measurement
% errors with linear state but nonlinear measurement
% state:   x=A+FX + e_t, cov(e)= Q
% measure: y=h(x_t,t) + e_t, cov(e)= R
% INPUTS   :  - y                : observation at k+1  
%             - xEst             : state mean estimate at time k  
%             - PEst             : state covariance at time k
%             -c1                : scaling parameters. Set to default values.  
%             - Q                : process noise covariance at time k  
%             - R                : measurement noise covariance at k+1  
%             - t                : time index
%             - ind              : index for non-missing values
%             - A, F             : X=A+FX 
%             - hfun             : observation model function  
%             - hfunpar          : parameters (in a structure) for the observation model function  
%
% OUTPUTS  :  - xEst             : updated estimate of state mean at time k+1
%             - PEst             : updated state covariance at time k+1
%             - xPred            : prediction of state mean at time k+1
%             - zPred            : prediction of state mean at time k+1
%	           - inovation        : prediction error
%	           - yPredVar                : prediction error variance
%	           - K                : Kalman gain
% @Liuren Wu, December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xEst,PEst,xPred,yPred,innovation,yPredVar,K] = ukf_lfnlh(y,xEst,PEst,Q,R,t,ind,A,F, hfun, hfunpar)
%nx = length(xEst(:)); %states
ny = length(y(:)); %observations
%Linear prediction on X
xPred=A+F*xEst;
xPredVar=F*PEst*F'+Q;

% Nonlinear Prediction of y
[xPredSigmaPts, wSigmaPts,wSigmaPtsc, nsp] = SigmaPoints(xPred,xPredVar);
yPredSigmaPts = feval(hfun,xPredSigmaPts,t,ind,hfunpar);
yPred = yPredSigmaPts*wSigmaPts';

% Prediction of covariances
wSigmaPts_ymat = repmat(wSigmaPtsc,ny,1);
exSigmaPt = xPredSigmaPts - repmat(xPred,1,nsp);
eySigmaPt = yPredSigmaPts - repmat(yPred,1,nsp);
yPredVar  = (wSigmaPts_ymat .* eySigmaPt) * eySigmaPt' + R;
xyPredVar = exSigmaPt * (wSigmaPts_ymat .* eySigmaPt)';

%K  = xyPredVar*pinv(yPredVar);
K  = xyPredVar/(yPredVar);
innovation = y - yPred;
xEst = xPred + K*innovation;
PEst = xPredVar - K*yPredVar*K';

return