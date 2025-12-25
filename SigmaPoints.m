function [xPts, wPts,wPtsc, nPts] = SigmaPoints(x,P)

% Outputs:
%        xPts	 The sigma points
%        wPts	 The weights on the points
%	 nPts	 The number of points
%

% Number of sigma points and scaling terms
n    = length(x);
nPts = 2*n+1;            

kappa=0;alpha=.01;beta=2;
alpha2=alpha^2;
delta=alpha2*(n+kappa)-n;

% Calculate matrix square root of weighted covariance matrix
[R,p]=chol(P);
Psqrtm = sqrt(n+delta)*R';

% Array of the sigma points
xPts=[zeros(n,1) -Psqrtm Psqrtm];
% Add mean back in
xPts = xPts + repmat(x,1,nPts);  

% Array of the weights for each sigma point
wPts=[delta, 0.5*ones(1,nPts-1)]/(n+delta);
wPtsc=wPts;
wPtsc(1)=wPtsc(1)+(1-alpha2+beta);

%wPtsc=[delta/(n+delta)+(1-alpha^2+beta), wP];