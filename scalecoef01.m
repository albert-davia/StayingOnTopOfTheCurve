%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = scalecoef01(kappav)
% Calculate the coefficients from the cascade structure, Laurent's alpha
% kappav is a vector of kappas with scaling properties
% Liuren Wu, liurenwu@gmail.com
% April, 2009 and after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=length(kappav);
A=eye(nx); A(1,1)=1/kappav(1);
for j=1:nx-1
    for  i=1:j
        A(i,j+1)=-kappav(j)*A(i,j)/(kappav(i)-kappav(j+1));
    end
    A(j+1,j+1)=(kappav(j)/kappav(j+1))*(A(1:j,j)'*(kappav(1:j)./(kappav(1:j)-kappav(j+1))));
end