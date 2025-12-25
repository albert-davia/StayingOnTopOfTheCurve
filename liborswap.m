%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = liborswap(x,t,ind,hfunpar)
% Measurement equation on libor and swap rates
% libor rates are based on actual/360; libor=100*(exp(yt)-1)/t
% swap rates are based on half year tenor; swap=200*(1-disc(mat))/cumsum(disc)
%  Liuren Wu, liurenwu@gmail.com, April 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nx nsigma]=size(x);

at_swap=   hfunpar.at_swap;
bt_swap=   hfunpar.bt_swap;
swapmat=   hfunpar.swapmat;
h =hfunpar.h; %number of payments per year

% dis=exp(-repmat(at_swap,1,nsigma)-bt_swap*x);
% swr=h*100*(1-dis)./cumsum(dis);
% y_swap=swr(swapmat*h,:);


negative_lnP = repmat(at_swap,1,nsigma)+bt_swap*x;
y_swap = 100*negative_lnP(h*swapmat,:)./repmat(swapmat,1,nsigma)  ;


libormat= hfunpar.libormat;
at_libor=   hfunpar.at_libor;
bt_libor=   hfunpar.bt_libor;

if size(libormat,1)==0
    y=[y_swap];
else
    y_libor=100*(exp(repmat(at_libor,1,nsigma)+bt_libor*x) -1)./repmat(libormat,1,nsigma);
    y=[y_libor;y_swap];
    y=y(ind,:);
end





