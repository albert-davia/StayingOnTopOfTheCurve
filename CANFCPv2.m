%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ffunpar, hfunpar,xEst,PEst,Q,R] = CANFCPv2(par,hfunpar)
% 2c: the model is about right;scaling in Q, constant market price 
% 2cr: constrain kappanp>0
% 2cr3: scale in p, constant g0+g1*X
% 2cr4: scale in P, constant g0+g1(x2t-x1t); estimates are better
% 2cr6: g0 only.
% Use analytical formula
% Liuren Wu, liurenwu@gmail.com
% April, 2009 and after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
swapmat=hfunpar.swapmat;
libormat=hfunpar.libormat;
dt=hfunpar.dt;
ny=hfunpar.ny;
nx=hfunpar.nx;

npar=length(par(:));
par=reshape(par, npar,1);
eyex=eye(nx);
epar=exp(par);
kappar=epar(1);
sigmar=epar(2);
thetarp=epar(3);
b=exp(epar(4));
gamma0=par(5);
R=epar(6)*eye(ny);
gamma1=par(7:6+nx);

gamma0v=gamma0*sigmar;
kappav=zeros(nx,1); kappav(nx)=kappar;
Kappa=zeros(nx,nx); Kappa(nx,nx)=kappar;
for n=nx-1:-1:1
    kappav(n)=kappav(n+1)*b; 
    Kappa(n,n:n+1)=[kappav(n),-kappav(n)];
end
Kappatheta=zeros(nx,1);Kappatheta(nx)=kappar*thetarp; % sembl
Kappas=Kappa-repmat(sigmar*gamma1(:)',nx,1);
Kappathetas=Kappatheta-gamma0v;
theta=Kappa\Kappatheta;  % Point de départ qui ne dépend que des paramètres. Semble correpsondre à choisir X0 poru annuler k(theta - X0) dans l'équation (2)
br=zeros(nx,1);br(1)=1;

SS2=sigmar^2*eye(nx);
Q=SS2*dt;
Phi = expm(-Kappa*dt); 
ffunpar.Phi=Phi;
A=(eyex-Phi)*theta;
ffunpar.A=A;
xEst = theta;
PEst = Q;
ffunpar.Q=Q;

% avec A, phi, racine(Q), on peut calculer le prochain Xt

%%%%Measurement
%%%%Measurement
h=2;hfunpar.h=h;
matv=[1/h:1/h:max(swapmat)]; nm=length(matv); at_swap=zeros(nm,1); bt_swap=zeros(nm,nx); 
[U,D]=eig(Kappas); invU=inv(U);d=diag(D);
epd=exp(d); lepdv=log(epd*epd');vvs=invU*invU'; 
invKappas=inv(Kappas); invKappaspbr=invKappas'*br;
for k=1:nm
    t=matv(k);
    epd=exp(-d*t); epdv=epd*epd'; 
    VVs=U*(vvs.*(1-epdv)./lepdv)*U';
    IKappas=eyex-expm(-Kappas*t); 
    IKappasp=eyex-expm(-Kappas'*t); 
    
    btv=IKappasp*invKappaspbr;
    atv=  Kappathetas'*invKappaspbr*t -Kappathetas'*invKappas'*IKappasp*invKappaspbr   ...
        -0.5*br'*invKappas*SS2*invKappas'*br*t  ...
        +0.5*br'*invKappas*SS2*invKappas'*IKappasp*invKappaspbr ...
        +0.5*br'*invKappas*SS2*invKappas*IKappas*invKappaspbr ...
        -0.5*br'*invKappas*SS2*VVs*invKappaspbr;
    bt_swap(k,:)=real(btv');    at_swap(k,1)=real(atv);
end

hfunpar.at_swap=at_swap;
hfunpar.bt_swap=bt_swap;
hfunpar.swapmatv=matv';

nlibor=length(libormat);

at_libor=zeros(nlibor,1); bt_libor=zeros(nlibor,nx); 
for k=1:nlibor
    t=libormat(k);
    epd=exp(-d*t); epdv=epd*epd';
    VVs=U*(vvs.*(1-epdv)./lepdv)*U';
    IKappas=eyex-expm(-Kappas*t);
    IKappasp=eyex-expm(-Kappas'*t);

    btv=IKappasp*invKappaspbr;
    atv=  Kappathetas'*invKappaspbr*t -Kappathetas'*invKappas'*IKappasp*invKappaspbr   ...
        -0.5*br'*invKappas*SS2*invKappas'*br*t  ...
        +0.5*br'*invKappas*SS2*invKappas'*IKappasp*invKappaspbr ...
        +0.5*br'*invKappas*SS2*invKappas*IKappas*invKappaspbr ...
        -0.5*br'*invKappas*SS2*VVs*invKappaspbr;
    bt_libor(k,:)=real(btv');    at_libor(k,1)=real(atv);
end

hfunpar.at_libor=at_libor;
hfunpar.bt_libor=bt_libor;