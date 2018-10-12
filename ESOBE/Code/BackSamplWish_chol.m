%% Forward Filtering and Backward Sampling in State-Space Models for Symmetric Positive-Definite Matrices
%
% $$r_t\in{R}^{m}, \, t=1\ldots T,$$
%
% $$r_t\sim N(0,X_t^{-1}), $$
%
% $$Y_t = r_tr_t^{\prime},$$
%
% $$Y_t\sim W_m(k,(kX_t)^{-1}),$$
%
% $$X_t=T_{t-1}^{\prime}\Psi_tT_{t-1}/\lambda,\qquad \Psi_t\sim \beta_m(\frac{n}{2},\frac{k}{2}),$$
%
% $$T_{t-1}=upper\,\, chol\,\, of\,\, X_{t-1}.$$
%
% See Windle & Carvalho (2014).
%

function [X,Sig]=BackSamplWish_chol(r_t,T,m,cSig0,lambda,n,k)

% % Fixa semente:
% s = RandStream('mt19937ar','Seed',1234); % o primeiro argumento fixa o algoritmo gerador de números pseudo-aleatórios, o segundo indica que iremos escolher o valor da semente e o terceiro é o próprio valor da semente.
% RandStream.setGlobalStream(s); % usa a semente s como semente para o gerador de números aleatórios de todas as distribuições.
%%%% DEBUG%%%%%%%%%%%%%%
%BackSamplWish_chol(ut,T,M,cSig0_draw,lam_Hdraw,nu_Hdraw,1);
% r_t = ut;
% m = M;
% cSig0 = cSig0_draw;
% lambda = lam_Hdraw;
% n = nu_Hdraw;
% k = 1;
% END OF DEBUG%%%%%%%%%%%%%%%%%%%%%
Sig        = NaN(m,m,T+1);
Sig(:,:,1) = cSig0;
sqrtlam    = sqrt(lambda);
sqrtk      = sqrt(k);
Im         = eye(m);
%% Forward Filter
for t=2:T+1
    cCt = sqrtlam*Sig(:,:,t-1);
    Sig(:,:,t)= cholupdate(cCt,r_t(:,t-1));
end

X = Sig;
cS  = sqrtk*Sig(:,:,end);
icS = cS\Im;
X(:,:,end) = (wishrnd(icS*icS',k+n));

%% Backward Sampling
for t=T:-1:2
    cS  = sqrtk*Sig(:,:,t);
    icS = cS\Im;
    iS = icS*icS';
    X(:,:,t) = (lambda*X(:,:,t+1) + wishrnd(iS,k));
end
X = X(:,:,2:end);


    
    