%% Computes the joint density of $\{Y_t\}_{t=1}^T$ as in Proposition 3 of Windle & Carvalho(2014).
%
% $$Y_t = r_t\cdot r_t'.$$
%
% $$p(\{Y_t\}_{t=1}^T|D_0) = \prod_{t=1}^T p(Y_t|D_{t-1}).$$
%
% $$P(Y_t|D_{t-1})=\pi^{\frac{mk-k^2}{2}}\frac{\Gamma(0.5\nu)|L_t|^{(k-m-1)*0.5}|C_t|^{0.5n}}{\Gamma(n0.5)\Gamma(0.5k)|C_t+Y_t|^{0.5\nu}}.$$
%
% $$\nu=n+k;\,C_t=\lambda\Sigma_{t-1};\, \Sigma_t=\lambda\Sigma_{t-1}+Y_t.$$
%
% 

function [sum_ln_pYt]=MarginalWish_n(ut,T,M,cSig0_n,n,k)
% MarginalWish_n(ut,T,M,cSig0_draw,nugrid,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG

%cSig0_n = cSig0_draw;
%n = nugrid;
%k = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%


% Some useful constants
grid   = length(n); % size of the grid on n
lambda = n./(n+1);  % uses Uhlig(1997) restriction presented in Appendix 3.

% Evaluate constants outside the loop
n  = repmat(n,M,1);
nu = n+k;
dimgrid = (1:M)'*ones(1,grid);
%% Multivariate Gamma Function
const = T*sum(gammaln(nu*0.5+(1-dimgrid)*0.5)-gammaln(0.5*n+(1-dimgrid)*0.5));% repeats in all p(Y_t|D_{t-1}).

% Constant values used in the loop.
sqrtlam = sqrt(lambda);  
lndetS0 = 2*sum(log(diag(cSig0_n)));
lnlam   = log(lambda);
%% Forward Filter
ln_pYt = NaN(T,grid);
parfor i=1:grid
% for i=1:grid
    lndetSig = lndetS0; % resets lndetSig to its initial value.
    cSig = cSig0_n;     % resets Sig to its initial value.
    for t=1:T
        lndetCt  = M*lnlam(i)+lndetSig; % log determinant of Ct.
        cSig     = cholupdate(sqrtlam(i)*cSig,ut(:,t));  % computes the cholesky decomposition of Sigma_t.
        lndetSig = 2 * sum(log(diag(cSig))); % log determinant of Sigma_t.
        %% n(1,i) nu(1,i)
        ln_pYt(t,i) = (n(1,i)*0.5)*lndetCt-(nu(1,i)*0.5)*lndetSig; % computes parts of the likelihood that depend on lambda and n.        
    end
end
    
sum_ln_pYt = const+sum(ln_pYt);    
    
    




