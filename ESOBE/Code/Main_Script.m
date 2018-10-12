%% Script for a Wishart TVP-VAR with homoscedastic Bs
%
% The model uses Uhlig (1997) restriction and is given by:
%
% $$y_t = B_tX_t+u_t \qquad u_t\sim N_m(0,H_t^{-1}),\\$$
%
% $$vec(B_{t}^{\prime}) = vec(B_{t-1}^{\prime})+\eta_t \qquad \eta_t\sim N_s(0,W^{-1}),\\$$
%
% $$H_{t} =\frac{1}{\lambda}\mathcal{U}(H_{t-1})^{\prime}\Psi_t\mathcal{U}(H_{t-1}),\qquad \Psi_t\sim B_m(n/2,0.5),\\$$
%
% $$\lambda=\frac{n}{n+1},$$  
%
% $$ H_{1}\sim W(n,\Sigma_0^{-1}/\lambda).$$ 
%

% Housekeeping
clearvars
format shortg
c = clock;
disp(table(c(3),c(2),c(1),c(4),c(5),c(6),'VariableNames',{'day' 'month' 'year' 'hour' 'minutes' 'seconds'}))

% Fix seed and algorithm for pseudo-random number generation
seed = 1234;
semente = RandStream('mt19937ar','Seed',seed); 
RandStream.setGlobalStream(semente); % use object s to fix seed and algorithm of the global random number generator.

%% Load data 

% Data is in this order: K/L, Inflation, GDP, Exchange rate and interest
% rate - data cleaning was made using R
load xdata2.dat
Y = xdata2;
% Year labels
load anos.dat
yearlab = anos;

% Dimensions of the data set
T  = size(Y,1); % T is the number of time periods.
M  = size(Y,2); % M is the number of time series.
p  = 2;         % p is the desired lag length.
tr = 50-p;      % size of traning sample for the prior. 

% Prepare data
% Generates matrix with p lags of the observables
for ii=1:p
    ylag(p+1:T,(M*(ii-1)+1):M*ii)=Y(p+1-ii:T-ii,1:M);
end
ylag = ylag(p+1:T,:);

s = M + p*(M^2); % number of VAR coefficients
Is = speye(s);

% Generates matrix with explanatory variables to go with vec(Y)
Z = zeros((T-p)*M,s);
for i = 1:T-p
    ztemp = eye(M); % constante
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];  
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end
% Data set for prior estimation
Zprior = Z(1:M*tr,:);
yprior = Y(1+p:tr+p,:)';
% Data set for parameter estimation
Z      = Z(M*tr+1:end,:);
y      = Y(tr+p+1:T,:)';
% New number of observations
T  = size(y,2);  % T is now T - p - tr.
    

%% Prior based on a training sample.
Bols    = (Zprior'*Zprior)\(Zprior'*yprior(:)); 
res_OLS = yprior(:)-Zprior*Bols;
uprior  = reshape(res_OLS,M,tr);
hbar  = (uprior*uprior')./(tr-M-1);
vbar  = zeros(length(Bols));
for i = 1:tr
    zhat1 = Zprior((i-1)*M+1:i*M,:);
    vbar  = vbar + zhat1'*(hbar\zhat1);
end
ic_vbar  = chol(vbar)\Is;
Bols_var = ic_vbar*ic_vbar';

% Paramenters of the prior for B0
c = Bols;
B_0_prmean = Bols;
B_0_prvar  = 4*Bols_var; % inflate OLS covariance matrix
num_Betas  = size(Bols,1);

%% Wishart conjugate prior for W obtained from training sample following Primiceri.
%
W_prdf     = tr;%
K_w        = 0.01;
iW_prscale = (0.25*B_0_prvar*W_prdf*(K_w^2));


%% Uniform prior for $\nu_H$ using equally spaced grid.
%
% $$H_1\sim W_M(\nu_H,\Sigma_0^-1/\lambda_H).$$
%
nugrid   = linspace(M-1+0.1,30,201); % equally spaced grid

% Initial value $\Sigma_0$ will be estimated based on a Wishart conjugate prior.
% Cov(A) - If A is a matrix whose columns represent random variables and whose rows represent observations, C is the covariance matrix with the corresponding column variances along the diagonal.
cSig0  = sqrt((tr-1)/(tr-M-1))*chol(cov(uprior'));
nu_H0  = 4;
Im     = speye(M);
cinvH0 = cSig0\Im;
invS0  = cinvH0*cinvH0';
cSig0_draw = cSig0;

% Set initial values for Bs, W and H.
iWdraw = 0.0001*eye(s);
Hdraw  = repmat(eye(M),1,1,T);
Btdraw     = B_0_prmean;

%% Starts the Gibbs sampler
nrep    = 5000;    % number of replications.
nburn   = 10000;    % burnin sample
Bt_M    = zeros(s,T);
Bt_S    = Bt_M;
iH_M    = zeros(M,M,T);
iH_S    = iH_M;
iHdraw  = iH_M;
Sig0_M  = iH_M;
Sig0_S  = Sig0_M;
iW_M    = zeros(s,s);
iW_S    = iW_M;
ut_M    = zeros(M,T);
ut_S    = zeros(M,T);
ut2_M   = zeros(M,T);
ut2_S   = zeros(M,T);
nuH_save = nan(nrep,1);
Bt_save  = nan(num_Betas, T+1, nrep);
Ht_save  = nan((M.^2+M)/2, T+1, nrep);
temp = nan(M,M);

% This is used for impulse response functions
nhor = 25; %horizont for impulse response function
imp00 = zeros(nrep,M,nhor);
imp08 = zeros(nrep,M,nhor);
imp17 = zeros(nrep,M,nhor);
bigj = zeros(M,M*p);
bigj(1:M,1:M) = eye(M);

%% Allow for parallelization.
if isempty(p)
    parpool 
end
for irep   = 1:nrep + nburn    
    %% Sample Bs using Carter & Kohn (1994).
    Btdraw = carter_kohn_MSV_1(y,Z,Hdraw,iWdraw,s,M,T,B_0_prmean,B_0_prvar);%,semente);
      
    
    %Btdraw2 = array2table(Btdraw);
    %writetable(Btdraw2);    
    %% Sample W from Wishart conditional posterior.
    %
    % $$W\sim W_s(T+\tau_0,(SSE+Q_0^{-1})^{-1}).$$
    %
    deltaB   = Btdraw(:,2:end) - Btdraw(:,1:end-1);
    sse_2W   = deltaB*deltaB';
    icQ      = chol(sse_2W + iW_prscale)\Is;
    Q        = icQ*icQ';
    Wdraw_CK = wishrnd(Q,T+W_prdf);
    Wdraw_CK2 = array2table(Wdraw_CK);
    writetable(Wdraw_CK2);
    iWdraw   = Wdraw_CK\Is;
    
    %% Sample $n_H$ and H jointly from a colapsed Gibbs.
    %
    %  $$H_{t} =\frac{1}{\lambda}\mathcal{U}(H_{t-1})^{\prime}\Psi_t\mathcal{U}(H_{t-1}),\qquad \Psi_t\sim B_m(n/2,0.5).$$
    % Note that all information from W is included on Btdraw, so we don't
    % need it here
    %
    
    ut = zeros(M,T);
    for i = 1:T
        ut(:,i)   = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i+1);
        %ut(:,i)   = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw; % for non TVP
        %coefficients
    end
    
    % Compute the marginal likelihood as in Proposition 3 of Windle &
    % Carvalho(2014) for a grid of lambda values.
    sum_ln_pYt = MarginalWish_n(ut,T,M,cSig0_draw,nugrid,1);
    pys        = exp(sum_ln_pYt-max(sum_ln_pYt))./sum(exp(sum_ln_pYt-max(sum_ln_pYt))); % normalized probabilities
    r          = mnrnd(1,pys); % draw from the multinomial
    ind        = r==1; % get the index of the grid point drawn.
    nu_Hdraw   = nugrid(ind); % defines the n draw.
    lam_Hdraw  = nu_Hdraw/(nu_Hdraw+1); % restriction in Appendix 3 of Uhlig (1997).

    % Sample H conditional on n using Widle & Carvalho (2014).
    Hdraw = BackSamplWish_chol(ut,T,M,cSig0_draw,lam_Hdraw,nu_Hdraw,1);
    
     %Hdraw2 = array2table(Hdraw); % I'm using this to debug my R code
     %writetable(Hdraw2); % I'm using this to debug my R code
    
    %% Sample H0 from Wishart posterior based on conjugate prior.  
    ciS_bar    = chol(lam_Hdraw*Hdraw(:,:,1)+invS0);
    cinvS      = ciS_bar\Im;
    S_bar      = cinvS*cinvS';
    cSig0_draw = chol(wishrnd(S_bar,nu_Hdraw+nu_H0));
    
    % This part (for fixed B) is not working yet
%     %% Sample Bdraw using Koop's algorithm
%     % I inverted because Carter and Kohn used an empty matrix Hdraw, but I
%     % need it with some values
%     % I got from Clark and Ravazzolo (2015) the following scheme
%     % B_0_prmean = Bols;
%     % B_0_prvar  = 4*Bols_var; % inflate OLS covariance matrix
%     soma_media = zeros(num_Betas,1);
%     soma_variancia = zeros(num_Betas, num_Betas);
%     for i = 1:T    
%         cH = chol(Hdraw(:,:,i));
%         icH = inv(cH);
%         R = icH*icH';
%         % Calculates the sum inside the posterior mean and variance
%         %soma_media = soma_media + Z(M*(i-1)+1:M*i,:)*y(:,i)*R;
%         % Bols    = (Zprior'*Zprior)\(Zprior'*yprior(:)); 
%         soma_media = soma_media + Z(M*(i-1)+1:M*i,:)'*R*y(:,i);
%         soma_variancia = soma_variancia + Z(M*(i-1)+1:M*i,:)'*R*Z(M*(i-1)+1:M*i,:);
%     end
%     cprvar = chol(B_0_prvar);
%     icprvar = inv(cprvar);
%     inprvar = icprvar*icprvar';
%     inv_variancia = inprvar + soma_variancia; 
%     cvariancia =  chol(inv_variancia);
%     inv_cvariancia = inv(cvariancia);
%     variancia = inv_cvariancia * inv_cvariancia';
%     media =  variancia*(soma_media + inprvar * B_0_prmean); 
%     % Sample the coefficients from the normal distribution
%     Btdraw = media + chol(variancia)'*randn(num_Betas,1); % Draw of alpha
    
%% Compute moments after burnin
    if irep > nburn
        nuH_save(irep-nburn) = nu_Hdraw;
        Bt_save(:,:,irep-nburn) = Btdraw;
        for t=1:T
            iHdraw(:,:,t) = Hdraw(:,:,t)\Im;
            temp = Hdraw(:,:,t);
            Ht_save(:,t,irep-nburn) = temp(tril(true(size(temp))));
        end
         if irep==nburn+1
             [deltaB2_M,deltaB2_S] = RunningStat(irep-nburn,deltaB.^2,[],[]);
             [ut_M,ut_S] = RunningStat(irep-nburn,ut,[],[]);
             [ut2_M,ut2_S] = RunningStat(irep-nburn,ut.^2,[],[]);
             [Bt_M,Bt_S] = RunningStat(irep-nburn,Btdraw,[],[]);
             [iW_M,iW_S] = RunningStat(irep-nburn,iWdraw,[],[]);
             [iH_M,iH_S] = RunningStat(irep-nburn,iHdraw,[],[]);
             [Sig0_M,Sig0_S] = RunningStat(irep-nburn,cSig0_draw'*cSig0_draw,[],[]);            
         else
             [Bt_M,Bt_S] = RunningStat(irep-nburn,Btdraw,Bt_M,Bt_S);
             [iW_M,iW_S] = RunningStat(irep-nburn,iWdraw,iW_M,iW_S);
             [iH_M,iH_S] = RunningStat(irep-nburn,iHdraw,iH_M,iH_S);
             [ut_M,ut_S] = RunningStat(irep-nburn,ut,ut_M,ut_S);
             [ut2_M,ut2_S] = RunningStat(irep-nburn,ut.^2,ut2_M,ut2_S);
             [deltaB2_M,deltaB2_S] = RunningStat(irep-nburn,deltaB.^2,deltaB2_M,deltaB2_S);
             [Sig0_M,Sig0_S] = RunningStat(irep-nburn,cSig0_draw'*cSig0_draw,Sig0_M,Sig0_S);            
         end
        
        %% Impulse response (got this from Koop)
    
         % Impulse response analysis. Note that Htsd contains the
         % structural error cov matrix
         % Set up things in VAR(1) format as in Lutkepohl (2005) page 51
            
             biga = zeros(M*p,M*p);
             for j = 1:p-1
                 biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
             end
             Htsd = Hdraw; % I don't want to change all names... This is the Omega matrix from Primiceri
             for i = 1:T % Get impulses recursively for each time period
                 bbtemp = Btdraw(M+1:s,i);  % get the draw of B(t) (exclude intercept) and repeat at all times #lazy_girl 
                 splace = 0;
                 for ii = 1:p
                     for iii = 1:M
                         biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                         splace = splace + M;
                     end
                 end
% 
%                 % ------------Identification code:                
%                 % St dev matrix for structural VAR
                 Hsd = chol(Htsd(:,:,i))';   % First shock is the Cholesky of the VAR covariance
                 diagonal = diag(diag(Hsd));
                 Hsd = inv(diagonal)*Hsd;    % Unit initial shock
%                 
%                 % Now get impulse responses for 1 through nhor future periods
                 impresp = zeros(M,M*nhor);
                 impresp(1:M,1:M) = Hsd; % First shock is the Cholesky of the VAR covariance
                 bigai = biga;
                 for j = 1:nhor-1
                     impresp(:,j*M+1:(j+1)*M) = bigj*bigai*bigj'*Hsd;
                     bigai = bigai*biga;
                 end
%                 
%                 % Only for specified periods - beginning, middle and end of
%                 % sample
                  if yearlab(i,2) == 20001;   % Jan 2000
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
%                     imp00(irep-nburn,:,:) = impf_m; % store draws of responses
                 end
                 if yearlab(i,2) == 20087;   % Jul 2008
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp08(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
                 if yearlab(i,2) == 201712;   % Dez 2017
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp17(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
% 
%                 % Only for specified periods - monetary policy             
                 if yearlab(i,2) == 20031;   % Jan 2003
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp03(irep-nburn,:,:) = impf_m; % store draws of responses
                 end
                 if yearlab(i,2) == 20063;   % March 2006
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp06(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
                 if yearlab(i,2) == 201111;   % Nov 2011
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp11(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
                 if yearlab(i,2) == 20151;   % Jan 2015
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp15(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
                 if yearlab(i,2) == 201512;   % Dez 2015
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp152(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
                 if yearlab(i,2) == 201611;   % Nov 2016
                     impf_m = zeros(M,nhor);
                     jj=0;
                     for ij = 1:nhor
                         jj = jj + M;    % restrict to the M-th equation, the interest rate
                         impf_m(:,ij) = impresp(:,jj);
                     end
                     imp16(irep-nburn,:,:) = impf_m;  % store draws of responses
                 end
            end %END geting impulses for each time period
         
     end
end

% Some plots from the output

figure(1)
plot(Bt_M')
ylim([-1 2])
title('Posterior mean of Bs')
snapnow

% Plot Hs
figure
subplot(3,1,1) 
plot(ut2_M(1,:))
hold on
plot(squeeze(iH_M(1,1,:)))
ylim([0 max(iH_M(1,1,:))*1.4])
xlim([1 T])
title('Posterior means for inv(H) and for u_t^2')
subplot(3,1,2) 
plot(ut2_M(2,:))
hold on
plot(squeeze(iH_M(2,2,:)))
ylim([0 max(iH_M(2,2,:))*1.4])
xlim([1 T])
subplot(3,1,3) 
plot(ut2_M(3,:))
hold on
plot(squeeze(iH_M(3,3,:)))
ylim([0 max(iH_M(3,3,:))*1.4])
xlim([1 T])
snapnow


figure
subplot(1,2,1)
ksdensity(nuH_save./(nuH_save+1))
title('Posterior density of \lambda_H')
subplot(1,2,2)
ksdensity(nuH_save)
title('Posterior density of n_H')
snapnow

% figure
% for i=1:s
%     subplot(7,3,i)
%     plot(deltaB2_M(i,:))
%     xlim([1 T])
% end
% suptitle('Posterior means of \eta^2')
% 
% % Compute mean and standard deviation of diagonal elements of posterior and prior distributions for W.
% Ordering = {'B0_1';'B0_2';'B0_3';'B11_1';'B12_1';'B13_1';'B21_1';'B22_1';'B23_1';'B31_1';'B32_1';'B33_1';'B11_2';'B12_2';'B13_2';'B21_2';'B22_2';'B23_2';'B31_2';'B32_2';'B33_2';'n_H'};
% iW_poststd=iW_S/(nrep-1); % posterior standard deviation.
% iW_primean = iW_prscale./(W_prdf-s-1); % prior mean
% diagiW_pristd = sqrt(2*diag(iW_prscale.^2)./((W_prdf-s-1)^2*(W_prdf-s-3))); % prior standard deviation
% % Write results in a table
% tabela_M4 = table([diag(iW_M);mean(nuH_save)],[diag(iW_poststd);std(nuH_save)],[diag(iW_primean);0.5*(max(nugrid)-min(nugrid))],[diagiW_pristd;sqrt(1/12*(max(nugrid)-min(nugrid))^2)],'RowNames',Ordering,'VariableNames',{'PostMean_W' 'PostStD_W' 'PriorMean_W' 'PriorStD_W'}) % table of means and standard deviations
% 
% disp(Sig0_M)
% 
qus = [.16, .5, .84];
qus2 = [.5, .25, .5, .75, .95];
    imp00XY=squeeze(quantile(imp00,qus));
    imp08XY=squeeze(quantile(imp08,qus));
    imp17XY=squeeze(quantile(imp17,qus));
    
    imp00XY2=squeeze(quantile(imp00,qus2));
    imp08XY2=squeeze(quantile(imp08,qus2));
    imp17XY2=squeeze(quantile(imp17,qus2));
    
    imp03XY=squeeze(quantile(imp03,qus));
    imp06XY=squeeze(quantile(imp06,qus));
    imp11XY=squeeze(quantile(imp11,qus));
    imp15XY=squeeze(quantile(imp15,qus));
    imp152XY=squeeze(quantile(imp152,qus));
    imp16XY=squeeze(quantile(imp16,qus));
    
    imp03XY2=squeeze(quantile(imp03,qus2));
    imp06XY2=squeeze(quantile(imp06,qus2));
    imp11XY2=squeeze(quantile(imp11,qus2));
    imp15XY2=squeeze(quantile(imp15,qus2));
    imp152XY2=squeeze(quantile(imp152,qus2));
    imp16XY2=squeeze(quantile(imp16,qus2));
    
 
    % Plot impulse responses
    figure       
    %set(0,'DefaultAxesColorOrder',[0 0 0],...
        %'DefaultAxesLineStyleOrder','--|-|--')
    subplot(3,5,1)
    plot(1:nhor,squeeze(imp03XY(:,1,:)))
    title('Impulse response of capital labor ratio, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,2)
    plot(1:nhor,squeeze(imp03XY(:,2,:)))
    title('Impulse response of inflation rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,3)
    plot(1:nhor,squeeze(imp03XY(:,3,:)))
    title('Impulse response of GDP, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,4)
    plot(1:nhor,squeeze(imp03XY(:,4,:)))
    title('Impulse response of interest rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,5)
    plot(1:nhor,squeeze(imp03XY(:,5,:)))
    title('Impulse response of exchange rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    %%%%%%
    subplot(3,5,6)
    plot(1:nhor,squeeze(imp06XY(:,1,:)))
    title('Impulse response of capital labor ratio, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,7)
    plot(1:nhor,squeeze(imp06XY(:,2,:)))
    title('Impulse response of inflation rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,8)
    plot(1:nhor,squeeze(imp06XY(:,3,:)))
    title('Impulse response of GDP, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,9)
    plot(1:nhor,squeeze(imp06XY(:,4,:)))
    title('Impulse response of interest rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,10)
    plot(1:nhor,squeeze(imp06XY(:,5,:)))
    title('Impulse response of exchange rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    %%%%%%
    subplot(3,5,11)
    plot(1:nhor,squeeze(imp11XY(:,1,:)))
    title('Impulse response of capital labor ratio, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,12)
    plot(1:nhor,squeeze(imp11XY(:,2,:)))
    title('Impulse response of inflation rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,13)
    plot(1:nhor,squeeze(imp11XY(:,3,:)))
    title('Impulse response of GDP, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,14)
    plot(1:nhor,squeeze(imp11XY(:,4,:)))
    title('Impulse response of interest rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,15)
    plot(1:nhor,squeeze(imp11XY(:,5,:)))
    title('Impulse response of exchange rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
snapnow


    % Plot impulse responses
    figure       
    %set(0,'DefaultAxesColorOrder',[0 0 0],...
        %'DefaultAxesLineStyleOrder','--|-|--')
    subplot(3,5,1)
    plot(1:nhor,squeeze(imp15XY(:,1,:)))
    title('Impulse response of capital labor ratio, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,2)
    plot(1:nhor,squeeze(imp15XY(:,2,:)))
    title('Impulse response of inflation rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,3)
    plot(1:nhor,squeeze(imp15XY(:,3,:)))
    title('Impulse response of GDP, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,4)
    plot(1:nhor,squeeze(imp15XY(:,4,:)))
    title('Impulse response of interest rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,5)
    plot(1:nhor,squeeze(imp15XY(:,5,:)))
    title('Impulse response of exchange rate, Jan/2000')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    %%%%%%
    subplot(3,5,6)
    plot(1:nhor,squeeze(imp152XY(:,1,:)))
    title('Impulse response of capital labor ratio, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,7)
    plot(1:nhor,squeeze(imp152XY(:,2,:)))
    title('Impulse response of inflation rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,8)
    plot(1:nhor,squeeze(imp152XY(:,3,:)))
    title('Impulse response of GDP, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,9)
    plot(1:nhor,squeeze(imp152XY(:,4,:)))
    title('Impulse response of interest rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,10)
    plot(1:nhor,squeeze(imp152XY(:,5,:)))
    title('Impulse response of exchange rate, Jan/2008')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    %%%%%%
    subplot(3,5,11)
    plot(1:nhor,squeeze(imp16XY(:,1,:)))
    title('Impulse response of capital labor ratio, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(3,5,12)
    plot(1:nhor,squeeze(imp16XY(:,2,:)))
    title('Impulse response of inflation rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(3,5,13)
    plot(1:nhor,squeeze(imp16XY(:,3,:)))
    title('Impulse response of GDP, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(3,5,14)
    plot(1:nhor,squeeze(imp16XY(:,4,:)))
    title('Impulse response of interest rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
    subplot(3,5,15)
    plot(1:nhor,squeeze(imp16XY(:,5,:)))
    title('Impulse response of exchange rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
snapnow


figure
    subplot(2,2,1)
    plot(1:nhor,squeeze(imp15XY(:,1,:)))
    %title('Impulse response of capital labor ratio, Dec/2017')
    title('Impulse response of K/L, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)
    subplot(2,2,2)
    plot(1:nhor,squeeze(imp15XY(:,2,:)))
    title('Impulse response of inflation rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)    
    subplot(2,2,3)
    plot(1:nhor,squeeze(imp15XY(:,3,:)))
    title('Impulse response of GDP, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor)  
    subplot(2,2,4)
    plot(1:nhor,squeeze(imp15XY(:,4,:)))
    title('Impulse response of exchange rate, Dec/2017')
    xlim([1 nhor])
    set(gca,'XTick',0:3:nhor) 
snapnow   


figure
     set(0,'DefaultAxesColorOrder',[0 0 0],...
        'DefaultAxesLineStyleOrder','-|:|--',...
        'DefaultLineLineWidth',1)
     plot(1:nhor,squeeze(imp00XY(2,1,:)))
   %  title('Impulso resposta da razão capital trabalho')
    hold on
    plot(1:nhor,squeeze(imp08XY(2,1,:)))
    hold on
    plot(1:nhor,squeeze(imp17XY(2,1,:)))
    xlim([1 nhor])
    %ylim([0 0.04])
  legend('Jan 2000','Jul 2008','Dec 2017')
    
   figure
     set(0,'DefaultAxesColorOrder',[0 0 0],...
        'DefaultAxesLineStyleOrder','-|:|--',...
        'DefaultLineLineWidth',1)
     plot(1:nhor,squeeze(imp00XY(2,3,:)))
   %  title('Impulso resposta do IPCA')
    hold on
    plot(1:nhor,squeeze(imp08XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(imp17XY(2,3,:)))
    xlim([1 nhor])
    legend('Jan 2000','Jul 2008','Dec 2017')
