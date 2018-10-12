%% Implementa Carter and Kohn (1994), On Gibbs sampling for state space models.
function [bdraw] = carter_kohn_MSV_1(y,Z,iRt,Q,m,p,T,B0,V0)%,semente)
% Função que implementa o procedimento proposto por Carter and Kohn (1994)
% para simular estados de uma normal T-variada para ser usado dentro de um
% algoritmo de Gibbs.
%
% $$y_t=Z_tb_t+e_t,  e_t\sim N(0,R_t),$$
%
% $$b_{t+1}=b_t+u_t, u_t\sim N(0,Q).$$
%
%

% Fixa semente:
% s = RandStream('mt19937ar','Seed',semente); % o primeiro argumento fixa o algoritmo gerador de números pseudo-aleatórios, o segundo indica que iremos escolher o valor da semente e o terceiro é o próprio valor da semente.
% RandStream.setGlobalStream(s); % usa a semente s como semente para o gerador de números aleatórios de todas as distribuições.

%%%%%%%%%%%%%% 
% DEBUG
% carter_kohn_MSV_1(y,Z,Hdraw,iWdraw,s,M,T,B_0_prmean,B_0_prvar)

% iRt = Hdraw;
% Q = iWdraw;
% m = s;
% p = M;
% B0 = B_0_prmean;
% V0 = B_0_prvar;

% END OF DEBUG
%%%%%%%%%%%%%%


% Valores iniciais
btt = B0; % média da a priori de beta em t=1.
Vtt = V0; % variância a priori de beta em t=1.
bt  = zeros(T,m); % guarda espaço para salvar as médias da distribuição filtrada.
Vt  = zeros(m^2,T); % guarda as variâncias da distribuição filtrada.
% log_lik = 0; % inicializa o log da verossimilhança
% Q = inv(iQ);

%% O Filtro de Kalman 
% Essa parte do código utiliza o filtro de Kalman para calcular e salvar
% algumas quantidades que serão necessárias mais tarde no suavizador de
% Kalman.
% 
for i=1:T
    bp = btt;
    Vp = Vtt+Q;
    cH = chol(iRt(:,:,i));
    icH = inv(cH);
    R = icH*icH';
    H = Z((i-1)*p+1:i*p,:);
    cfe = y(:,i) - H*bp;   % erro de previsão um passo a frente
    VpH = Vp*H';
    f   = H*VpH + R;     % variância do erro de previsão um passo a frente
    Mt  = VpH/f;
    btt = bp + Mt*cfe;
    Vtt = Vp-Mt*f*Mt';
%     if ~issymmetric(Vtt)
%         Vtt = (Vtt+Vtt')/2;
%     end
    bt(i,:) = btt';
    Vt(:,i) = Vtt(:);
end

% Amostre beta de t=T da densidade filtrada beta(T|T) ~ N(bt(:,T),Vt(:,T))
bdraw = zeros(T,m); % guarda espaço para os betasque serão amostrados.
if ~issymmetric(Vtt)       
    Vtt = (Vtt+Vtt')/2;
end

bdraw(T,:) = mvnrnd(btt,Vtt); % amostra beta de t=T

%% Amostragem de trás para a frente
% Usa as recursões do suavizador para amostrar das a posterioris de todos
% os betas para todos os t.

for i=1:T-1
    bf  = bdraw(T-i+1,:)'; % beta amostrado em t+1.
    btt = bt(T-i,:)'; % média da distribuição filtrada em t.
    Vtt = reshape(Vt(:,T-i),m,m); % variância da distribuição filtrada em t.
    f   = Vtt + Q;
    cfe = bf - btt;
    VF = Vtt/f;
    bmean = btt + VF*cfe; % média a posteriori
    bvar = Vtt - VF*Vtt;  % variância a posterirori (VF*VF)
    if ~issymmetric(bvar)       
        bvar = (bvar+bvar')/2;
    end
    bdraw(T-i,:) = mvnrnd(bmean,bvar); % amostra beta de sua a posteriori.
end

%bdraw2 = array2table(bdraw);
%writetable(bdraw2);


%% Amostra b0|T
bf  = bdraw(1,:)'; % beta amostrado em t=1.
btt = B0;          % média da distribuição inicial para b0.
Vtt = V0;          % variância da distribuição inicial para V0.
f   = Vtt + Q;
cfe = bf - B0;
VF = Vtt/f;
bmean = btt + VF*cfe; % média a posteriori
bvar = Vtt - VF*Vtt;  % variância a posterirori (VF*VF)   
if ~issymmetric(bvar)       
    bvar = (bvar+bvar')/2;
end
b0draw = mvnrnd(bmean,bvar); % amostra beta de sua a posteriori.

bdraw = [b0draw;bdraw]'; % includes initial condition on the vector of betas.