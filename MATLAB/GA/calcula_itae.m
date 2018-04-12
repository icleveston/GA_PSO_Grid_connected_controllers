function [itae] = calcula_itae(K)
% -----------------------------------------------------------------------
% SIMULA��O DLQR com PWM
% Inversor monof�sico conectado a rede - filtro L
% -----------------------------------------------------------------------
% Controle de Vers�es
% Caio Os�rio - V5
%
% -----------------------------------------------------------------------
% Altera��es e Observa��es:
%
% (1) � poss�vel variar a indutancia durante a simula��o alterando os valores 
% de Lg e Lg2 e o temo em que s�o alterados
%
% (2) Agora fs w fsw s�o diferentes, � poss�vel alter�-las. 
%
% (3) O momento em que � aplicada a varia��o de indut�ncia influi no
% comportamento
% -----------------------------------------------------------------------
% clear all
% close all
% clc

%% -----------------------------------------------------------------------
% $$$$$$$$$$$$$$$$$   Projeto Controlador Robusto   $$$$$$$$$$$$$$$$$$$$$
% ------------------------------------------------------------------------

wn=2*pi*60;
csi=0.001;
R=[0 1;-wn^2 -2*csi*wn];
S=[0;1];
Ts=1/10000;  % Frequencia de amostragem !!!

% Implementacao mais favoravel do ressonante
s     = tf('s'); 
z     = tf('z',Ts);
w = 2*pi*60;
zeta    = 1*0.0001;

G_res    = w^2/(s^2 + 2*zeta*s + w^2);% trocamos o s do numerador por w^2
G_res_discreto    = c2d(G_res,Ts,'tustin');

[Rd,Sd,U1,V1] = ssdata(G_res_discreto);


% Par�metros para projeto robusto ----------------------------------------

R1=0.1;
Lg_nom=5e-3;
delta_Lg=3e-3;
delta_Rg=0.1;

Lmin=Lg_nom-delta_Lg;
Lmax=Lg_nom+delta_Lg;
Rlmin=R1-delta_Rg;
Rlmax=R1+delta_Rg;

C=[1 0];
D=0;
a11=(-Rlmin/Lmin)*Ts +1;
b11=Ts/Lmin;

a12=(-Rlmax/Lmin)*Ts +1;
b12=Ts/Lmin;

a21=(-Rlmin/Lmax)*Ts +1;
b21=Ts/Lmax;

a22=(-Rlmax/Lmax)*Ts +1;
b22=Ts/Lmax;

A1=[a11 b11 zeros(1,2); 0 0 zeros(1,2); -Sd zeros(2,1) Rd];
B1=[0;1;zeros(2,1)];
Br1=[0;0;Sd];
Bd1=[-b11;0;zeros(2,1)];

A2=[a12 b12 zeros(1,2); 0 0 zeros(1,2); -Sd zeros(2,1) Rd];
B2=[0;1;zeros(2,1)];
Br2=[0;0;Sd];
Bd2=[-b12;0;zeros(2,1)];

A3=[a21 b21 zeros(1,2); 0 0 zeros(1,2); -Sd zeros(2,1) Rd];
B3=[0;1;zeros(2,1)];
Br3=[0;0;Sd];
Bd3=[-b21;0;zeros(2,1)];

A4=[a22 b22 zeros(1,2); 0 0 zeros(1,2); -Sd zeros(2,1) Rd];
B4=[0;1;zeros(2,1)];
Br4=[0;0;Sd];
Bd4=[-b22;0;zeros(2,1)];

Ai=[A1 A2 A3 A4];
Bi=[B1 B2 B3 B4];
Bri=[Br1 Br2 Br3 Br4];
Bdi=[Bd1 Bd2 Bd3 Bd4];

% r=1  % Raio para aloca��o de polos
% out = ssf_stab_K_d_mtb(Ai/r,Bi/r)
% K=out.K;


% $$$$$$$$$$$$$$$$$$$$$$$$ FIM DO PROJETO ROBUSTO  $$$$$$$$$$$$$$$$$$$$$$$

%% -----------------------------------------------------------------------
%                        PARAMETROS DE SIMULACAO
%  -----------------------------------------------------------------------

Cic= 2;                                     % N�mero de c�clos de rede simulados
fr = 60;                                    % Frequ�ncia de sa�da
Tr = 1/fr;                                  % Per�odo da tens�o da rede
fs = 10000;                                  % Frequencia de amostragem
Ts = 1/fs;
fsw = 10000;
Tsw = 1/fsw;

% Tempo
PerT = Cic*Tr;                              % Per�odo de simula��o
dT = 1e-6;                                  % Passo da simula��o
t = (0:dT:PerT)';                           % Vetor do tempo total da simula��o
tsim = (0:Ts:PerT)';

% Pontos
Pontos_fs = 1/(dT*fs);                      % Pontos da simula��o (por ciclo de Ts)
Pontos_fsw = 1/(dT*fsw);
Pontos_fr = floor(Pontos_fs*Tr/Ts);         % Pontos em um ciclo da rede
Pontos_t = length(t);                       % Total de pontos

% Par�metros do sistema
Vcc = 200;                                  % Tens�o do barramento CC (V)
Lg = 2e-3;                                  % Indut�ncia da rede inicial (valor real)
Lg2 = 2e-3;                                 % Indut�ncia da rede ap�s varia��o (valor real)
Lg_nom = 5e-3;                              % Indut�ncia da rede (valor nominal - utilizado para projeto)
Rf = 0.1;                                   % Resist�ncia do filtro de sa�da (valor real)
Rf2 = 0.1;
Rf_nom = 0.1;                               % Resist�ncia do filtro de sa�da (valor nominal - projeto)
vg_pk = 127*sqrt(2);                        % Tens�o da rede (disturbio)
Ma = vg_pk/Vcc;                             % �ndice de modula��o de amplitude
w = 2*pi*fr;                                % Frequ�ncia angular
ig_ref_pk = 10;                             % Corrente de referencia (peak)

% Inicializa��es
t_k=0;
t_ks=0;
vtr_tk=zeros(Pontos_t,1);

upwm=[];
x(1)=0;
theta(1)=0;
rho1(1)=0;
rho2(1)=0;
xc(1)=0;
u(1)=0;
upwm_k=0;
ref(1)=0;
cont=1;
u_ks=0;
ref_ks=0;
ig_amost=0;
ref_amost=0;

ks=1;

%% Controlador Ressonante

s     = tf('s');
z     = tf('z',Ts);
zeta  = 1*0.0001;

G_res    = w^2/(s^2 + 2*zeta*s + w^2);% trocamos o s do numerador por w^2
G_res_discreto    = c2d(G_res,Ts,'tustin');
[Rd,Sd,U1,V1] = ssdata(G_res_discreto);

        
%% Modelo do conversor em espa�o de estados

% Planta utilizada para projeto (valores nominais)
a11 = (-Rf/Lg_nom)*Ts +1;
b11 = Ts/Lg_nom;

A1=[a11 b11 zeros(1,2); 0 0 zeros(1,2); -Sd zeros(2,1) Rd];
B1=[0;1;zeros(2,1)];

%% Controlador dlqr

% K=dlqr(A1,B1,diag([10 1 10 10]),10);
% K=-K % implementa realimentacao positiva de estados

%% Projeto de um Deadbeat Convencional

% K=acker(A1,B1,0.75*[1 1 1 1]);
% %K = [299.2453 2.9966 -6.3773e+003 4.7909e+003]
% K=-K;



%% La�o da varia��o do tempo da simula��o

        for k=1:(Pontos_t),   % k � o tempo "continuo"  
            
            if t_k>0.035 %0.03474%0.020833
                % Planta real (continua)
                an=(-Rf2/Lg2)*(dT) +1;
                bn=dT/Lg2;
            else 
                % Planta real (continua)
                an=(-Rf/Lg)*(dT) +1;
                bn=dT/Lg;
            end
            
            ref_k(k)= ig_ref_pk*sin(w*t_k); % corrente de referencia em t_k
            
            % Amostragem, ks � o tempo discreto
            if (mod(k,floor(Pontos_fs))==0),
                
             
                ref(ks) = ig_ref_pk*sin(w*t_ks); % corrente de referencia
                u(ks)=K(1)*x(ks)+K(2)*theta(ks)+K(3)*rho1(ks)+K(4)*rho2(ks);
                u_ks=u(ks);
       
                x(ks+1) = xc(k); %amostragem da variavel de saida c atraso
                theta(ks+1)=u(ks); 
                rho1(ks+1)=-Sd(1,1)*xc(k)+0*theta(ks)+Rd(1,1)*rho1(ks)+Rd(1,2)*rho2(ks)+0*u(ks)+Sd(1,1)*ref(ks);
                rho2(ks+1)=-Sd(2,1)*xc(k)+0*theta(ks)+Rd(2,1)*rho1(ks)+Rd(2,2)*rho2(ks)+0*u(ks)+Sd(2,1)*ref(ks);
                
                ig_amost = xc(k);
                ref_amost = ref(ks);
                
                
                ks=ks+1;         
                t_ks = t_ks + Ts;
                
                
            end

            % Modula��o phase-shift ---------------------------------------
            % Nesta t�cnica, frequencia efetiva = 2*fsw.
            
            v_tri = 2*asin(sin(2*pi*t_k/Tsw-pi/2))/pi;
   
            if (u_ks/Vcc > v_tri) Sa=1;
            else Sa=0;
            end
            
            if (-u_ks/Vcc > v_tri) Sb=1;
            else Sb=0;
            end
            
            upwm_k = (Sa-Sb)*Vcc;
            upwm(k) = upwm_k;
            
            % ------------------------------------------------------------- 
            
            
            % Gera��o do PWM centralizado com a triangular
%             if cont<=round(((Vcc-abs(u_ks))/(2*Vcc))*floor(Pontos_fsw))
%                upwm(cont)=0;
%             end
%             if cont>round(((Vcc-abs(u_ks))/(2*Vcc))*floor(Pontos_fsw)) & cont<floor(Pontos_fsw)-(round(((Vcc-abs(u_ks))/(2*Vcc))*floor(Pontos_fsw)))
%                upwm(cont)=sign(u_ks)*Vcc;
%             end
%             if cont>=floor(Pontos_fsw)-(round(((Vcc-abs(u_ks))/(2*Vcc))*floor(Pontos_fsw))) & cont<=floor(Pontos_fsw)
%                upwm(cont)=0;
%             end
%         
%             upwm_k=upwm(cont);
            
        
            % Disturbio - tensao da rede
            vg(k) = vg_pk*sin(w*t_k);
            vtr_tk(k) = t_k;
            
            % Modelo do Conversor (real)
            xc(k+1) = an*xc(k) + bn*upwm_k - bn*vg(k);
            
             
            % Atualiza��es
            
            vtr_ref_amost(k)=ref_amost;
            vtr_u_amost(k)=u_ks;
            vtr_ig_amost(k)=ig_amost;
            
            t_k = t_k + dT;
            k = k+1;
            
            if cont<floor(Pontos_fsw)
                cont=cont+1;
            else
                    cont=1;
            end
            
        end
        
        
        
%% M�tricas para qualidade da resposta

itae=0;
for cont=1:8334*2
    erro(cont) = ref_k(cont)-xc(cont);
    itae = itae + cont*abs(erro(cont));
end






% Gr�ficos
%         
plot(vtr_tk,xc(1:Pontos_t))
hold on
plot(vtr_tk,ref_k(1:Pontos_t),'r')
plot(vtr_tk,vtr_ref_amost(1:Pontos_t),'r')
grid on

figure
plot(vtr_tk,vtr_ig_amost(1:Pontos_t))
hold on
plot(vtr_tk,vtr_ref_amost(1:Pontos_t),'r')
legend('medida','referencia')

figure
plot(vtr_tk,vtr_u_amost(1:Pontos_t),'*')
