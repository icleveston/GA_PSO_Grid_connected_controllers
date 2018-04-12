function resultados = projeta_K_deadbeat_GA
A_til_1 = 0;
B_til = 0;
Br_til = 0;
Bdist1_til = 0;
C_til = 0;
Ts = 0;

load('deadbeat_ga.mat', 'Acli', 'A1', 'B1', 'A2', 'B2','A3', 'B3', ...
	'A4', 'B4','Ts')
% ---------------------------------------------------------------------------------------------//-------------------------------------------------------

%ganhos H2 - produz menor sobrecorrente na partida
% Kfull2=[-9.747423809401397  -1.851225100912732  -0.350767371688946  -0.446663368860505  50.868719652179720 -50.184356792361314  16.232137022975223 -15.756612317015055   8.734346658571226  -8.843404878241969 3.712318389421682  -4.140705855357965];

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Simula��o
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% to = 0;
% tf = 0.03;
% 
% t = to:Ts:tf;       
% [n1 n2] = size(t);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Dist�rbio
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% disturbio = 1*(311)*sin(2*pi*60*t);%+3*sin(2*pi*180*t); 
% 
% 
% for nsample = 1:n2
%           
%           if (nsample < 2*334)
%       	       ref(nsample) = 0;
%           end
%           if (nsample >= 2*334 && nsample < 4.75*334)
%       	       ref(nsample) = 10.*sin(60*2*pi*Ts*nsample);
%           end
%           if (nsample >= 4.75*334 && nsample < 9*334)
%       	       ref(nsample) = -10.*sin(60*2*pi*Ts*nsample);
%           end
%           if (nsample >= 9*334 && nsample < 13*334)
%       	       ref(nsample) = 10.*cos(60*2*pi*Ts*nsample);
%           end
%           if (nsample >= 13*334)
%       	       ref(nsample) = 20.*cos(60*2*pi*Ts*nsample);
%           end
% 
% end
% 
% u = [ref; disturbio;];
% 
% % MF1full2 = ss(A_til_1 + B_til*Kfull2, [Br_til Bdist1_til], C_til, 0, Ts ); 
% % MF2full2 = ss(A_til_2 + B_til*Kfull2, [Br_til Bdist2_til], C_til, 0, Ts ); 



%%
% N�mero de vari�veis
nvar = 4;
% Limites

% limsup = [100 100 100 100];
% liminf = [-100 -100 -100 -100];

% limsup = [-20.6270 -0.15614 250.1633 -150.7911];
% liminf = [-60.6270 -2.05614 150.1633 -250.7911];

% K(r=1)  -8.935078398724579  -0.211068302121902  22.187574090816270 -21.543163811119925

% limsup = [-6 0.1 25 -15];
% liminf = [-10 -2 15 -25];

% K(0.95)  -22.6270   -0.5614   79.1633  -76.7911 

limsup = [-20.6270 0.15614 85.1633 -65.7911];
liminf = [-25.6270 -2 65.1633 -85.7911];

% K(0.92) [-43.0775 -1.1616 201.6442 -192.7824]

% limsup = [-42 -1 202 -191];
% liminf = [-44 -2 200 -193];

% K para Lnom: K=acker(A1,B1,[exp(-500*Ts) exp(-1000*Ts) exp(-1500*Ts) exp(-2000*Ts)])
% K = [-3.651455891318733  0.538926310262791  6.974440469643639   -6.955045668469140

% limsup = [-1 1 10 -4];
% liminf = [-6 -1 4 -8];

% inicializando limites com ganhos maximos e minimos de deadbeat dos 4
% vertices nao funciona
% limsup =   1.0e+03 * [    0.4794    0.0030   -2.5509    7.6655];
% liminf=   1.0e+04 *  [    0.0119    0.0003   -1.0204    0.1916]

% DB LMIs

%realimentacao parcial de estados
% limsup = [50 0 500.500];
% liminf = [-50 0 -500 -500];

% falha no sensor de corrente => nao consegue estabilizar
% limsup = [0 50 500.500];
% liminf = [0 -50 -500 -500];
% limsup = [0 2 2 2];
% liminf = [0 -2 -2 -2];

% 
% Configura��es do ga
%setupga = gaoptimset(...
%	'EliteCount', 0,...
%	'CrossoverFraction', 0.3,...
%	'Generations', 10,...
%	'StallGenLimit', 20,...
%	'Display', 'iter',...
%	'PopulationSize', 100,...
%	'SelectionFcn', @selectiontournament,...
%	'CrossoverFcn', @crossoverscattered, ...
%	'MutationFcn', @mutationadaptfeasible);
%
%% Executa o GA
%qr_dlqr = ga(@(x)fobj(x), nvar,...
%	[], [], ... % Desigualdades lineares
%	[], [],... % Igualdades
%	liminf, limsup,... % Limites de busca
%	[],... % Restri��es
%	setupga);
%
%resultados = qr_dlqr;

% % Fecha a malha e simula
% Kfull = -dlqr(A_til_1, B_til, diag(qr_dlqr(1:12)), qr_dlqr(13));
% MF1full = ss(A_til_1 + B_til*Kfull, [Br_til Bdist1_til], C_til, 0, Ts); 
% yh21full2 =  lsim(MF1full, u, t, zeros(1,12));
% plot(t, yh21full2)

Kdb = [-22.6270   -0.5614   79.1633  -76.7911 ];
Acli=[A1+B1*Kdb A2+B2*Kdb A3+B3*Kdb A4+B4*Kdb];
raio = nuvem_d(Acli,0)
res = calcula_itae(Kdb)

result = res^raio

%%
function res = fobj(genes)
	% Projeta o controlador
	Kdb = genes(1:4);
    Acli=[A1+B1*Kdb A2+B2*Kdb A3+B3*Kdb A4+B4*Kdb];
	% Calcula a fun��o objetivo
 	
    %res = calcula_itae(Kdb);
    
    raio = nuvem_d(Acli,0);
    res = calcula_itae(Kdb)^raio;
    
%     raio = nuvem_d(Acli,0);
%     if raio<1
%        res = calcula_itae(Kdb) + 1e7*raio;
%     else
%         res = calcula_itae(Kdb) + 1e10*raio;
%     end
        
%     if res >1
%         res=2;
%     end
end % Fim da fun��o objetivo
end