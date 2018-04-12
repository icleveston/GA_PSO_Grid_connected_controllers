function resultados = projeta_K_deadbeat_GA_FF_v1
A_til_1 = 0;
B_til = 0;
Br_til = 0;
Bdist1_til = 0;
C_til = 0;
Ts = 0;

%load('matrizes.mat', 'A1', 'B1', 'A2', 'B2','L1', 'L2')
load('matrizes_orig.mat', 'A1', 'B1', 'A2', 'B2','L1', 'L2')


% ---------------------------------------------------------------------------------------------//-------------------------------------------------------

%ganhos H2 - produz menor sobrecorrente na partida
% Kfull2=[-9.747423809401397  -1.851225100912732  -0.350767371688946  -0.446663368860505  50.868719652179720 -50.184356792361314  16.232137022975223 -15.756612317015055   8.734346658571226  -8.843404878241969 3.712318389421682  -4.140705855357965];

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Simulação
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% to = 0;
% tf = 0.03;
% 
% t = to:Ts:tf;       
% [n1 n2] = size(t);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Distúrbio
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
% Número de variáveis
nvar = 5;
% Limites
%limsup = [-20.6270 -0.15614 250.1633 -150.7911];
%liminf = [-60.6270 -2.05614 150.1633 -250.7911];

% Inicialização para o caso matrizes.mat
% limsup = [-30 1 200 -120 10];
% liminf = [-50 -2 100 -200 -10];
% [-37.3717 -0.4763 169.9207 -162.5174] gnahos para r=0.92

% Inicialização para o caso matrizes_orig.mat
limsup = [-20 2 300 -120 2];
liminf = [-50 -3 100 -300 -2];
% [-42.4746 -1.1353 198.9936 -190.2412] 
% A otimização encontrada foi:
%   100         10100          0.8459          0.8508        1 
% Ganhos =
%   -24.3462   -1.3701  150.4616 -140.3833    2.1072

% Configurações do ga
setupga = gaoptimset(...
	'EliteCount', 0,...
	'CrossoverFraction', 0.3,...
	'Generations', 100,...
	'StallGenLimit', 20,...
	'Display', 'iter',...
	'PopulationSize', 100,...
	'SelectionFcn', @selectiontournament,...
	'CrossoverFcn', @crossoverscattered, ...
	'MutationFcn', @mutationadaptfeasible);

% Executa o GA
qr_dlqr = ga(@(x)fobj(x), nvar,...
	[], [], ... % Desigualdades lineares
	[], [],... % Igualdades
	liminf, limsup,... % Limites de busca
	[],... % Restrições
	setupga);

resultados = qr_dlqr;

% % Fecha a malha e simula
% Kfull = -dlqr(A_til_1, B_til, diag(qr_dlqr(1:12)), qr_dlqr(13));
% MF1full = ss(A_til_1 + B_til*Kfull, [Br_til Bdist1_til], C_til, 0, Ts); 
% yh21full2 =  lsim(MF1full, u, t, zeros(1,12));
% plot(t, yh21full2)



%%
function res = fobj(genes)
	% Projeta o controlador
	Kdb = genes(1:5);
    Acli=[A1+B1*Kdb(1,1:4)+B1*L1*Kdb(1,5) A2+B2*Kdb(1,1:4)+B2*L2*Kdb(1,5)];
	% Calcula a função objetivo
	res = nuvem_d(Acli,0);
%     if res >1
%         res=2;
%     end
end % Fim da função objetivo
end