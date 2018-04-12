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

% ganhos para enviar para o PSIM
%[-2.3413    0.7883    1.1917   -1.2580]

%%
% Número de variáveis
nvar = 4;
% Limites
%limsup = [-20.6270 -0.15614 250.1633 -150.7911];
%liminf = [-60.6270 -2.05614 150.1633 -250.7911];

limsup = [100 100 100 100];
liminf = [-100 -100 -100 -100];

% Configurações do ga
setupga = gaoptimset(...
	'EliteCount', 0,...
	'CrossoverFraction', 0.3,...
	'Generations', 30,...
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




%%
function res = fobj(genes)
	% Projeta o controlador
	Kdb = genes(1:4);
    Acli=[A1+B1*Kdb A2+B2*Kdb A3+B3*Kdb A4+B4*Kdb];
	% Calcula a função objetivo
	res = nuvem_d(Acli,0);
%     if res >1
%         res=2;
%     end
end % Fim da função objetivo
end