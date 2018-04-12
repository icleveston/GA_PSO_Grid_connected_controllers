function resultados = projeta_K_deadbeat_GA_FF_v1

%load('matrizes.mat', 'A1', 'B1', 'A2', 'B2','L1', 'L2')
load('matrizes_orig.mat', 'A1', 'B1', 'A2', 'B2','L1', 'L2')

% N?mero de vari?veis
nvar = 5;

% Inicializa??o para o caso matrizes_orig.mat
limsup = [-20 2 300 -120 2];
liminf = [-50 -3 100 -300 -2];

% Configura??es do ga
setupga = gaoptimset(...
	'EliteCount', 0,...
	'CrossoverFraction', 0.3,...
	'Generations', 10,... % 100
	'StallGenLimit', 20,...
	'Display', 'iter',...
	'PopulationSize', 10,... % 100
	'SelectionFcn', @selectiontournament,...
	'CrossoverFcn', @crossoverscattered, ...
	'MutationFcn', @mutationadaptfeasible);

% Executa o GA
qr_dlqr = ga(@(x)fitness(x), nvar,...
	[], [], ... % Desigualdades lineares
	[], [],... % Igualdades
	liminf, limsup,... % Limites de busca
	[],... % Restri??es
	setupga);

function res = fitness(genes)
    
    i = sym('i');
    i1 = sym('i1');
    i2 = sym('i2');
    vout = sym('vout');
       
    % Adiciona os genes para o workspace
    assignin('base', 'k1', genes(1));
    assignin('base', 'k2', genes(2));
    assignin('base', 'k3', genes(3));
    assignin('base', 'k4', genes(4));
    assignin('base', 'k5', genes(5));
    
    assi
    
	% Executa a simulacao
	simOut = sim('inversor_full_bridge', 'SimulationMode', 'normal');
    
    % Obtem os valores de vOut
    vout = simOut.get('vout');
    
    % Obtem o maximo valor
    res = max(vout(1, :));

end % Fim da fun??o objetivo
end