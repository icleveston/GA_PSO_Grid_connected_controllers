
xRobust = [-9.7474, -1.8512, -0.3507, -0.4466, 50.8687, -50.1843, 16.2321, -15.7566, 8.7343, -8.8434, 3.7123, -4.1407];

% Limites
limsup = xRobust + 10;
liminf = xRobust - 10;

options = gaoptimset('PopulationSize',100, 'HybridFcn', @fgoalattain, 'display', 'iter', 'UseParallel', false);

%Create the object 
o = Controle;

%Define a fitness function
fitness = @(x)o.testa(x);

% Executa o GA
[x, fval] = gamultiobj(fitness, 12,[],[],[],[],liminf,limsup, [], options);

[val, indexMin] = find(fval(:, 2) < 0.05 & fval(:, 3) < 1);

fval(indexMin, :)

x(indexMin, :)