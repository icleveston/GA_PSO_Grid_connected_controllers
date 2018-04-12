function resultados = otimiza_dlqr
A_til_1 = 0;
A_til_2 = 0;
B_til = 0;
Br_til = 0;
Bdist1_til = 0;
Bdist2_til = 0;
C_til = 0;
Ts = 0;

load('parametros_projeto.mat', 'A_til_1', 'A_til_2', 'B_til', 'Bdist1_til', 'Bdist2_til', 'Br_til', 'C_til', 'Ts')

%ganhos H2 - produz menor sobrecorrente na partida
Kfull2=[-9.747423809401397  -1.851225100912732  -0.350767371688946  -0.446663368860505  50.868719652179720 -50.184356792361314  16.232137022975223 -15.756612317015055   8.734346658571226  -8.843404878241969 3.712318389421682  -4.140705855357965];
%Kfull2 = -dlqr(A_til_1, B_til, diag([1 1 1 1 50*ones(1,8)]), 10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simula��o
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

to = 0;
tf = 0.1;

t = to:Ts:tf;       
[n1 n2] = size(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dist�rbio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disturbio = 1*(311)*sin(2*pi*60*t);%+3*sin(2*pi*180*t); 

for nsample = 1:n2
           
          if (nsample < 2*334)
      	       ref(nsample) = 0;
          end
          if (nsample >= 2*334 && nsample < 4.75*334)
      	       ref(nsample) = 10.*sin(60*2*pi*Ts*nsample);
          end
          if (nsample >= 4.75*334 && nsample < 9*334)
      	       ref(nsample) = -10.*sin(60*2*pi*Ts*nsample);
          end
          if (nsample >= 9*334 && nsample < 13*334)
      	       ref(nsample) = 10.*cos(60*2*pi*Ts*nsample);
          end
          if (nsample >= 13*334)
      	       ref(nsample) = 20.*cos(60*2*pi*Ts*nsample);
          end

end

u = [ref; disturbio;];

MF1full2 = ss(A_til_1 + B_til*Kfull2, [Br_til Bdist1_til], C_til, 0, Ts ); 
MF2full2 = ss(A_til_2 + B_til*Kfull2, [Br_til Bdist2_til], C_til, 0, Ts ); 

%%
% N�mero de vari�veis
nvar = 8;
% Limites
limsup = 100*[1 1 1 1 1 1 1 1];
liminf = -100*[1 1 1 1 1 1 1 1];

% Configura��es do ga
%setupga = gaoptimset(...
%	'EliteCount', 0,...
%	'CrossoverFraction', 0.5,...
%	'Generations', 100,...
%	'StallGenLimit', 20,...
%	'Display', 'iter',...
%	'PopulationSize', 100,...
%	'SelectionFcn', @selectiontournament,...
%	'CrossoverFcn', @crossoverscattered, ...
%	'MutationFcn', @mutationadaptfeasible);
%
%% Executa o GA
x0 = [1 1 1 1 1 1 1 1];
%ga(@(x)fobj(x), nvar,...
%	[], [], ... % Desigualdades lineares
%	[], [],... % Igualdades
%	liminf, limsup,... % Limites de busca
%	[],... % Restri��es
%setupga);

resultados = x0;

plot(lsim(MF1full2, u, t, [0 0 0 0 1 1 1 1 1 1 1 1]));
%hold on
%plot(lsim(MF1full2, u, t, [0 0 0 0 0*x0]),'r')

function res = fobj(genes)
	% Simula o sistema
	yh21full1 =  lsim(MF1full2, u, t, [0 0 0 0 genes]);
	yh21full2 =  lsim(MF2full2, u, t, [0 0 0 0 genes]);
  
	% Fitness function
	res = max([abs(yh21full1); abs(yh21full2)]);
end 
end