%=================================================================
%Funcao que calcula nuvem do politopo Ai discreto com grid fino!
%
function [max_avl] = nuvem_d(Ai,desenha)

[nAi, mAi] = size(Ai);
vertices = mAi/nAi;
max_avl = -10000;
eig_ca=[];
for i = 1:vertices-1
    for j = i+1:vertices
        for alfa=0:0.005:1%Combinacao convexa dois a dois
            aux1 = alfa*Ai(:,nAi*(i-1)+1:i*nAi) + (1-alfa)*Ai(:,nAi*(j-1)+1:j*nAi);
            avls = eig(aux1);
            eig_ca=[eig_ca;avls];
            max_eig = max(abs(avls));%Pega o maior autovalor
            if  max_eig > max_avl
                max_avl = max_eig;%Vai guardando o maior autovalor
            end
        end
    end
end
eig_cb=[];
for k = 1:1000;
    i = 1;
    j = 2;
    comb = rand(1,2);
    comb = comb/sum(comb');
    aux1 = comb(1)*Ai(:,nAi*(i-1)+1:i*nAi) + comb(2)*Ai(:,nAi*(j-1)+1:j*nAi);
    avls = eig(aux1);
    eig_cb = [eig_cb;avls];
    max_eig = max(abs(avls));%Pega o maior autovalor
    if  max_eig > max_avl
        max_avl = max_eig;%Vai guardando o maior autovalor
    end
end

%if nargin == 2
%    plot(real(eig_ca(1:size(eig_ca,1))),imag(eig_ca(1:size(eig_ca,1))),'.b')
%    hold on;
%    plot(real(eig_cb(1:size(eig_cb,1))),imag(eig_cb(1:size(eig_cb,1))),'.g')
%end
