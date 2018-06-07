classdef Controle < handle
   properties
     ref,
     u,
     t,
     Ts,
     A_til_1,
     A_til_2,
     B_til,
     Br_til,
     Bdist1_til,
     Bdist2_til,
     C_til,
     n2
   end
   methods
       
       function obj = Controle()

            set(0,'DefaultFigureWindowStyle','docked') 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Parâmetros para o projeto
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Parâmetros Rede

            Lg_2_min = 0.5e-3;   % Indutância mínima da rede
            Lg_2_max = 1e-3;   % Indutância máxima da rede

            % Discretização
            fs    = 2*10020; %Amostragem na prática 20040 Hz e comutação 10020 Hz; 
            obj.Ts    = 1/fs;
            s     = tf('s'); 
            z     = tf('z',obj.Ts);
            f_sw  = 10020;         % Freqüência de comutação
            Tsw   = 1/f_sw;        % Período de comutação
            w_sw  = 2*pi*f_sw;     %

            % Parâmetros do filtro
            Lc   = 1.e-3;        % Indutância do lado do conversor
            rc   = 1e-10;       % Resistência série do indutor do conversor
            Cap   = 62e-6;       % Capacitor do filtro LCL
            Lg_1 = 0.3e-3;      % Indutância mínima da rede
            rg_1 = 1e-10;       % Resistência série do indutor do lado da rede
            rg_2 = 0.;          % Resistência equivalente série da rede

            Lg_min = Lg_1+Lg_2_min;       % Indutância TOTAL da rede
            Lg_max = Lg_1+Lg_2_max;       % Indutância TOTAL da rede
            rg = rg_1+rg_2;

            %% Espaço de estados
            % PASSO 1
            Ap_1 = [-rc/Lc -1/Lc 0; 1/Cap 0 -1/Cap; 0 1/Lg_min -rg/Lg_min];
            Ap_2 = [-rc/Lc -1/Lc 0; 1/Cap 0 -1/Cap; 0 1/Lg_max -rg/Lg_max];
            Bp   = [1/Lc; 0; 0]; 
            Fp_1 = [0; 0; -1/Lg_min]; 
            Fp_2 = [0; 0; -1/Lg_max]; 
            Cp_g = [0 0 1];  
            Dp   = 0;

            % Discretização (corrente da rede)

            % Discretização ZOH
            % PASSO 2
            [Ad_1,Bd_1,Cd_g,Dd]   = ssdata(c2d(ss(Ap_1,Bp,Cp_g,Dp),obj.Ts));
            [Ad_2,Bd_2,Cd_g,Dd]   = ssdata(c2d(ss(Ap_2,Bp,Cp_g,Dp),obj.Ts));
            [Ad_1,Fd_1,Cd_g,Dd]   = ssdata(c2d(ss(Ap_1,Fp_1,Cp_g,Dp),obj.Ts));
            [Ad_2,Fd_2,Cd_g,Dd]   = ssdata(c2d(ss(Ap_2,Fp_2,Cp_g,Dp),obj.Ts));

            %% CORRENTE DA REDE => Inclusão do atraso de transporte em espaço de estados
            % PASSO 3
            Gd_1     = [Ad_1 Bd_1; 0 0 0 0]; 
            Gd_2     = [Ad_2 Bd_2; 0 0 0 0];
            Hd       = [0; 0; 0; 1];         
            Hd_dist1 = [Fd_1; 0];             
            Hd_dist2 = [Fd_2; 0];
            Cd_grid  = [Cd_g 0];
            Dd       = 0;


            %% Controlador ressonante fundamental
            % PASSO 4a
            w = 2*pi*60;

            zeta    = 1*0.0001;
            zeta_3a = 1*0.0001;
            zeta_5a = 1*0.0001;
            zeta_7a = 1*0.0001;

            G_res    = s^1/(s^2 + 2*zeta*s + w^2);
            G_res_3a = s^1/(s^2 + 2*zeta_3a*s + (3*w)^2);
            G_res_5a = s^1/(s^2 + 2*zeta_5a*s + (5*w)^2);
            G_res_7a = s^1/(s^2 + 2*zeta_7a*s + (7*w)^2);

            % PASSO 4b
            G_res_discreto    = c2d(G_res,obj.Ts,'tustin');
            G_res_discreto_3a = c2d(G_res_3a,obj.Ts,'tustin');
            G_res_discreto_5a = c2d(G_res_5a,obj.Ts,'tustin');
            G_res_discreto_7a = c2d(G_res_7a,obj.Ts,'tustin');


            % Forma 1: Matlab
            % PASSO 4c
            [R1,T1,U1,V1] = ssdata(G_res_discreto);
            [R3,T3,U3,V3] = ssdata(G_res_discreto_3a);
            [R5,T5,U5,V5] = ssdata(G_res_discreto_5a);
            [R7,T7,U7,V7] = ssdata(G_res_discreto_7a);

            R = [R1 zeros(2,2) zeros(2,2) zeros(2,2);
                 zeros(2,2) R3 zeros(2,2) zeros(2,2);
                 zeros(2,2) zeros(2,2) R5 zeros(2,2);
                 zeros(2,2) zeros(2,2) zeros(2,2) R7];
            T = [T1;T3;T5;T7];
            U = [U1 U3 U5 U7];
            V = V1+V3+V5+V7;


            % Ressonante espaço de estados
            % PASSO 5
            obj.A_til_1 = [  Gd_1      zeros(4,8);
                       -T*Cd_grid      R];
            obj.A_til_2 = [  Gd_2      zeros(4,8);
                       -T*Cd_grid      R];
            obj.B_til = [Hd;
                    zeros(8,1)];

            obj.Br_til = [zeros(4,1);
                     T];
            obj.Bdist1_til = [Hd_dist1;zeros(8,1)];
            obj.Bdist2_til = [Hd_dist2;zeros(8,1)];

            obj.C_til = [Cd_grid zeros(1,8)];



            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Simulação
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            to = 0;
            tf1 = 0.3;

            obj.t = to:obj.Ts:tf1;       
            [n1 n2] = size(obj.t);
            
            obj.n2 = n2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Distúrbio
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disturbio = 1*(311)*sin(2*pi*60*obj.t);%+3*sin(2*pi*180*t); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Referência
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for nsample = 1:n2

                  if (nsample < 2*334)
                       obj.ref(nsample) = 0;
                  end
                  if (nsample >= 2*334 && nsample < 4.75*334)
                       obj.ref(nsample) = 10.*sin(60*2*pi*obj.Ts*nsample);
                  end
                  if (nsample >= 4.75*334 && nsample < 9*334)
                       obj.ref(nsample) = -10.*sin(60*2*pi*obj.Ts*nsample);
                  end
                  if (nsample >= 9*334 && nsample < 13*334)
                       obj.ref(nsample) = 10.*cos(60*2*pi*obj.Ts*nsample);
                  end
                  if (nsample >= 13*334)
                       obj.ref(nsample) = 20.*cos(60*2*pi*obj.Ts*nsample);
                  end

            end

            obj.u = [obj.ref; disturbio;];

       end;
       
       function res = testa(obj, K)
                                
            raio = obj.nuvem(K);
                        
            MFdu1h2 = ss(obj.A_til_1 + obj.B_til*K, obj.B_til, obj.C_til, 0, obj.Ts ); %ig/u
            MFdu2h2 = ss(obj.A_til_2 + obj.B_til*K, obj.B_til, obj.C_til, 0, obj.Ts ); %ig/u
            
            [mag1,phase1]=bode(MFdu1h2);
            [mag2,phase2]=bode(MFdu2h2);

            bodeRes = max(max(mag1),max(mag2));

            %if raio < 1
            %
            %    MF1full2 = ss(obj.A_til_1 + obj.B_til*K, [obj.Br_til obj.Bdist1_til], obj.C_til, 0, obj.Ts);
            %    MF2full2 = ss(obj.A_til_2 + obj.B_til*K, [obj.Br_til obj.Bdist2_til], obj.C_til, 0, obj.Ts);

            %    [yh21full2,t,xh21s1] = lsim(MF1full2,obj.u,obj.t,'b');
            %    [yh22full2,t,xh22s1] = lsim(MF2full2,obj.u,obj.t,'c');

            %    e1 = obj.ref(1000:6013)'-yh21full2(1000:6013);
            %    e2 = obj.ref(1000:6013)'-yh22full2(1000:6013);

            %    ise1 = e1'*e1;
            %    ise2 = e2'*e2;
                
            %    ise = max(ise1,ise2);
               
            %else
                ise = 1e10;
            %end

            %k=1;
            %k2=1;
            %x1(12,1)=0;
            %y1(1)=0;
            %x2(12,1)=0;
            %y2(1)=0;
            %tp(1)=0;

            %for nsample = 1:obj.n2

            %    disturb(k) = 1*(311)*sin(2*pi*60*k*obj.Ts);
            %    u1(k)=K*x1(:,k);
            %    x1(:,k+1)=obj.A_til_1*x1(:,k)+obj.B_til*u1(k)+obj.Bdist1_til*disturb(k)+obj.Br_til*obj.ref(nsample);
            %    y1(k)=obj.C_til*x1(:,k);
            %    k=k+1;
            %end

            %for nsample = 1:obj.n2
            %    u2(k2)=K*x2(:,k2);
            %    x2(:,k2+1)=obj.A_til_2*x2(:,k2)+obj.B_til*u2(k2)+obj.Bdist2_til*disturb(k2)+obj.Br_til*obj.ref(nsample);
            %    y2(k2)=obj.C_til*x2(:,k2);
            %    k2=k2+1;
            %end

            %for nsample = 1:obj.n2-1
            %    derivada1(nsample) = u1(nsample+1)-u1(nsample);
            %    derivada2(nsample) = u2(nsample+1)-u2(nsample);
            %end

            %max_derivada = max(max(abs(derivada1)), max(abs(derivada2)));

            %if max_derivada > 1e10
            %    max_derivada = 1e10;
            %end

            res = [raio bodeRes ise];

       end
       
       function max_avl = nuvem(obj, K)

            Ai=[obj.A_til_1 + obj.B_til*K obj.A_til_2 + obj.B_til*K];

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
       end
        
       function [t, y1, y2, ref] = plot(obj, K)
           
            MF1full2 = ss(obj.A_til_1 + obj.B_til*K, [obj.Br_til obj.Bdist1_til], obj.C_til, 0, obj.Ts); 
            MF2full2 = ss(obj.A_til_2 + obj.B_til*K, [obj.Br_til obj.Bdist2_til], obj.C_til, 0, obj.Ts);

            [yh21full2,t,xh21s1] = lsim(MF1full2,obj.u,obj.t,'b');
            [yh22full2,t,xh22s1] = lsim(MF2full2,obj.u,obj.t,'c');
           
            t = obj.t';
            y1 = yh21full2;
            y2 = yh22full2;
            ref = obj.ref';
       end
       
    end
end

