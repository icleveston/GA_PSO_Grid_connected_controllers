
from math import pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import itertools
import scipy


class ConverterPlant():


    def __init__(self):

        Lg_2_min = 0e-3  # Indutância    mínima    da    rede
        Lg_2_max = 1e-3  # Indutância    máxima    da    rede

        # Discretização
        fs = 2 * 10020  # Amostragem    na    prática    20040    Hz    e    comutação    10020Hz
        self.Ts = 1 / fs
        Ts = self.Ts
        f_sw = 10020  # Freqüência    de    comutação
        Tsw = 1 / f_sw  # Período    de    comutação
        w_sw = 2 * pi * f_sw

        # Parâmetros    do    filtro
        Lc = 1.e-3  # Indutância    do    lado    do    conversor
        rc = 1e-10  # Resistência    série    do    indutor    do    conversor
        Cap = 62e-6  # Capacitor do filtro LCL
        Lg_1 = 0.3e-3  # Indutância    mínima    da    rede
        rg_1 = 1e-10  # Resistência    série    do    indutor    do    lado    da    rede
        rg_2 = 0.  # Resistência    equivalente    série    da    rede

        Lg_min = Lg_1 + Lg_2_min  # Indutância    TOTAL    da    rede
        Lg_max = Lg_1 + Lg_2_max  # Indutância    TOTAL    da    rede
        rg = rg_1 + rg_2

        # Espaço    de    estados
        # PASSO    1
        Ap_1 = [[-rc / Lc, -1 / Lc, 0], [1 / Cap, 0, -1 / Cap], [0, 1 / Lg_min, -rg / Lg_min]]
        Ap_2 = [[-rc / Lc, - 1 / Lc, 0], [1 / Cap, 0, - 1 / Cap], [0, 1 / Lg_max, - rg / Lg_max]]
        Bp = [[1 / Lc], [0], [0]]
        Fp_1 = [[0], [0], [-1 / Lg_min]]
        Fp_2 = [[0], [0], [-1 / Lg_max]]
        Cp_g = [0, 0, 1]
        Dp = 0

        # Discretização(corrente da    rede)

        # Discretização    ZOH
        # PASSO    2

        x1 = ss(Ap_1, Bp, Cp_g, Dp)
        x2 = ss(Ap_2, Bp, Cp_g, Dp)
        x3 = ss(Ap_1, Fp_1, Cp_g, Dp)
        x4 = ss(Ap_2, Fp_2, Cp_g, Dp)

        Ad_1, Bd_1, Cd_g, Dd = ssdata(matlab.c2d(x1, Ts))
        Ad_2, Bd_2, Cd_g, Dd = ssdata(matlab.c2d(x2, Ts))
        Ad_1, Fd_1, Cd_g, Dd = ssdata(matlab.c2d(x3, Ts))
        Ad_2, Fd_2, Cd_g, Dd = ssdata(matlab.c2d(x4, Ts))

        # CORRENTE  DA REDE = > Inclusão  do atraso  de transporte  em  espaço  de  estados
        # PASSO  3
        Gd_1 = np.concatenate((Ad_1, Bd_1), axis=1)
        Gd_1 = np.concatenate((Gd_1, np.zeros((1, 4))), axis=0)
        Gd_2 = np.concatenate((Ad_2, Bd_2), axis=1)
        Gd_2 = np.concatenate((Gd_2, np.zeros((1, 4))), axis=0)
        Hd = [[0], [0], [0], [1]]
        Hd_dist1 = np.concatenate((Fd_1, [[0]]))
        Hd_dist2 = np.concatenate((Fd_2, [[0]]))
        Cd_grid = np.asarray(np.concatenate((Cd_g, [[0]]), axis=1))
        Dd = 0

        # Controlador   ressonante    fundamental
        # PASSO    4a
        # w = 2 * pi * 60
        #
        # zeta = 1 * 0.0001
        # zeta_3a = 1 * 0.0001
        # zeta_5a = 1 * 0.0001
        # zeta_7a = 1 * 0.0001
        #
        # G_res = control.tf([1, 0], [1, 2*zeta, w**2]) #s ^ 1 / (s ^ 2 + 2 * zeta * s + w ^ 2)
        # G_res_3a = control.tf([1, 0], [1, 2*zeta_3a, 3*w**2]) #s ^ 1 / (s ^ 2 + 2 * zeta_3a * s + (3 * w) ^ 2)
        # G_res_5a = control.tf([1, 0], [1, 2*zeta_5a, 5*w**2]) #s ^ 1 / (s ^ 2 + 2 * zeta_5a * s + (5 * w) ^ 2)
        # G_res_7a = control.tf([1, 0], [1, 2*zeta_7a, 7*w**2]) #s ^ 1 / (s ^ 2 + 2 * zeta_7a * s + (7 * w) ^ 2)
        #
        # #PASSO   4    b
        # G_res_discreto = control.matlab.c2d(G_res, Ts, 'tustin')
        # G_res_discreto_3a = control.matlab.c2d(G_res_3a, Ts, 'tustin')
        # G_res_discreto_5a = control.matlab.c2d(G_res_5a, Ts, 'tustin')
        # G_res_discreto_7a = control.matlab.c2d(G_res_7a, Ts, 'tustin')
        #
        # #Forma   1: Matlab
        # #PASSO 4c
        # R1, T1, U1, V1 = control.ssdata(G_res_discreto)
        # R3, T3, U3, V3 = control.ssdata(G_res_discreto_3a)
        # R5, T5, U5, V5 = control.ssdata(G_res_discreto_5a)
        # R7, T7, U7, V7 = control.ssdata(G_res_discreto_7a)
        #
        #
        # R_a = np.concatenate((R1, np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))))
        # R_b = np.concatenate((np.zeros((2, 2)), R3, np.zeros((2, 2)), np.zeros((2, 2))))
        # R_c = np.concatenate((np.zeros((2, 2)), np.zeros((2, 2)), R5 , np.zeros((2, 2))))
        # R_d = np.concatenate((np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2)), R7))
        #
        # R = np.concatenate((R_a, R_b, R_c, R_d), axis=1)
        #
        # T = np.concatenate((T1, T3, T5, T7))
        # U = [U1, U3, U5, U7]
        # V = V1 + V3 + V5 + V7

        R = np.array([[1.9996,   -1.0000,         0,         0,         0,         0,         0,         0],
    [1.0000  ,       0  ,       0  ,       0  ,       0,         0,         0,         0],
     [    0   ,      0   , 1.9968   ,-1.0000   ,      0,         0,         0,         0],
      [   0    ,     0    ,1.0000    ,     0    ,     0,         0,         0,         0],
       [  0     ,    0     ,    0     ,    0   , 1.9912,   -1.0000,         0,         0],
      [   0      ,   0      ,   0      ,   0    ,1.0000,         0,         0,         0],
        [ 0       ,  0       ,  0       ,  0     ,    0,         0,    1.9827,   -1.0000],
       [  0        , 0        , 0        , 0      ,   0,         0,    1.0000,         0]])



        T = np.array([[0.0078],
             [0],
             [0.0078],
             [0],
             [0.0078],
             [0],
             [0.0078],
             [0]])


        # Ressonante    espaço    de    estados
        # PASSO    5
        A_til_1_a = np.concatenate((Gd_1, np.zeros((4, 8))), axis=1)
        A_til_1_b = np.concatenate((-T*Cd_grid, R), axis=1)
        self.A_til_1 = np.asarray(np.concatenate((A_til_1_a, A_til_1_b)))

        A_til_2_a = np.concatenate((Gd_2, np.zeros((4, 8))), axis=1)
        A_til_2_b = np.concatenate((-T*Cd_grid, R), axis=1)
        self.A_til_2 = np.asarray(np.concatenate((A_til_2_a, A_til_2_b)))

        self.B_til = np.concatenate((Hd, np.zeros((8, 1))))

        Br_til = np.concatenate((np.zeros((4, 1)), T))
        Bdist1_til = np.concatenate((Hd_dist1, np.zeros((8, 1))))
        Bdist2_til = np.concatenate((Hd_dist2, np.zeros((8, 1))))

        self.C_til = np.asarray(np.concatenate((Cd_grid, np.zeros((1, 8))), axis=1))

        # Adi = [A_til_1, A_til_2] # Matriz    dinâmica    do    sistema(x)
        # Bdi = [B_til, B_til] # Matriz    de    controle(u)
        # Edi = [Bdist1_til, Bdist2_til] #Matriz    de    distúrbio(vd)
        # Br = [Br_til, Br_til] #Matriz    de   ref(ref)
        # Ci = [C_til, C_til] #Matriz    de    saída
        # D = [0, 0]
        # Di = [0, 0]
        # Dd = [0, 0]
        #
        # #H2
        # MF1h2 = control.ss(A_til_1 + B_til * Kh2, [Br_til], C_til, 0, Ts) #ig / iref
        # MF2h2 = control.ss(A_til_2 + B_til * Kh2, [Br_til], C_til, 0, Ts) #ig / iref
        #
        # MFdu1h2 = control.ss(A_til_1 + B_til * Kh2, B_til, C_til, 0, Ts) # ig / u
        # MFdu2h2 = control.ss(A_til_2 + B_til * Kh2, B_til, C_til, 0, Ts) # ig / u
        #
        # MFdist1h2 = control.ss(A_til_1 + B_til * Kh2, Bdist1_til, C_til, 0, Ts) # ig / vd
        # MFdist2h2 = control.ss(A_til_2 + B_til * Kh2, Bdist2_til, C_til, 0, Ts) #ig / vd
        #
        # MFx1h2 = control.ss(A_til_1 + B_til * Kh2, A_til_1, C_til, 0, Ts) # ig / dx
        # MFx2h2 = control.ss(A_til_2 + B_til * Kh2, A_til_2, C_til, 0, Ts) #ig / dx

        # Simulação

        to = 0
        tf1 = 0.3 + Ts

        self.t = np.arange(to, tf1, float(Ts))

        # Distúrbio
        disturbio = [1*311*sin(2*pi*60*x) for x in self.t]  # +3 * sin(2 * pi * 180 * t);

        # Referência
        # r = 10 * sin(2 * pi * 60 * t); % 1 * (1 - exp(-(1 / 0.1) * t)). * sin(2 * pi * 60 * t);

        self.ref = []

        for nsample in range(1, len(self.t) + 1):

            if nsample < 2 * 334:
                self.ref.append(0)
            elif nsample >= 2 * 334 and nsample < 4.75 * 334:
                self.ref.append(10 * sin(60 * 2 * pi * Ts * nsample))
            elif nsample >= 4.75 * 334 and nsample < 9 * 334:
                self.ref.append(-10 * sin(60 * 2 * pi * Ts * nsample))
            elif nsample >= 9 * 334 and nsample < 13 * 334:
                self.ref.append(10 * cos(60 * 2 * pi * Ts * nsample))
            elif nsample >= 13 * 334:
                self.ref.append(20 * cos(60 * 2 * pi * Ts * nsample))

        self.u = list(zip(self.ref, disturbio))

        self.B1 = np.asarray(list(zip(Br_til, Bdist1_til))).reshape(-1, 2)
        self.B2 = np.asarray(list(zip(Br_til, Bdist2_til))).reshape(-1, 2)

        d1 = np.diff(self.ref[1000:])

        self.diffRef = len(list(itertools.groupby(d1, lambda x: x > 0)))


    def run(self, Kh2, raio = 10, plot=False):

        if 0 in Kh2:
            return 1e10, 1e5, 1e5, 1e5

        if raio < 1 or plot:

            MF1full2 = ss(self.A_til_1 + self.B_til * Kh2, self.B1, self.C_til, [0, 0], self.Ts)
            MF2full2 = ss(self.A_til_2 + self.B_til * Kh2, self.B2, self.C_til, [0, 0], self.Ts)

            stateSpace1 = scipy.signal.StateSpace(self.A_til_1 + self.B_til * Kh2, self.B1, self.C_til, [0, 0]).to_tf()
            stateSpace2 = scipy.signal.StateSpace(self.A_til_2 + self.B_til * Kh2, self.B2, self.C_til, [0, 0]).to_tf()

            #t, yh21full2, xh21s1 = scipy.signal.lsim(stateSpace1, self.u, self.t, 0, True)
            #t, yh22full2, xh21s1 = scipy.signal.lsim(stateSpace2, self.u, self.t, 0, True)

            yh21full2, t, xh21s1 = matlab.lsim(MF1full2, self.u, self.t)
            yh22full2, t, xh22s1 = matlab.lsim(MF2full2, self.u, self.t)

            # print(sum(yh21full2[0]))
            #
            # import scipy
            # Tt, yout, xout = scipy.signal.lsim((self.A_til_1 + self.B_til * Kh2, self.B1, self.C_til, [0, 0]), self.u, self.t)
            #
            # print(sum(yout))

            if plot:
                return t, yh21full2[0], yh22full2[0], self.ref
            else:

                absYh21 = np.nan_to_num(yh21full2[0])[1000:]
                absYh22 = np.nan_to_num(yh22full2[0])[1000:]

                # e1 = np.transpose(self.ref[1000:6013]) - absYh21
                # e2 = np.transpose(self.ref[1000:6013]) - absYh22
                #
                # ise1 = e1**2
                # ise2 = e2**2
                # resultados1 = sum(ise1)
                # resultados2 = sum(ise2)

                error1 = mean_squared_error(self.ref[1000:], absYh21)
                error2 = mean_squared_error(self.ref[1000:], absYh22)

                # d1 = np.diff(absYh21)
                # d2 = np.diff(absYh22)
                #
                # signalChanges1 = len(list(itertools.groupby(d1, lambda x: x > 0)))
                # signalChanges2 = len(list(itertools.groupby(d2, lambda x: x > 0)))
                #
                # osc1 = abs(signalChanges1 - self.diffRef)
                # ocs2 = abs(signalChanges2 - self.diffRef)

                w1, mag1, phase1 = stateSpace1.bode()
                w2, mag2, phase2 = stateSpace2.bode()

                bode = max(max(mag1), max(mag2))

                if error1 > 1e10:
                    error1 = 1e10

                if error2 > 1e10:
                    error2 = 1e10

                #ISE_Robusto = 1.71e4

                return max(error1, error2), bode

        else:
            return 1e10, 1e5
