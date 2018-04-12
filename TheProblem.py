
from platypus import *
import matlab.engine
from math import pi, sin, cos, asin
import numpy as np
from scipy import io


class TheProblem(Problem):

    def __init__(self, eng, c):

        self.nobjs = 3
        self.nvars = 12

        self.eng = eng
        self.c = c

        #Initialize the parent
        super(TheProblem, self).__init__(self.nvars, self.nobjs)

        self.types[0] = Real(-15, 0)
        self.types[1] = Real(-15, 0)
        self.types[2] = Real(-15, 0)
        self.types[3] = Real(-15, 0)
        self.types[4] = Real(0, 100)
        self.types[5] = Real(-100, 0)
        self.types[6] = Real(0, 50)
        self.types[7] = Real(-50, 0)
        self.types[8] = Real(0, 50)
        self.types[9] = Real(-50, 0)
        self.types[10] = Real(0, 50)
        self.types[11] = Real(-50, 0)


    def evaluate(self, solution):

        a = matlab.double(solution.variables)

        res = self.eng.testa(self.c, a, nargout=1)

        solution.objectives[:] = [res[0][0], res[0][1], res[0][2]]

    def evaluateRaio(self, solution):

        # Import the parameters
        mat = io.loadmat('parametros_projeto.mat')

        # Get the variables
        A_til_1 = mat['A_til_1']
        A_til_2 = mat['A_til_2']
        B_til = mat['B_til']
        Bdist1_til = mat['Bdist1_til']
        Br_til = mat['Br_til']
        C_til = mat['C_til']
        Ts = mat['Ts']

        MF1full2 = A_til_1 + B_til * solution
        MF2full2 = A_til_2 + B_til * solution

        ai = np.concatenate((MF1full2, MF2full2), axis=1)

        return self._nuvem(ai)[0],


    def _nuvem(self, ai):

        (nAi, mAi) = np.shape(ai)
        vertices = int(mAi / nAi)
        max_avl = -10000
        eig_ca = []
        eig_cb = []

        for i in range(1, vertices):
            for j in range(i + 1, vertices+1):
                for alfa in np.arange(0, 1.005, 0.005):

                    aux1 = alfa * ai[:, nAi*(i-1): i*nAi] + (1 - alfa) * ai[:, nAi*(j-1): j*nAi]
                    avls = np.linalg.eig(aux1)[0].reshape(-1, 1)

                    #np.concatenate((eig_ca, avls))

                    max_eig = max(abs(avls))

                    if max_eig > max_avl:
                        max_avl = max_eig

        for k in range(1, 1001):
            i = 1
            j = 2

            comb = np.random.rand(1, 2)[0]
            comb = list(map(lambda x: x/sum(comb), comb))

            aux1 = comb[0] * ai[:, nAi*(i-1):i * nAi] + comb[1] * ai[:, nAi*(j-1):j*nAi]
            avls = np.linalg.eig(aux1)[0].reshape(-1, 1)

            #eig_cb = zip(eig_cb, avls)

            max_eig = max(abs(avls))

            if max_eig > max_avl:
                max_avl = max_eig

        return max_avl


    def _itae(self, K, plot=False):

        wn = 2 * pi * 60
        csi = 0.001
        R = [0, 1, -wn ** 2 - 2 * csi * wn]
        S = np.array([[0], [1]])
        Ts = 1 / 10000  # Frequencia deamostragem !!!

        # Implementacao mais favoravel do ressonante
        w = 2 * pi * 60
        zeta = 1 * 0.0001

        # G_res = control.tf([w**2], [1, 2*zeta, w**2]) #Trocamos o s do numerador por w ^ 2

        # G_res_discreto = control.matlab.c2d(G_res, Ts, 'tustin')

        # (Rd, Sd, U1, V1) = control.ssdata(G_res)

        Rd = np.array([[1.19986, -1.0000], [1.0000, 0]])
        Sd = np.array([[0.0313], [0]])

        # Par�metros para projeto robusto - ---------------------------------------

        R1 = 0.1
        Lg_nom = 5e-3
        delta_Lg = 3e-3
        delta_Rg = 0.1

        Lmin = Lg_nom - delta_Lg
        Lmax = Lg_nom + delta_Lg
        Rlmin = R1 - delta_Rg
        Rlmax = R1 + delta_Rg

        C = [1, 0]
        D = 0
        a11 = (-Rlmin / Lmin) * Ts + 1
        b11 = Ts / Lmin

        a12 = (-Rlmax / Lmin) * Ts + 1
        b12 = Ts / Lmin

        a21 = (-Rlmin / Lmax) * Ts + 1
        b21 = Ts / Lmax

        a22 = (-Rlmax / Lmax) * Ts + 1
        b22 = Ts / Lmax

        A1 = np.concatenate((([[a11, b11, 0, 0], [0, 0, 0, 0]]), np.concatenate((-Sd, np.zeros((2, 1)), Rd), axis=1)))

        B1 = np.array([[0], [1], [0], [0]])
        Br1 = np.concatenate(([[0], [0]], Sd))
        Bd1 = np.array([[-b11], [0], [0], [0]])

        A2 = np.concatenate((([[a12, b12, 0, 0], [0, 0, 0, 0]]), np.concatenate((-Sd, np.zeros((2, 1)), Rd), axis=1)))
        B2 = np.array([[0], [1], [0], [0]])
        Br2 = np.concatenate(([[0], [0]], Sd))
        Bd2 = np.array([[-b12], [0], [0], [0]])

        A3 = np.concatenate((([[a21, b21, 0, 0], [0, 0, 0, 0]]), np.concatenate((-Sd, np.zeros((2, 1)), Rd), axis=1)))
        B3 = np.array([[0], [1], [0], [0]])
        Br3 = np.concatenate(([[0], [0]], Sd))
        Bd3 = np.array([[-b21], [0], [0], [0]])

        A4 = np.concatenate((([[a22, b22, 0, 0], [0, 0, 0, 0]]), np.concatenate((-Sd, np.zeros((2, 1)), Rd), axis=1)))
        B4 = np.array([[0], [1], [0], [0]])
        Br4 = np.concatenate(([[0], [0]], Sd))
        Bd4 = np.array([[-b22], [0], [0], [0]])

        Ai = np.array([A1, A2, A3, A4])
        Bi = np.array([B1, B2, B3, B4])
        Bri = np.array([Br1, Br2, Br3, Br4])
        Bdi = np.array([Bd1, Bd2, Bd3, Bd4])

        # r = 1 % Raio para  aloca��o   de   polos
        # out = ssf_stab_K_d_mtb(Ai / r, Bi / r)
        # K = out.K;

        # $$$$$$$$$$$$$$$$$$$$$$$$ FIM DO PROJETO ROBUSTO  $$$$$$$$$$$$$$$$$$$$$$$

        # ------------------------------------------------------------------------
        # PARAMETROS DE SIMULACAO
        # ------------------------------------------------------------------------

        Cic = 2  # N�mero  de  c�clos    de   rede    simulados
        fr = 60  # Frequ�ncia de sa�da
        Tr = 1 / fr  # Per�odo da tens�o da rede
        fs = 10000  # Frequencia de amostragem
        Ts = 1 / fs
        fsw = 10000
        Tsw = 1 / fsw

        # Tempo
        PerT = Cic * Tr  # Per�odo de simula��o
        dT = 1E-6  # Passo     da    simula��o

        t = np.transpose([x for x in np.arange(0, PerT, dT)])  # Vetor do tempo total da simula��o
        # tsim = np.transpose([x for x in np.arange(0, PerT, dT)])

        # Pontos
        Pontos_fs = 1 / (dT * fs)  # Pontos da simula��o(por ciclo de Ts)
        # Pontos_fsw = 1 / (dT * fsw)
        # Pontos_fr = np.floor(Pontos_fs * Tr / Ts) #Pontos em um ciclo da rede
        Pontos_t = len(t)  # Total de pontos

        # Par�metros do sistema
        Vcc = 200  # Tens�o do barramento CC(V)
        Lg = 2e-3  # Indut�ncia da rede inicial(valor real)
        Lg2 = 2e-3  # Indut�ncia da rede ap�s varia��o(valor real)
        Lg_nom = 5e-3  # Indut�ncia  da rede(valor nominal - utilizado para projeto)
        Rf = 0.1  # Resist�ncia do filtro de sa�da(valor real)
        Rf2 = 0.1
        Rf_nom = 0.1  # Resist�ncia do filtro de sa�da(valor nominal - projeto)
        vg_pk = 127 * np.sqrt(2)  # Tens�o da rede(disturbio)
        # Ma = vg_pk / Vcc # �ndice de modula��o de amplitude
        w = 2 * pi * fr  # Frequ�ncia angular
        ig_ref_pk = 10  # Corrente de referencia(peak)

        # Inicializa��es
        t_k = 0
        t_ks = 0
        vtr_tk = np.zeros((Pontos_t, 1))

        upwm = []
        x = {}
        theta = {}
        rho1 = {}
        rho2 = {}
        xc = {}
        ref = {}
        ref_k = {}
        u = {}
        x[0] = 0
        theta[0] = 0
        rho1[0] = 0
        rho2[0] = 0
        xc[0] = 0
        u[0] = 0
        upwm_k = 0
        ref[0] = 0
        cont = 1
        u_ks = 0
        ref_ks = 0
        ig_amost = 0
        ref_amost = 0
        vtr_ref_amost = {}
        vtr_u_amost = {}
        vtr_ig_amost = {}
        vg = {}

        ks = 0

        # Modelo do conversor emespa�o de estados

        # Planta utilizada para projeto(valores nominais)
        # a11 = (-Rf / Lg_nom) * Ts + 1
        # b11 = Ts / Lg_nom

        # A1 = np.concatenate((([[a11, b11, 0, 0], [0, 0, 0, 0]]), np.concatenate((-Sd, np.zeros((2, 1)), Rd), axis=1)))
        # B1 = np.array([[0], [1], [0], [0]])

        for k in range(0, Pontos_t):  # k � o tempo "continuo"

            if t_k > 0.035:  # 0.03474 % 0.020833:
                # Planta real(continua)
                an = (-Rf2 / Lg2) * (dT) + 1
                bn = dT / Lg2
            else:
                # Planta real(continua)
                an = (-Rf / Lg) * (dT) + 1
                bn = dT / Lg

            ref_k[k] = ig_ref_pk * sin(w * t_k)  # corrente de referencia em t_k

            # Amostragem, ks � o tempo discreto
            if (k % np.floor(Pontos_fs) == 0):
                ref[ks] = ig_ref_pk * sin(w * t_ks)  # corrente    de   referencia
                u[ks] = K[0] * x[ks] + K[1] * theta[ks] + K[2] * rho1[ks] + K[3] * rho2[ks]
                u_ks = u[ks]

                x[ks + 1] = xc[k]
                theta[ks + 1] = u[ks]
                rho1[ks + 1] = -Sd[0, 0] * xc[k] + 0 * theta[ks] + Rd[0, 0] * rho1[ks] + Rd[0, 1] * rho2[ks] + 0 * u[
                    ks] + Sd[0, 0] * ref[ks]
                rho2[ks + 1] = -Sd[1, 0] * xc[k] + 0 * theta[ks] + Rd[1, 0] * rho1[ks] + Rd[1, 1] * rho2[ks] + 0 * u[
                    ks] + Sd[1, 0] * ref[ks]

                ks = ks + 1
                t_ks = t_ks + Ts

            # Modula��o phase - shift - --------------------------------------
            # Nesta t�cnica, frequencia efetiva = 2 * fsw.

            v_tri = 2 * asin(sin(2 * pi * t_k / Tsw - pi / 2)) / pi

            if (u_ks / Vcc > v_tri):
                sa = 1
            else:
                sa = 0

            if (-u_ks / Vcc > v_tri):
                sb = 1
            else:
                sb = 0

            upwm_k = (sa - sb) * Vcc

            # Disturbio - tensao da rede
            vg[k] = vg_pk * sin(w * t_k)
            vtr_tk[k] = t_k

            # Modelo do Conversor(real)
            xc[k + 1] = an * xc[k] + bn * upwm_k - bn * vg[k]

            t_k = t_k + dT

        itae = 0

        for cont in range(0, 8334 * 2):
            itae = itae + cont * abs(ref_k[cont] - xc[cont])

        return itae