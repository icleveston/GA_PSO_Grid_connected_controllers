

import numpy as np


def nuvem(ai):

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