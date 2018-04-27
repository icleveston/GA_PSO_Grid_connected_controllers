#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from scipy import io
from matplotlib import cm
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.preprocessing import normalize
import numpy as np
import pickle
from scipy.interpolate import griddata
import matlab.engine
from mpl_toolkits.mplot3d import Axes3D


def normalizeResults(data, attenuation=False):

    # Filter radius < 1
    a = list(filter(lambda x: x['val'][0] < 1, data))

    if attenuation:
        a = filter(lambda x: x['val'][1] < 0.2, a)

    v = [[i['val'][0], i['val'][1], i['val'][2]] for i in a]

    if len(v) > 0:

        normalizedV = normalize(v, axis=0)

        g = 0

        for i in a:

            i['norm'] = sqrt((normalizedV[g][1]) ** 2 + (normalizedV[g][2]) ** 2)

            g += 1

        # Sort by the norm
        return sorted(a, key=lambda x: x['norm'])

    else:
        return []


def plotHistoric(geracao, raio=[], bode=[], ise=[], together=False, title='', savePath=None):

    if together:

        f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=False, sharex=False)

        if len(raio) > 0:

            ax1.plot(geracao, raio)

            # Major ticks every 20, minor ticks every 5
            majorX_ticks = np.arange(0, max(geracao), 20)
            minorX_ticks = np.arange(0, max(geracao), 5)
            majorY_ticks = np.arange(min(raio), max(raio), max(raio)/600)
            minorY_ticks = np.arange(min(raio), max(raio), max(raio)/2400)

            ax1.set_xticks(majorX_ticks)
            ax1.set_xticks(minorX_ticks, minor=True)
            ax1.set_yticks(majorY_ticks)
            ax1.set_yticks(minorY_ticks, minor=True)

            # Or if you want different settings for the grids:
            ax1.grid(which='minor', alpha=0.2)
            ax1.grid(which='major', alpha=0.5)

            ax1.set_xlabel("Geracao")
            ax1.set_ylabel("Raio")
            ax1.set_xlim(0, len(geracao))
            ax1.set_title('Raio')

        if len(bode) > 0:

            ax2.plot(geracao, bode)
            ax2.set_xticks(majorX_ticks)
            ax2.set_xticks(minorX_ticks, minor=True)

            majorY_ticks = np.arange(min(bode), max(bode), max(bode)/60)
            minorY_ticks = np.arange(min(bode), max(bode), max(bode)/240)
            ax2.set_yticks(majorY_ticks)
            ax2.set_yticks(minorY_ticks, minor=True)
            ax2.grid(which='minor', alpha=0.2)
            ax2.grid(which='major', alpha=0.5)
            ax2.set_xlabel("Geracao")
            ax2.set_ylabel("Bode")
            ax2.set_xlim(0, len(geracao))
            ax2.set_title('Bode')

        if len(ise) > 0:

            ise = [x if x < 1e10 else 1.8e5 for x in ise]

            ax3.plot(geracao, ise)
            majorY_ticks = np.arange(min(ise), max(ise), max(ise)/15)
            minorY_ticks = np.arange(min(ise), max(ise), max(ise)/60)
            ax3.set_yticks(majorY_ticks)
            ax3.set_yticks(minorY_ticks, minor=True)
            ax3.set_xticks(majorX_ticks)
            ax3.set_xticks(minorX_ticks, minor=True)
            ax3.grid(which='minor', alpha=0.2)
            ax3.grid(which='major', alpha=0.5)
            ax3.set_xlabel("Geracao")
            ax3.set_ylabel("ISE")
            ax3.set_xlim(0, len(geracao))
            ax3.set_title('ISE')

        f.show()

    else:

        plt.plot(geracao, raio)
        plt.grid()
        plt.xlabel("Geração")
        plt.ylabel("Fitness")
        plt.title(title)

        if savePath is not None:
            plt.savefig(savePath + "/fitness.png", dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                        transparent=False, bbox_inches=None, pad_inches=0.1,
                        frameon=None)
        else:
            plt.show()


# Plot the function 3D
def plotFunction3D(fitness, limitMin=0, limitMax=10, step=0.005, xlabel='x', ylabel='y', zlabel='Custo', title=""):

    from matplotlib.colors import LogNorm

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = y = np.arange(limitMin, limitMax, step)
    X, Y = np.meshgrid(x, y)
    zs = np.array([fitness([x, y])[0] for x, y in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)

    CS = ax.plot_surface(X, Y, Z, rstride=10, cstride=10,  norm=LogNorm(), cmap=cm.jet, linewidth=0.2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    plt.title(title)

    clb = plt.colorbar(CS, orientation='horizontal', shrink=0.5)
    clb.ax.set_title('Custo')

    plt.show()

# Plot function contour
def plotFunctionCountour(fitness, limitMin=0, limitMax=10, step=0.005, all=[], best=[], title="", legend=False):

    if not os.path.exists('Validation/plotContourData.pkl'):

        x = y = np.arange(limitMin, limitMax, step)
        X, Y = np.meshgrid(x, y)
        zs = np.array([fitness([x, y])[0] for x, y in zip(np.ravel(X), np.ravel(Y))])
        Z = zs.reshape(X.shape)

        output = open('Validation/plotContourData.pkl', 'wb')
        pickle.dump([X, Y, Z], output)
        output.close()

    else:
        file = open('Validation/plotContourData.pkl', 'rb')

        data = pickle.load(file)

        X = data[0]
        Y = data[1]
        Z = data[2]

    fig = plt.figure()

    if all:

        i = 0
        for e in all:

            k1 = [part[0] for part in e['pop']]
            k2 = [part[1] for part in e['pop']]

            ax = fig.add_subplot(2, 3, i + 1)
            ax.set_title(str(title[i]))
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_xlim(0, 10)
            ax.set_ylim(0, 10)

            ax.contour(X, Y, Z, cmap=cm.coolwarm)
            ind = ax.scatter(k1, k2)

            if not isinstance(e['best'], bool):
                best = ax.scatter(e['best'][0], e['best'][1], c='red')
                fig.legend([ind, best], [legend, 'Melhor Posição do Enxame'], ncol=2, loc=8)
            else:
                fig.legend([ind], [legend], ncol=1, loc=8)

            i = i + 1

        plt.show()

    else:

        CS = plt.contour(X, Y, Z, cmap=cm.coolwarm)

        if legend:
            plt.clabel(CS, fontsize=9, inline=1)

        if best:
            plt.scatter(best[0], best[1], c='red')

        if title != "":
            plt.title(title)

        plt.xlim(-500, 500)
        plt.ylim(-500, 500)
        plt.show()


# Plot the Pareto Frontier 3D
def plotPareto3D(path):

    file = open(path + '/results.pkl', 'rb')

    data = pickle.load(file)

    b = normalizeResults(data, True)

    if len(b) > 0:

        #Plot Pareto Frontier
        paretoX = []
        paretoY = []
        paretoZ = []

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in b:

            if i['val'][0] == b[0]['val'][0]:
                continue

            paretoX.append(i['val'][0])
            paretoY.append(i['val'][1])
            paretoZ.append(i['val'][2])

        io.savemat('pareto.mat', {'pareto': [paretoY, paretoZ, paretoX]})

        exit()

        # re-create the 2D-arrays
        x1 = np.linspace(min(paretoX), max(paretoX), len(paretoX))
        y1 = np.linspace(min(paretoY), max(paretoY), len(paretoY))
        x2, y2 = np.meshgrid(x1, y1)
        z2 = griddata((paretoX, paretoY), paretoZ, (x2, y2), method='cubic')

        # Plot the frontier
        #ax.plot_surface(y2, z2, x2, rstride=1, cstride=1, linewidth=0, antialiased=False)
        #ax.plot_trisurf(paretoY, paretoZ, paretoX, edgecolor='none')
        ax.plot_wireframe(y2, z2, x2)
        #ax.scatter(paretoY, paretoZ, paretoX)

        #Plot the best gene
        ax.scatter(b[0]['val'][1], b[0]['val'][2], b[0]['val'][0], c='red')

        print(b[0]['ind'])
        print(b[0]['val'])

        ax.set_xlabel('Bode')
        ax.set_ylabel('ISE')
        ax.set_zlabel('Raio')

        plt.show()


# Plot the Pareto Frontier 2D
def plotPareto2D(path):

    file = open(path + '/results.pkl', 'rb')

    data = pickle.load(file)

    print(len(data))

    # Normalize data
    a = normalizeResults(data)

    if len(a) > 0:

        print(len(a))
        print(a[0]['ind'])
        print(evaluate(a[0]['ind']))

        #Plot Pareto Frontier
        paretoX = []
        paretoY = []

        for i in a:
            if i['val'][0] == a[0]['val'][0]:
                continue

            #paretoX.append(i['normX'])
            #paretoY.append(i['normY'])
            paretoX.append(i['val'][0])
            paretoY.append(i['val'][1])

        #Plot the frontier
        plt.scatter(paretoY, paretoX)

        #Plot the best gene
        #plt.scatter(a[0]['normY'], a[0]['normX'], c='red')
        plt.scatter(a[0]['val'][1], a[0]['val'][0], c='red')

        plt.xlabel('Bode')
        plt.ylabel('Raio')

        #plt.xlim(0, 0.2)
        #plt.ylim(0, 0.2)
        #plt.gca().set_aspect('equal', adjustable='box')

        plt.show()


#Print stats for every simulation
def plotStats(path):

    file = open(path + '/paretoFrontier.pkl', 'rb')

    data = pickle.load(file)

    for i in data:

        results = []

        for ind in i:

            res = {
                'ind': list(ind),
                'val': evaluate(list(ind))
            }

            results.append(res)

        # Normalize data
        a = normalizeResults(results)

        if len(a) > 0:

            print(a[0]['ind'])
            print(evaluate(a[0]['ind']))

            plot(a[0]['ind'])


# Plot pareto for each generation
def plotParetoForEachGeneration(path):

    file = open(path + '/paretoFrontier.pkl', 'rb')

    data = pickle.load(file)

    resultsGeneration = []

    for i in data:

        results = []

        for ind in i:
            res = {
                'ind': list(ind),
                'val': evaluate(list(ind))
            }

            results.append(res)

        resultsGeneration.append(min(results, key=lambda x:x['val'][0]))


    geracao = list(range(1, len(resultsGeneration) + 1))

    raio = []
    bode =[]
    ise = []

    count = 1
    for i in resultsGeneration:

        # Print the data
        print(str(count) + "\t\t" + str(i))
        count = count + 1

        # Add to the array
        raio.append(i['val'][0])

        if len(i['val']) >= 2:
            bode.append(i['val'][1])

        if len(i['val']) >= 3:
            ise.append(i['val'][2])

    print("Total de Gerações: " + str(len(geracao)))

    # Plot the graphs
    plotGeracao(geracao, raio, bode, ise)


def correlation(path):

    file = open(path + '/results.pkl', 'rb')

    data = pickle.load(file)

    print(len(data))

    # Normalize data
   #a = normalizeResults(data)
    a = filter(lambda x: x['val'][0] < 1, data)
    a = sorted(a, key=lambda x: x['val'][2])

    raio = []
    bode = []
    ise = []

    for i in a:

        raio.append(i['val'][0])
        bode.append(i['val'][1])
        ise.append(i['val'][2])


    geracao = list(range(1, len(a) + 1))

    # Plot the graphs
    plotHistoric(geracao, raio, bode, ise)


def plotConversor(solution, eng, c):

    a = matlab.double(solution)

    t, y1, y2, ref = eng.plot(c, a, nargout=4)

    f = plt.figure()
    ax = f.add_subplot(2, 1, 1)
    ax.plot(t, y1)
    ax.plot(t, ref)
    ax.grid()
    ax.set_xlabel("Tempo de Simulação (ms)")
    ax.set_ylabel("Tensão (V)")
    ax.set_title("Indutância Mínima")

    ax = f.add_subplot(2, 1, 2)
    ax.plot(t, y2)
    ax.plot(t, ref)
    ax.grid()
    ax.set_xlabel("Tempo de Simulação (ms)")
    ax.set_ylabel("Tensão (V)")
    ax.set_title("Indutância Máxima")
    f.legend(['Sinal de Controle', 'Sinal de Referência'], ncol=2, loc=8)
    plt.show()
