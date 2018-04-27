#!/usr/bin/env python
# -*- coding: utf-8 -*-

from helpers import *
import os.path
from ga import *
from pso import *
from platypus import *
from datetime import datetime
from tabulate import tabulate
from Benchmark import Benchmark


def executeGA(info, population_size, ngen, crossover_rate, mutation_rate, path, method='other', multiprocessing=False):
    now = datetime.now()

    # The new simulation path
    pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - GA'

    # Create the simulation directory
    os.mkdir(pathResult)

    ga = Ga(info.fitness, limInf=info.limitInf, limSup=info.limitSup, populationSize=population_size,
            path=pathResult,
            weights=(-1,), multiprocessing=multiprocessing)

    hof = ga.run(nGenerations=ngen, crossOver=crossover_rate, mutation=mutation_rate,
                 method=method, saveGeneration=1, verbose=False)


def executePSO(info, population_size, ngen, phi1, phi2, path, method='other', multiprocessing=False):
    now = datetime.now()

    # The new simulation path
    pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - PSO'

    # Create the simulation directory
    os.mkdir(pathResult)

    # Set the PSO
    pso = PSO(info.fitness, limInf=info.limitInf, limSup=info.limitSup, populationSize=population_size, path=pathResult,
              weights=(-1,), phi1=phi1, phi2=phi2, multiprocessing=multiprocessing)

    # Run the PSO
    best = pso.run(nGenerations=ngen, saveEpoch=1, verbose=False, method=method)


def validateNSGA2():
    # 3 objetivos, 12 variaveis
    problem = DTLZ2(3, 12)

    # Initialize the algorithm
    algorithm = NSGAII(problem, population_size=100)

    # Run the algorithm
    algorithm.run(100 * 100)

    # Get the pareto front
    results = nondominated(algorithm.result)

    # display the results
    fig = plt.figure()

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.scatter([s.objectives[0] for s in results],
               [s.objectives[1] for s in results],
               [s.objectives[2] for s in results])
    ax.set_title("NSGA-II - Visão Superior")
    ax.set_xlim([0, 1.1])
    ax.set_ylim([0, 1.1])
    ax.set_zlim([0, 1.1])
    ax.view_init(elev=30.0, azim=45.0)
    ax.locator_params(nbins=4)

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ind = ax.scatter([s.objectives[0] for s in results],
                     [s.objectives[1] for s in results],
                     [s.objectives[2] for s in results])
    ax.set_title("NSGA-II - Visão Lateral")
    ax.set_xlim([0, 1.1])
    ax.set_ylim([0, 1.1])
    ax.set_zlim([0, 1.1])
    ax.view_init(elev=0, azim=140)
    ax.locator_params(nbins=4)
    fig.legend([ind], ["Indivíduo"], ncol=1, loc=8)
    plt.show()


def validateMOPSO():
    # 3 objetivos, 12 variaveis
    problem = DTLZ2(3, 12)

    # Initialize the algorithm
    algorithm = OMOPSO(problem, swarm_size=100, epsilons=[0.05])

    # Run the algorithm
    algorithm.run(100 * 100)

    # Get the pareto front
    results = nondominated(algorithm.result)

    # display the results
    fig = plt.figure()

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.scatter([s.objectives[0] for s in results],
               [s.objectives[1] for s in results],
               [s.objectives[2] for s in results])
    ax.set_title("OMOPSO - Visão Superior")
    ax.set_xlim([0, 1.1])
    ax.set_ylim([0, 1.1])
    ax.set_zlim([0, 1.1])
    ax.view_init(elev=30.0, azim=45.0)
    ax.locator_params(nbins=4)

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    part = ax.scatter([s.objectives[0] for s in results],
                      [s.objectives[1] for s in results],
                      [s.objectives[2] for s in results])
    ax.set_title("OMOPSO - Visão Lateral")
    ax.set_xlim([0, 1.1])
    ax.set_ylim([0, 1.1])
    ax.set_zlim([0, 1.1])
    ax.view_init(elev=0, azim=140)
    ax.locator_params(nbins=4)
    fig.legend([part], ["Partícula"], ncol=1, loc=8)
    plt.show()


def compareMO():
    # 3 objetivos, 12 variaveis
    problem = DTLZ2(3, 12)

    algorithms = [
        NSGAII,
        (OMOPSO, {"epsilons": [0.05]})
    ]

    # Run the experimenter
    results = experiment(algorithms, problem, seeds=1, nfe=100 * 100)

    # display the results
    fig = plt.figure()

    for i, algorithm in enumerate(six.iterkeys(results)):
        result = results[algorithm]["DTLZ2"][0]

        ax = fig.add_subplot(2, 2, i + 1, projection='3d')
        ax.scatter([s.objectives[0] for s in result],
                   [s.objectives[1] for s in result],
                   [s.objectives[2] for s in result])
        ax.set_title(algorithm + " - Visão Superior")
        ax.set_xlim([0, 1.1])
        ax.set_ylim([0, 1.1])
        ax.set_zlim([0, 1.1])
        ax.view_init(elev=30.0, azim=45.0)
        ax.locator_params(nbins=4)

    for i, algorithm in enumerate(six.iterkeys(results)):
        result = results[algorithm]["DTLZ2"][0]

        ax = fig.add_subplot(2, 2, i + 3, projection='3d')
        ax.scatter([s.objectives[0] for s in result],
                   [s.objectives[1] for s in result],
                   [s.objectives[2] for s in result])
        ax.set_title(algorithm + " - Visão Lateral")
        ax.set_xlim([0, 1.1])
        ax.set_ylim([0, 1.1])
        ax.set_zlim([0, 1.1])
        ax.view_init(elev=0, azim=140)
        ax.locator_params(nbins=4)

    plt.show()

def summary(path, plot=False):
    results = []
    sucessoGa = []
    sucessoPSO = []

    for i in os.listdir(path):

        infoFile = open(path + '/' + i + '/info.pkl', 'rb')
        file = open(path + '/' + i + '/historic.pkl', 'rb')

        method = i.split('-')[-1].strip()

        data = pickle.load(file)
        info = pickle.load(infoFile)

        dataPrint = [
            method,
            str(round(info['bestVal'][0], 3)),
            str(list(map(lambda x: round(x, 3), info['bestInd']))),
            str(info['elapsedTime']),
            info['evalTotal'],
            info['nGeneration'],
            'Sim' if round(info['bestVal'][0], 3) < 0.01 else 'Não'
        ]

        if method == 'GA':
            if round(info['bestVal'][0], 3) < 0.01:
                sucessoGa.append(1)
            else:
                sucessoGa.append(0)

        if method == 'PSO':
            if round(info['bestVal'][0], 3) < 0.01:
                sucessoPSO.append(1)
            else:
                sucessoPSO.append(0)

        results.append(dataPrint)

    print(tabulate(results, headers=[
        'Method',
        'Best Val',
        'Best Ind',
        'Elapsed Time',
        'Total Eval',
        'Total Gen',
        'Sucesso'
    ], tablefmt="pipe"))

    print('Sucesso GA: ' + str(sum(sucessoGa) / len(sucessoGa) * 100) + '%')
    print('Sucesso PSO: ' + str(sum(sucessoPSO) / len(sucessoPSO) * 100) + '%')

    if plot:

        for i in os.listdir(path):
            method = i.split('-')[-1].strip()

            statsFile = open(path + '/' + i + '/statsPopulation.pkl', 'rb')

            stats = pickle.load(statsFile)

            geracao = list(range(1, len(stats) + 1))

            plotHistoric(geracao, raio=stats, together=False, title=method)


def plotProgress(fitness, path):
    file = open(path + '/paretoFrontier.pkl', 'rb')

    data = pickle.load(file)

    results = []

    for i in data:
        res = {'pop': i, 'best': False}
        results.append(res)

    print(list(data[0]))
    print(fitness(data[0]))

    title = ['Geração ' + str(2 * i) for i in range(0, len(results))]

    # Plot the evolution
    plotFunctionCountour(results, title=title, legend='Indivíduo')


if __name__ == '__main__':

    now = datetime.now()

    a = Benchmark(nvar=2).getBohachevskyInfo()
    b = Benchmark(nvar=2).getSchafferInfo()
    c = Benchmark(nvar=2).getSchwefelInfo()

    infos = [a, b, c]

    for info in infos:

        # The new simulation path
        path = 'Validation/Results/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - Experiment ' + info.name

        # Create the simulation directory
        os.mkdir(path)

        print("\n\n---------------------------------- " + info.name + "-----------------------------------------------")

        for i in range(10):
            # Execute the ga
            # executeGA(info, population_size=1000, ngen=1000, crossover_rate=0.8, mutation_rate=0.2, path=path)
            executeGA(info, population_size=1000, ngen=300, crossover_rate=0.8, mutation_rate=0.2, path=path,
                      method='modified', multiprocessing=False)

            # Execute the pso
            # executePSO(info, population_size=1000, ngen=1000, phi1=0.5, phi2=0.5, path=path)
            executePSO(info, population_size=1000, ngen=300, phi1=0.5, phi2=0.5, path=path,
                       method='modified', multiprocessing=False)

        # Show the summary for the experiment
        summary(path, plot=False)
