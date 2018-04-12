#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
from math import sin
import scipy
from helpers import *
import os.path
from ga import *
from pso import *
from platypus import *
from deap.benchmarks.tools import convergence, diversity


class Validation:
    # Search space limits
    limSup = [10, 10]
    limInf = [0, 0]

    def validateGA(self):
        # The new simulation path
        path = 'Validation/' + time.strftime("%m-%d-%Y-%H-%M-%S") + ' - GA'

        # Create the simulation directory
        os.mkdir(path)

        ga = Ga(self._evaluate, limInf=Validation.limInf, limSup=Validation.limSup, populationSize=40, path=path,
                weights=(-1,))

        hof = ga.run(nGenerations=10, crossOver=0.6, mutation=0.2, method='other', saveGeneration=2)

        file = open(path + '/paretoFrontier.pkl', 'rb')

        data = pickle.load(file)

        results = []

        for i in data:
            res = {'pop': i, 'best': False}
            results.append(res)
            print(sum(i))

        print(list(hof[0]))
        print(self._evaluate(hof[0]))

        title = ['Geração ' + str(2 * i) for i in range(0, len(results))]

        # Plot the evolution
        self._plotFitnessEvolutionContour(results, title=title, legend='Indivíduo')

    def validatePSO(self):
        # The new simulation path
        path = 'Validation/' + time.strftime("%m-%d-%Y-%H-%M-%S") + ' - PSO'

        # Create the simulation directory
        os.mkdir(path)

        # Set the PSO
        pso = PSO(self._evaluate, limInf=Validation.limInf, limSup=Validation.limSup, populationSize=20, path=path,
                  weights=(-1,), phi1=0.2, phi2=0.8)

        # Run the PSO
        _, _, best, results = pso.run(nGenerations=25, saveEpoch=5)

        print(best)
        print(self._evaluate(best))

        title = ['Época ' + str(5 * i) for i in range(0, len(results))]

        # Plot the evolution
        self._plotFitnessEvolutionContour(results, title=title, legend='Partícula')

    def validateNSGA2(self):

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
        ax.scatter([s.objectives[0] for s in results],
                   [s.objectives[1] for s in results],
                   [s.objectives[2] for s in results])
        ax.set_title("NSGA-II - Visão Lateral")
        ax.set_xlim([0, 1.1])
        ax.set_ylim([0, 1.1])
        ax.set_zlim([0, 1.1])
        ax.view_init(elev=0, azim=140)
        ax.locator_params(nbins=4)

        plt.show()

        # convergence()
        # diversity()
        # scipy.stats.mannwhitneyu()

    def validateMOPSO(self):

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
        ax.scatter([s.objectives[0] for s in results],
                   [s.objectives[1] for s in results],
                   [s.objectives[2] for s in results])
        ax.set_title("OMOPSO - Visão Lateral")
        ax.set_xlim([0, 1.1])
        ax.set_ylim([0, 1.1])
        ax.set_zlim([0, 1.1])
        ax.view_init(elev=0, azim=140)
        ax.locator_params(nbins=4)

        plt.show()

    def compareMO(self):

        # 3 objetivos, 12 variaveis
        problem = DTLZ2(3, 12)

        algorithms = [
            # SPEA2,
            NSGAII,
            (OMOPSO, {"epsilons": [0.05]})
            # SMPSO
        ]

        # Run the experimenter
        results = experiment(algorithms, problem, seeds=1, nfe=100*100)

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

    def plotFitness3D(self):
        plotFunction3D(self._evaluate, limitMin=0, limitMax=10, step=0.005)

    def plotFitnessContour(self):
        plotFunctionCountour(self._evaluate, limitMin=0, limitMax=10, step=0.005, legend=True)

    def _plotFitnessEvolutionContour(self, all, title, legend):

        plotFunctionCountour(self._evaluate, limitMin=0, limitMax=10, step=0.005, all=all, title=title, legend=legend)

    def _evaluate(self, individual):
        # y = x1*sin(4*x1) + 1.1*x2*sin(2*x2)
        return individual[1] * sin(4 * individual[1]) + 1.1 * individual[0] * sin(2 * individual[0]),
