#!/usr/bin/env python
# -*- coding: utf-8 -*-

from helpers import *
import os.path
from ga import *
from pso import *
from platypus import *
from Validation.Benchmark import Benchmark
from datetime import datetime
from tabulate import tabulate


class Validation:

    def runExperiment(self, executionTimes=10):

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

            for i in range(executionTimes):
                # Execute the ga
                self._executeGA(info, population_size=1000, ngen=1000, crossover_rate=0.8, mutation_rate=0.2, path=path)
                # self._executeGA(info, population_size=100, ngen=300, crossover_rate=0.8, mutation_rate=0.2, path=path,
                #                 method='modified')

                # Execute the pso
                self._executePSO(info, population_size=1000, ngen=1000, phi1=0.5, phi2=0.5, path=path)
                # self._executePSO(info, population_size=100, ngen=300, phi1=0.5, phi2=0.5, path=path,
                #                  method='modified')

            # Show the summary for the experiment
            self.summary(path, plot=False)

    def _executeGA(self, info, population_size, ngen, crossover_rate, mutation_rate, path, method='other'):

        now = datetime.now()

        # The new simulation path
        pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - GA'

        # Create the simulation directory
        os.mkdir(pathResult)

        # Start time
        start_time = datetime.now()

        ga = Ga(info.fitness, limInf=info.limitInf, limSup=info.limitSup, populationSize=population_size,
                path=pathResult,
                weights=(-1,))

        hof, _, logbook, statsPopulation, evaluationTotal = ga.run(nGenerations=ngen, crossOver=crossover_rate,
                                                                   mutation=mutation_rate,
                                                                   method=method, saveGeneration=1, verbose=False)

        # End time
        end_time = datetime.now()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(pathResult + '/parameters.txt', "w")
        file.write("Method: GA\n")
        file.write("Fitness Function: " + info.name + "\n")
        file.write("N Population: " + str(population_size) + "\n")
        file.write("N Generation: " + str(ngen) + "\n")
        file.write("Crossover Rate: " + str(crossover_rate) + "\n")
        file.write("Mutation Rate: " + str(mutation_rate) + "\n")
        file.write("Elapsed Time: " + str(elapsed_time) + "\n")
        file.close()

        # Save generation
        output = open(pathResult + "/statsPopulation.pkl", 'wb')
        pickle.dump(statsPopulation, output)
        output.close()

        # Select the historic
        genMin, genMean, genMax = logbook.select("min", 'mean', 'max')

        # Save historic
        output = open(pathResult + "/historic.pkl", 'wb')
        pickle.dump([genMin, genMean, genMax], output)
        output.close()

        info = {
            'bestInd': list(hof[0]),
            'bestVal': info.fitness(hof[0]),
            'evalTotal': evaluationTotal,
            'elapsedTime': str(elapsed_time),
            'nGeneration': len(genMin)
        }

        # Save additional information
        output = open(pathResult + "/info.pkl", 'wb')
        pickle.dump(info, output)
        output.close()

    def _executePSO(self, info, population_size, ngen, phi1, phi2, path, method='other'):

        now = datetime.now()

        # The new simulation path
        pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - PSO'

        # Create the simulation directory
        os.mkdir(pathResult)

        # Start time
        start_time = datetime.now()

        # Set the PSO
        pso = PSO(info.fitness, limInf=info.limitInf, limSup=info.limitSup, populationSize=population_size, path=path,
                  weights=(-1,), phi1=phi1, phi2=phi2)

        # Run the PSO
        pop, logbook, best, statsPopulation, evaluationTotal = pso.run(nGenerations=ngen, saveEpoch=1, verbose=False, method=method)

        # End time
        end_time = datetime.now()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(pathResult + '/parameters.txt', "w")
        file.write("Method: PSO\n")
        file.write("Fitness Function: " + info.name + "\n")
        file.write("N Population: " + str(population_size) + "\n")
        file.write("N Generation: " + str(ngen) + "\n")
        file.write("phi1: " + str(phi1) + "\n")
        file.write("phi2: " + str(phi2) + "\n")
        file.write("Elapsed Time: " + str(elapsed_time) + "\n")
        file.close()

        # Save generation
        output = open(pathResult + "/statsPopulation.pkl", 'wb')
        pickle.dump(statsPopulation, output)
        output.close()

        # Select the historic
        genMin, genMean, genMax = logbook.select("min", 'mean', 'max')

        # Save historic
        output = open(pathResult + "/historic.pkl", 'wb')
        pickle.dump([genMin, genMean, genMax], output)
        output.close()

        info = {
            'bestInd': list(best),
            'bestVal': info.fitness(best),
            'evalTotal': evaluationTotal,
            'elapsedTime': str(elapsed_time),
            'nGeneration': len(genMin)
        }

        # Save additional information
        output = open(pathResult + "/info.pkl", 'wb')
        pickle.dump(info, output)
        output.close()

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

    def plotFitness3D(self):
        plotFunction3D(self.fitness, limitMin=0, limitMax=10, step=0.005)

    def plotFitnessContour(self):
        plotFunctionCountour(self.fitness, limitMin=0, limitMax=10, step=0.005, legend=True)

    def _plotFitnessEvolutionContour(self, all, title, legend):

        plotFunctionCountour(self.fitness, limitMin=0, limitMax=10, step=0.005, all=all, title=title, legend=legend)

    def summary(self, path, plot=False):

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
                str(list(map(lambda x:round(x, 3), info['bestInd']))),
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

        print('Sucesso GA: ' + str(sum(sucessoGa)/len(sucessoGa)*100) + '%')
        print('Sucesso PSO: ' + str(sum(sucessoPSO)/len(sucessoPSO)*100) + '%')

        if plot:

            for i in os.listdir(path):
                method = i.split('-')[-1].strip()

                statsFile = open(path + '/' + i + '/statsPopulation.pkl', 'rb')

                stats = pickle.load(statsFile)

                geracao = list(range(1, len(data[1])))

                plotHistoric(geracao, raio=stats, together=False, title=method)

    def plotProgress(self, path):

        file = open(path + '/paretoFrontier.pkl', 'rb')

        data = pickle.load(file)

        results = []

        for i in data:
            res = {'pop': i, 'best': False}
            results.append(res)

        print(list(data[0]))
        print(self.fitness(data[0]))

        title = ['Geração ' + str(2 * i) for i in range(0, len(results))]

        # Plot the evolution
        self._plotFitnessEvolutionContour(results, title=title, legend='Indivíduo')
