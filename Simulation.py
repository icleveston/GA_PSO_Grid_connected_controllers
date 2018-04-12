#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from platypus import *
from helpers import *
from TheProblem import TheProblem
import pickle
import time
import matlab.engine

limSup = [0,    0,    0,   0, 100,  0, 50,  0, 50,  0, 50, 0]
limInf = [-15, -15, -15, -15, 0, -100, 0, -50, 0, -50, 0, -50]


class Simulation:

    def __init__(self):

        self.eng = matlab.engine.start_matlab()

        self.c = self.eng.Controle()

        self.problem = TheProblem(self.eng, self.c)

    def executeNSGA2(self, population_size=10, ngen=10, crossover_rate=0.8, mutation_rate=0.2):

        # The new simulation path
        path = 'Simulation/' + time.strftime("%m-%d-%Y-%H-%M-%S") + ' - NSGA2'

        # Create the simulation directory
        os.mkdir(path)

        # Initialize the algorithm
        algorithm = NSGAII(self.problem, population_size=population_size,
                           variator=GAOperator(SBX(probability=crossover_rate), PM(probability=mutation_rate)))

        # Start time
        start_time = time.time()

        # Run the algorithm
        algorithm.run(ngen * population_size)

        # Get the pareto front
        results = nondominated(algorithm.result)

        # End time
        end_time = time.time()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(path + '/parameters.txt', "w")
        file.write("Method: NSGA-II\n")
        file.write("N Population: " + str(population_size) + "\n")
        file.write("N Generation: " + str(ngen) + "\n")
        file.write("Crossover Rate: " + str(crossover_rate) + "\n")
        file.write("Mutation Rate: " + str(mutation_rate) + "\n")
        file.write("Elapsed Time: " + str(time.strftime("%H:%M:%S", time.gmtime(elapsed_time))) + "\n")
        file.close()

        pareto = []

        for ind in results:
            res = {'ind': list(ind.variables), 'val': list(ind.objectives)}
            pareto.append(res)

        # Save the pareto front
        output = open(path + '/pareto.pkl', 'wb')
        pickle.dump(pareto, output)
        output.close()

        # Save the pareto front
        output = open(path + '/historic.pkl', 'wb')
        pickle.dump(algorithm.historic, output)
        output.close()

        return path


    def executeMOPSO(self, swarm_size=10, ngen=10, crossover_rate=0.8, mutation_rate=0.2):

        # The new simulation path
        path = 'Simulation/' + time.strftime("%m-%d-%Y-%H-%M-%S") + ' - MOPSO'

        # Create the simulation directory
        os.mkdir(path)

        # Initialize the algorithm
        algorithm = OMOPSO(self.problem, swarm_size=100, epsilons=[0.05])

        # Start time
        start_time = time.time()

        # Run the algorithm
        algorithm.run(ngen * swarm_size)

        # Get the pareto front
        results = nondominated(algorithm.result)

        # End time
        end_time = time.time()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(path + '/parameters.txt', "w")
        file.write("Method: MOPSO\n")
        file.write("N Swarm: " + str(swarm_size) + "\n")
        file.write("N Generation: " + str(ngen) + "\n")
        file.write("Crossover Rate: " + str(crossover_rate) + "\n")
        file.write("Mutation Rate: " + str(mutation_rate) + "\n")
        file.write("Elapsed Time: " + str(time.strftime("%H:%M:%S", time.gmtime(elapsed_time))) + "\n")
        file.close()

        pareto = []

        for ind in results:
            res = {'ind': list(ind.variables), 'val': list(ind.objectives)}
            pareto.append(res)

        # Save the pareto front in Solution Class
        output = open(path + '/pareto.pkl', 'wb')
        pickle.dump(results, output)
        output.close()

        # Save the pareto front
        output = open(path + '/historic.pkl', 'wb')
        pickle.dump(algorithm.historic, output)
        output.close()

        return path


    def executeNSGA3(self):
        pass
        # algorithm = NSGAIII(problem, population_size=100, divisions_outer=12, variator=GAOperator(SBX(probability=0.8), PM(probability=0.2)))


    def getReport(self, path):

        self.plotHistoric(path)

        file = open(path + '/pareto.pkl', 'rb')

        data = pickle.load(file)

        # Normalize data
        a = normalizeResults(data, attenuation=False)

        if len(a) > 0:

            for i in a[0:3]:

                print(i['ind'])
                print(i['val'])

                plotConversor(i['ind'], self.eng, self.c)


    def plotHistoric(self, path):

        file = open(path + '/historic.pkl', 'rb')

        data = pickle.load(file)

        raio = []
        bode = []
        ise = []

        for x in data:
            raio.append(x[0])
            bode.append(x[1])
            ise.append(x[2])

        geracao = [x for x in range(0, len(raio))]

        # Plot the historic
        plotHistoric(geracao, raio, bode, ise)


    def compare(self, path1, path2):

        file = open(path1 + '/pareto.pkl', 'rb')

        data = pickle.load(file)

        data_solution = []

        problem = TheProblem(self.eng, self.c)

        for e in data:
            s = Solution(problem, e['ind'], e['val'])
            data_solution.append(s)


        # Create the hypervolume
        hyp = Hypervolume(minimum=[0, 0, 0], maximum=[2, 2, 2])

        hyp_result = hyp.calculate(data_solution)

        print(hyp_result)

