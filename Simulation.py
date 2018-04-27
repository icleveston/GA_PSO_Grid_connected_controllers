#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from platypus import *
from helpers import *
from nuvem import *
from ga import *
from pso import *
import pickle
import time
from math import sqrt, exp
import matlab.engine

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

limSup = [0, 0, 0, 0, 100, 0, 50, 0, 50, 0, 50, 0]
limInf = [-15, -15, -15, -15, 0, -100, 0, -50, 0, -50, 0, -50]

eng = matlab.engine.start_matlab()

c = eng.Controle()


def evaluateRaio(individual):
    MF1full2 = A_til_1 + B_til * individual
    MF2full2 = A_til_2 + B_til * individual

    ai = np.concatenate((MF1full2, MF2full2), axis=1)

    return nuvem(ai)[0],


def evaluate(individual):
    a = matlab.double(individual)

    res = eng.testa(c, a, nargout=1)

    x = res[0][0]
    y = res[0][1]

    # if res[0][0] > 0.997:
    #     return 10 ** res[0][0],
    # else:
    #     return res[0][1]*res[0][0],

    # return y*x + (1 / (1 + exp(-1000*x + 1000)))*(20*(x+10)**2 - 20*120),
    # return y*x + (1 / (1 + exp(-5000*x + 5000)))*(450*x - 440),
    # return y*x*(-(1 / (1 + exp(-5000*x + 5000))) +1) + (1 / (1 + exp(-5000*x + 5000)))*10*x,
    return y*x + (1 / (1 + exp(-1000*x + 1000)))*(380*x - 360),


def evaluateAll(individual):
    a = matlab.double(individual)

    res = eng.testa(c, a, nargout=1)

    return res[0][0], res[0][1] #, res[0][2]


def executeGA(info, population_size, ngen, crossover_rate, mutation_rate, path, method='other', multiprocessing=False):
    now = datetime.now()

    # The new simulation path
    pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - GA'

    # Create the simulation directory
    os.mkdir(pathResult)

    ga = Ga(info['fitness'], limInf=info['limInf'], limSup=info['limSup'], populationSize=population_size,
            path=pathResult,
            weights=(-1,), multiprocessing=multiprocessing)

    hof = ga.run(nGenerations=ngen, crossOver=crossover_rate, mutation=mutation_rate,
                 method=method, saveGeneration=1, verbose=True)

    return hof[0]


def executePSO(info, population_size, ngen, phi1, phi2, path, method='other', multiprocessing=False):
    now = datetime.now()

    # The new simulation path
    pathResult = path + '/' + now.strftime("%m-%d-%Y-%H-%M-%S.%f") + ' - PSO'

    # Create the simulation directory
    os.mkdir(pathResult)

    # Set the PSO
    pso = PSO(info['fitness'], limInf=info['limInf'], limSup=info['limSup'], populationSize=population_size,
              path=pathResult,
              weights=(-1,), phi1=phi1, phi2=phi2, multiprocessing=multiprocessing)

    # Run the PSO
    best = pso.run(nGenerations=ngen, saveEpoch=1, verbose=True, method=method)

    return best


def executeNSGA2(self, population_size=100, ngen=100, crossover_rate=0.8, mutation_rate=0.2):
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


def executeMOPSO(self, swarm_size=500, ngen=300):
    # The new simulation path
    path = 'Simulation/' + time.strftime("%m-%d-%Y-%H-%M-%S") + ' - MOPSO'

    # Create the simulation directory
    os.mkdir(path)

    # Initialize the algorithm
    algorithm = OMOPSO(self.problem, swarm_size=swarm_size, epsilons=[0.05])

    # Start time
    start_time = time.time()

    # Run the algorithm
    algorithm.run(ngen * swarm_size)

    # Get the pareto front
    results = algorithm.result

    # End time
    end_time = time.time()

    elapsed_time = end_time - start_time

    # Save parameters
    file = open(path + '/parameters.txt', "w")
    file.write("Method: MOPSO\n")
    file.write("N Swarm: " + str(swarm_size) + "\n")
    file.write("N Generation: " + str(ngen) + "\n")
    # file.write("Phi1: " + str(crossover_rate) + "\n")
    # file.write("Phi2: " + str(mutation_rate) + "\n")
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


def getReport(self, path):
    # self.plotHistoric(path)

    file = open(path + '/pareto.pkl', 'rb')

    data = pickle.load(file)

    print(len(data))

    # Normalize data
    a = normalizeResults(data, attenuation=False)

    if len(a) > 0:

        for i in a[0:3]:
            print(i['ind'])
            print(i['val'])

            plotConversor(i['ind'], self.eng, self.c)


def _plotHistoric(path):
    file = open(path + '/historic.pkl', 'rb')

    data = pickle.load(file)

    raio = []
    bode = []
    ise = []

    for x in data:
        raio.append(x[0])

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


if __name__ == '__main__':

    info = {'fitness': evaluate, 'limInf': limInf, 'limSup': limSup}

    path = 'Simulation/Results/' + time.strftime("%d-%m-%Y-%H-%M-%S") + ' - Experiment'

    # Create the simulation directory
    os.mkdir(path)

    for i in range(0, 10):

         best = executePSO(info, population_size=50, ngen=200, phi1=0.5, phi2=0.5, path=path, method='other', multiprocessing=True)

         print(str(best))
         print(evaluateAll(best))

    # for i in range(0, 1):
    #     hof = executeGA(info, population_size=50, ngen=10, crossover_rate=0.8, mutation_rate=0.2, path=path,
    #                     method='other', multiprocessing=True)
    #
    #     print(str(hof))
    #     # print(str(best).replace(',', ' '))
    #     print(evaluateAll(hof))
