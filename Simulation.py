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

#limSup = [0, 0, 0, 0, 100, 0, 50, 0, 50, 0, 50, 0]
#limInf = [-15, -15, -15, -15, 0, -100, 0, -50, 0, -50, 0, -50]
limSup = [0, 0, 0, 0.5, 100, 0, 50, 0, 50, 0, 50, 0]
limInf = [-15, 1e-10, -15, -0.5, 0, -100, 0, -50, 0, -50, 0, -50]
#limSup = [-15, 0.0, -2, -1.292859528039262, 78, -76.00080872558523, 34.336290645512605, -33.407788915312466, 25.6114464208353,
# #  -26, 9.352959226656882, -8.283920091847143]
# limSup = [-14, 0, -1, -1, 85, -60, 38, -28, 30, -20, 15, -4]
# limInf = [-16.6, 0, -2.1, -1.5, 70, -85, 30, -38, 20, -31, 5, -20]


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

    return y * x + (1 / (1 + exp(-1000 * x + 1000))) * (380 * x - 360),


def evaluateAll(individual):
    a = matlab.double(individual)

    res = eng.testa(c, a, nargout=1)

    return res[0][0], res[0][1]  # , res[0][2]


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

def fitness(x, y):

    return y*x + (1 / (1 + exp(-1000 * x + 1000))) * (380 * x - 360)


if __name__ == '__main__':

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # x = np.arange(0.9, 1.5, 0.005)
    # y = np.arange(0, 5, 0.005)
    # X, Y = np.meshgrid(x, y)
    # zs = np.array([fitness(x, y) for x, y in zip(np.ravel(X), np.ravel(Y))])
    # Z = zs.reshape(X.shape)
    #
    # CS = ax.plot_surface(X, Y, Z,  cmap=cm.coolwarm)
    #
    # ax.set_xlabel("Raio")
    # ax.set_ylabel("Bode")
    # ax.set_zlabel("Custo")
    # #plt.title(title)
    #
    # clb = plt.colorbar(CS, orientation='horizontal', shrink=0.5)
    # clb.ax.set_title('Custo')
    #
    # plt.show()
    #
    #
    # exit()

    # path = 'Simulation/Results/29-05-2018-18-34-40 - Experiment'
    #
    # fitnessArray = []
    #
    # for i in os.listdir(path):
    #     infoFile = open(path + "/" + i + "/info.pkl", 'rb')
    #     statsPopulationFile = open(path + "/" + i + "/statsPopulation.pkl", 'rb')
    #     historicFile = open(path + "/" + i + "/historic.pkl", 'rb')
    #     # distanceFile = open(path + "/" + i + "/distance.pkl", 'rb')
    #
    #     info = pickle.load(infoFile)
    #     results = pickle.load(statsPopulationFile)
    #     fitnessArray.append(results)
    #     historic = pickle.load(historicFile)
    #     # distance = pickle.load(distanceFile)
    #     #
    #     # distanceModified = []
    #     #
    #     # print(len(distance))
    #     #
    #     # c = 0
    #     # for d in distance:
    #     #
    #     #     if c == 10:
    #     #         distanceModified.append(d)
    #     #         c = 0
    #     #
    #     #     c = c + 1
    #
    #
    #     # # Create the boxplot
    #     # plt.boxplot(distanceModified)
    #     # plt.title('Distância do enxame ao Gbest')
    #     # plt.grid()
    #     # plt.xlabel('Época')
    #     # plt.ylabel('Distância')
    #     # plt.xticks(np.arange(1, 21), [1] + list(range(10, 200, 10)))
    #     # plt.show()
    #     #
    #     #
    #     # Plot the historic
    #     # plotHistoric(range(0, len(results)), results, savePath=path + "/" + i)
    #     # plotHistoric(range(0, len(results)), results)
    #     #
    #     # plt.plot(range(0, len(results)+1), historic[0], label="Min")
    #     # plt.plot(range(0, len(results) + 1), historic[1], label="Mean")
    #     # plt.plot(range(0, len(results) + 1), historic[2], label="Max")
    #     # # plt.plot(range(0, len(results)), results)
    #     # plt.grid()
    #     # plt.xlabel("Geração")
    #     # plt.ylabel("Fitness")
    #     # plt.legend()
    #     # plt.show()
    #
    #
    #     print(info['bestInd'])
    #     print(evaluateAll(info['bestInd']))
    #     print(evaluate(info['bestInd']))
    #
    #
    # exit()

    info = {'fitness': evaluate, 'limInf': limInf, 'limSup': limSup}

    path = 'Simulation/Results/' + time.strftime("%d-%m-%Y-%H-%M-%S") + ' - Experiment'

    # Create the simulation directory
    os.mkdir(path)

    # for i in range(0, 10):
    #     best = executePSO(info, population_size=50, ngen=200, phi1=0.5, phi2=0.5, path=path, method='other',
    #                       multiprocessing=True)
    #
    #     print(str(best))
    #     print(evaluateAll(best))

    for i in range(0, 1):
        hof = executeGA(info, population_size=100, ngen=200, crossover_rate=0.8, mutation_rate=0.2, path=path,
                        method='other', multiprocessing=True)

        print(str(hof))
        print(evaluateAll(hof))

