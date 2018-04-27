#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing as mp
from scipy.stats.distributions import norm
from datetime import datetime
from deap import creator, base, tools, algorithms
from pyDOE import *
import copy
from helpers import *


class Ga:

    def __init__(self, fitnessFunction, limInf, limSup, path, weights, x0=[], populationSize=100, multiprocessing=False):

        self.fitnessFunction = fitnessFunction
        self.limInf = limInf
        self.limSup = limSup
        self.x0 = x0
        self.populationSize = populationSize
        self.path = path
        self.weights = weights
        self.multiprocessing = multiprocessing

        creator.create("FitnessMulti", base.Fitness, weights=self.weights)
        creator.create("Individual", list, fitness=creator.FitnessMulti)

        self.toolbox = base.Toolbox()

        if self.multiprocessing:

            self.pool = mp.Pool()
            self.toolbox.register("map", self.pool.map)

    def _generate(self, size):

        # Create the individual
        return creator.Individual([np.random.uniform(self.limInf[x], self.limSup[x]) for x in range(size)])

    # def _initPopulationLHS(self, pcls, ind_init):
    #
    #     design = lhs(len(self.x0), samples=self.populationSize, criterion='center')
    #
    #     means = self.x0
    #
    #     for i in range(0, len(self.x0)):
    #         design[:, i] = list(map(lambda x: round(x, 4), norm(loc=means[i], scale=self.stdv[i]).ppf(design[:, i])))
    #
    #     design = np.concatenate((design, [self.x0]))
    #
    #     return pcls(ind_init(list(c)) for c in design)


    def run(self, method='modified', nGenerations=10, crossOver=0.5, mutation=0.1, initPop=None, saveGeneration=20, verbose=True):

        # Start time
        start_time = datetime.now()

        self.toolbox.register("individual", self._generate, size=len(self.limSup))
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self.fitnessFunction)

        # Initialize the population
        pop = self.toolbox.population(n=self.populationSize)

        if len(self.weights) > 2:
            self.toolbox.register("mate", tools.cxSimulatedBinaryBounded, eta=5, low=self.limInf, up=self.limSup)
            self.toolbox.register("mutate", tools.mutPolynomialBounded, eta=5, low=self.limInf, up=self.limSup,
                                  indpb=0.05)
            self.toolbox.register("select", tools.selNSGA2)
            hof = tools.ParetoFront()
        else:
            self.toolbox.register("mate", tools.cxOnePoint)
            self.toolbox.register("mutate", tools.mutUniformInt, low=self.limInf, up=self.limSup, indpb=0.01)
            self.toolbox.register("select", tools.selTournament, tournsize=2)
            hof = tools.HallOfFame(3)

        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("min", np.min, axis=0)
        stats.register("mean", np.mean, axis=0)
        stats.register("max", np.max, axis=0)

        #population, logbook, statsPopulation = self._eaSimple(pop, self.toolbox, cxpb=crossOver, mutpb=mutation, stats=stats, ngen=nGenerations, halloffame=hof,
        #                saveGeneration=saveGeneration, verbose=verbose)

        population, logbook, statsPopulation, evaluationTotal = self._eaMuPlusLambdaModified(pop, self.toolbox, mu=self.populationSize, lambda_=self.populationSize,
                             cxpb=crossOver,
                             mutpb=mutation, nGeneration=nGenerations, method=method, halloffame=hof, stats=stats,
                             saveGeneration=saveGeneration, verbose=verbose)

        # End time
        end_time = datetime.now()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(self.path + '/parameters.txt', "w")
        file.write("Method: GA - " + method + "\n")

        if method == 'modified':
            file.write("N Stall: " + str(nGenerations) + "\n")
        else:
            file.write("N Generation: " + str(nGenerations) + "\n")

        file.write("N Population: " + str(self.populationSize) + "\n")
        file.write("Population Init: " + str(initPop) + "\n")
        file.write("Crossover Rate: " + str(crossOver) + "\n")
        file.write("Mutation Rate: " + str(mutation) + "\n")
        file.write("Limit Max: " + str(self.limSup) + "\n")
        file.write("Limit Min: " + str(self.limInf) + "\n")
        file.write("Weights: " + str(self.weights) + "\n")
        file.write("Elapsed Time: " + str(elapsed_time) + "\n")
        file.close()

        # Save generation
        output = open(self.path + "/statsPopulation.pkl", 'wb')
        pickle.dump(statsPopulation, output)
        output.close()

        # Select the historic
        genMin, genMean, genMax = logbook.select("min", 'mean', 'max')

        # Save historic
        output = open(self.path + "/historic.pkl", 'wb')
        pickle.dump([genMin, genMean, genMax], output)
        output.close()

        info = {
            'bestInd': list(hof[0]),
            'bestVal': self.fitnessFunction(hof[0]),
            'evalTotal': evaluationTotal,
            'elapsedTime': str(elapsed_time),
            'nGeneration': len(genMin)
        }

        # Save additional information
        output = open(self.path + "/info.pkl", 'wb')
        pickle.dump(info, output)
        output.close()

        # Plot the historic
        plotHistoric(range(0, len(statsPopulation)), statsPopulation, savePath=self.path)

        return hof


    def _eaMuPlusLambdaModified(self, population, toolbox, mu, lambda_, cxpb, mutpb, nGeneration,
                                method='modified',
                                stats=None, halloffame=None, saveGeneration=0, verbose=True):

        gen = 1
        saveGen = saveGeneration
        counter = 0
        minRaio = 0
        statsPopulation = []
        evaluationTotal = 0

        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]

        popFit = map(lambda x: list(x), invalid_ind)

        fitnesses = toolbox.map(toolbox.evaluate, popFit)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        evaluationTotal = evaluationTotal + len(invalid_ind)

        if halloffame is not None:
            halloffame.update(population)

        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)

        if verbose:
            print(logbook.stream)

        # Begin the generational process
        while (method == 'modified' and gen < 5000) or (gen < nGeneration + 2 and method != 'modified'):

            # Vary the population
            offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]

            popFit = map(lambda x: list(x), invalid_ind)

            fitnesses = list(map(toolbox.evaluate, popFit))
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            evaluationTotal = evaluationTotal + len(invalid_ind)

            # Update the hall of fame with the generated individuals
            if halloffame is not None:
                halloffame.update(offspring)

            # Select the next generation population
            population[:] = toolbox.select(population + offspring, mu)

            # Update the statistics with the new population
            record = stats.compile(population) if stats is not None else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)

            if verbose:
                print(logbook.stream)

            if saveGen == saveGeneration:
                statsPopulation.append(halloffame[0].fitness.values[0])
                saveGen = 0

            saveGen = saveGen + 1

            gen = gen + 1

            if method == 'modified':

                # Verifica estagnacao
                if round(record['min'][0], 3) <= minRaio:
                    counter = counter + 1
                else:
                    counter = 0

                # Atualiza as referencias
                minRaio = round(record['min'][0], 3)

                # Se estourou o limite, termina execucao
                if counter >= nGeneration:
                    break

        return population, logbook, statsPopulation, evaluationTotal

    def _eaSimple(self, population, toolbox, cxpb, mutpb, ngen, stats=None,
                  halloffame=None, saveGeneration=0, verbose=True):

        statsPopulation = []
        saveGen = saveGeneration

        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        if halloffame is not None:
            halloffame.update(population)

        record = stats.compile(population) if stats else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)

        if verbose:
            print(logbook.stream)

        statsPopulation.append(np.asarray(copy.deepcopy(population)))

        # Begin the generational process
        for gen in range(1, ngen + 1):
            # Select the next generation individuals
            offspring = toolbox.select(population, len(population))

            # Vary the pool of individuals
            offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            # Update the hall of fame with the generated individuals
            if halloffame is not None:
                halloffame.update(offspring)

            # Replace the current population by the offspring
            population[:] = offspring

            if verbose:
                print(logbook.stream)

            if saveGen == saveGeneration:
                statsPopulation.append(np.asarray(copy.deepcopy(population)))
                saveGen = 0

            saveGen = saveGen + 1

            # Append the current generation statistics to the logbook
            record = stats.compile(population) if stats else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)


        return population, logbook, statsPopulation
