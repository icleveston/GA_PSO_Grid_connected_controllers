#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
from scipy.stats.distributions import norm
from deap import creator, base, tools, algorithms
from pyDOE import *
import copy
import pickle
import numpy as np

np.set_printoptions(threshold=np.nan)


class Ga:

    def __init__(self, fitnessFunction, limInf, limSup, path, weights, x0=[], populationSize=100):

        self.fitnessFunction = fitnessFunction
        self.limInf = limInf
        self.limSup = limSup
        self.x0 = x0
        self.populationSize = populationSize
        self.path = path
        self.weights = weights

        creator.create("FitnessMulti", base.Fitness, weights=self.weights)
        creator.create("Individual", list, fitness=creator.FitnessMulti)

        self.toolbox = base.Toolbox()
        # self.pool = multiprocessing.Pool()
        # self.toolbox.register("map", self.pool.map)

    def _initIndividual(self, icls, content):

        return icls(list(map(lambda x: round(x, 4), content)))

    def _initPopulation(self, pcls, ind_init):

        design = []

        for j in range(0, self.populationSize):

            gene = []

            for i in range(0, len(self.limInf)):
                gene.append(round(np.random.uniform(self.limInf[i], self.limSup[i]), 5))

            design.append(gene)

        return pcls(ind_init(c) for c in design)

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

    # def _initPopulationPseudoRandom(self, pcls, ind_init):
    #
    #     population = []
    #
    #     for i in range(0, round(self.populationSize)):
    #
    #         k1 = round(np.random.uniform(-20, 0), 5)
    #         k2 = round(np.random.uniform(-1, 0.5)*abs(k1), 5)
    #         k3 = round(np.random.uniform(-1, 0.5)*abs(k1), 5)
    #         k4 = round(np.random.uniform(-1, 0.5)*abs(k1), 5)
    #
    #         k5 = round(np.random.uniform(0, 100), 5)
    #         k6 = round(np.random.uniform(-1, -0.5)*k5, 5)
    #
    #         k7 = round(np.random.uniform(abs(max(k5, k6))/3, abs(max(k5, k6))), 5)
    #         k8 = round(np.random.uniform(-1, -0.5)*k7, 5)
    #
    #         k9 = round(np.random.uniform(abs(max(k8, k7))/3, abs(max(k8, k7))), 5)
    #         k10 = round(np.random.uniform(-1, -0.5)*k9, 5)
    #
    #         k11 = round(np.random.uniform(abs(max(k10, k9))/3, abs(max(k10, k9))), 5)
    #         k12 = round(np.random.uniform(-1, -0.5)*k11, 5)
    #
    #         population.append([k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12])
    #
    #     population.append(self.x0)
    #
    #     return pcls(ind_init(c) for c in population)

    def run(self, method='modified', nGenerations=10, crossOver=0.5, mutation=0.1, initPop=None, saveGeneration=20):

        self.toolbox.register("individual_guess", self._initIndividual, creator.Individual)
        self.toolbox.register("population_guess", self._initPopulation, list, self.toolbox.individual_guess)

        # Initialize the population
        if initPop is None:
            pop = self.toolbox.population_guess()
        else:
            pop = [self.toolbox.individual_guess(c) for c in initPop]

        self.toolbox.register("evaluate", self.fitnessFunction)

        if len(self.weights) > 2:
            self.toolbox.register("mate", tools.cxSimulatedBinaryBounded, eta=5, low=self.limInf, up=self.limSup)
            self.toolbox.register("mutate", tools.mutPolynomialBounded, eta=5, low=self.limInf, up=self.limSup,
                                  indpb=0.05)
            self.toolbox.register("select", tools.selNSGA2)
            hof = tools.ParetoFront()
            print("Multi Objetive")
        else:
            self.toolbox.register("mate", tools.cxOnePoint)
            self.toolbox.register("mutate", tools.mutUniformInt, low=0, up=10, indpb=0.01)
            self.toolbox.register("select", tools.selTournament, tournsize=2)
            hof = tools.HallOfFame(3)
            print("Mono Objetive")

        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("min", np.min, axis=0)

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
        file.close()

        self._eaSimple(pop, self.toolbox, cxpb=crossOver, mutpb=mutation, ngen=nGenerations, halloffame=hof,
                       saveGeneration=saveGeneration, path=self.path)

        # self._eaMuPlusLambdaModified(pop, self.toolbox, mu=round(self.populationSize), lambda_=10,
        #                       cxpb=crossOver,
        #                       mutpb=mutation, nGeneration=nGenerations, method=method, halloffame=hof, stats=stats,
        #                       saveGeneration=saveGeneration, verbose=True, path=self.path)

        return hof

    def _eaMuPlusLambdaModified(self, population, toolbox, mu, lambda_, cxpb, mutpb, nGeneration, path,
                                method='modified',
                                stats=None, halloffame=None, saveGeneration=0, verbose=True):

        gen = 1
        saveGen = saveGeneration
        counter = 0
        minRaio = 0
        minBode = 0
        minIse = 0
        minDerivada = 0
        statsPopulation = []

        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        if halloffame is not None:
            halloffame.update(population)

        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)

        if verbose:
            print(logbook.stream)

        # Begin the generational process
        while (method == 'modified') or (gen < nGeneration and method != 'modified'):

            # Vary the population
            offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

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
                statsPopulation.append(np.asarray(copy.deepcopy(population)))
                saveGen = 0
                print("Salvando Geracao")
                print(population)

            saveGen = saveGen + 1

            gen = gen + 1

            if method == 'modified':

                # Verifica estagnacao
                if round(record['min'][0], 4) == minRaio \
                        and round(record['min'][1], 4) == minBode \
                        and round(record['min'][2] / 1e4, 3) == minIse:
                    counter = counter + 1
                else:
                    counter = 0

                # Atualiza as referencias
                minRaio = round(record['min'][0], 4)
                minBode = round(record['min'][1], 4)
                minIse = round(record['min'][2] / 1e4, 3)
                # minDerivada = round(record['min'][3], 4)

                # Se estourou o limite, termina execucao
                if counter >= nGeneration:
                    break

        # Select the historic
        genGen, genMin = logbook.select("gen", "min")

        # Save historic
        output = open(path + "/historic.pkl", 'wb')
        pickle.dump(genMin, output)
        output.close()

        if saveGeneration:
            # Save generation
            output = open(path + "/paretoFrontier.pkl", 'wb')
            pickle.dump(statsPopulation, output)
            output.close()

        return population, logbook

    def _eaSimple(self, population, toolbox, cxpb, mutpb, ngen, stats=None,
                  halloffame=None, path='', saveGeneration=0, ):

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

            if saveGen == saveGeneration:
                statsPopulation.append(np.asarray(copy.deepcopy(population)))
                saveGen = 0
                print("Salvando Geracao")

            saveGen = saveGen + 1

            # Append the current generation statistics to the logbook
            record = stats.compile(population) if stats else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)

        if saveGeneration:
            # Save generation
            output = open(path + "/paretoFrontier.pkl", 'wb')
            pickle.dump(statsPopulation, output)
            output.close()

        return population, logbook
